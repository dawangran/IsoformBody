#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pipeline_tx_qc.py

RNA direct-seq transcript-level QC pipeline with dual metrics:
  presence (bin touch 0/1) and depth (avg base coverage per read).

Pipeline:
  1) Alignment (minimap2 map-ont, --secondary=no) -> OUTDIR/aln.tx.bam (+ .bai)
  2) Forward-only BAM (drop FLAG 0x10) -> OUTDIR/aln.tx.forward.bam (+ .bai)
  3) Counts & exclusion breakdown on forward-only BAM
  4) QC visualizations on chosen BAM (ALL or forward):
       - Classification bar
       - Gene body curves (choose metric: presence/depth), both read-weighted & transcript-weighted
       - Transcript-level heatmap (presence/depth), sorted by completeness
         + right-side capped n_reads strip
         + compact vertical colorbars at far right
       - Read-level heatmap (optional; presence/depth)
       - Venn (detected vs DB)
       - Tables: per-read / per-transcript / heatmap matrix

Robustness knobs:
  --min-tx-len           : filter very short transcripts (default 200 bp)
  --bin-cover-min-frac   : minimal fraction overlap to mark a bin covered in presence mode
  --no-polya-relax       : disable 3' polyA soft-clip relaxation (default enabled)

Deps: numpy, pandas, matplotlib, pysam, tqdm
Optional: matplotlib-venn (for Venn)
External: minimap2, samtools
"""

import argparse, os, sys, subprocess
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import BoundaryNorm
from matplotlib import colorbar as mcolorbar

# Optional Venn
try:
    from matplotlib_venn import venn2
    HAS_VENN = True
except Exception:
    HAS_VENN = False


# ---------------- Shell helpers ----------------

def run_pipe_alignment(minimap2, samtools, ref, fastq, threads, out_bam):
    print(f"[INFO] Aligning -> {out_bam}")
    p1 = subprocess.Popen([minimap2, "-t", str(threads), "-ax", "map-ont", "--secondary=no", ref, fastq],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen([samtools, "view", "-@", str(threads), "-bS", "-"],
                          stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p3 = subprocess.Popen([samtools, "sort", "-@", str(threads), "-o", out_bam],
                          stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1.stdout.close(); p2.stdout.close()
    _, err3 = p3.communicate()
    if p1.wait()!=0 or p2.wait()!=0 or p3.wait()!=0:
        sys.stderr.write(err3.decode("utf-8"))
        raise RuntimeError("Alignment pipeline failed")

def run_cmd(cmd, desc=None):
    if desc: print(f"[INFO] {desc}: {' '.join(cmd)}")
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if res.returncode != 0:
        sys.stderr.write(res.stderr.decode("utf-8"))
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return res.stdout.decode("utf-8")

def count_fastq_reads(fq_path):
    import gzip
    if fq_path.endswith(".gz"):
        with gzip.open(fq_path, "rt") as fh:
            return sum(1 for i,_ in enumerate(fh) if i%4==0)
    else:
        with open(fq_path) as fh:
            return sum(1 for i,_ in enumerate(fh) if i%4==0)


# ---------------- QC helpers ----------------

def load_tx_len(fa_path):
    lens, tid = {}, None
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                tid = line[1:].strip().split()[0]
                lens[tid] = 0
            else:
                if tid:
                    lens[tid] += len(line.strip())
    return lens

def cigar_aligned_ref_bases(cigs):
    aligned = 0
    for op, ln in (cigs or []):
        if op in (0,7,8):  # M,=,X
            aligned += ln
    return aligned

def polyA_softclip_fraction(aln):
    """
    Inspect rightmost soft-clip as polyA (for DRS 3' relaxation).
    Returns (A_fraction, clip_len).
    """
    cigs = aln.cigartuples or []
    seq = aln.query_sequence or ""
    if not cigs or not seq: return 0.0, 0
    if cigs[-1][0] == 4:  # soft-clip at right
        ln = cigs[-1][1]
        if ln > 0:
            sc = seq[-ln:]
            a = sc.count('A') + sc.count('a')
            return a/ln, ln
    return 0.0, 0

def classify_read(aln, tx_len, tol5, tol3, mapq_min, min_cov_frac=0.30,
                  enable_polya_relax=True, polya_frac=0.7, polya_minlen=15,
                  min_tx_len=200):
    # Primary only; quality & SA filtered here
    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary: return None
    if aln.mapping_quality < mapq_min: return None
    if aln.has_tag("SA"): return None

    tid = aln.reference_name
    L = tx_len.get(tid, 0)
    if L <= 0 or L < min_tx_len: return None

    ref_start = aln.reference_start + 1
    ref_end   = aln.reference_end
    five_ok  = (ref_start <= tol5)
    three_ok = ((L - ref_end) <= tol3)

    # Optional 3' relaxation for direct RNA
    if enable_polya_relax and not three_ok:
        frac_r, ln_r = polyA_softclip_fraction(aln)
        if ln_r >= polya_minlen and frac_r >= polya_frac:
            three_ok = True

    cov_frac = cigar_aligned_ref_bases(aln.cigartuples)/L if L>0 else 0.0

    if five_ok and three_ok: cls = "full_length"
    elif five_ok and not three_ok: cls = "3prime_trunc"
    elif (not five_ok) and three_ok: cls = "5prime_trunc"
    else: cls = "internal" if cov_frac >= min_cov_frac else "low_cov"

    return {"read_id": aln.query_name, "transcript_id": tid, "tx_len": L, "class": cls}

def aligned_blocks_from_cigar(aln):
    ref_pos = aln.reference_start
    blocks, cur = [], None
    for op, ln in (aln.cigartuples or []):
        if op in (0,7,8):  # M,=,X
            if cur is None: cur = ref_pos
            ref_pos += ln
        elif op in (2,3):  # D,N -> split on reference
            if cur is not None:
                blocks.append((cur, ref_pos)); cur = None
            ref_pos += ln
        else:              # I,S,H,P -> split
            if cur is not None:
                blocks.append((cur, ref_pos)); cur = None
    if cur is not None:
        blocks.append((cur, ref_pos))
    return blocks

def bins_for_interval(L, nbins, s, e):
    if L <= 0: return ()
    s = max(0, s); e = min(L, e)
    if e <= s: return ()
    # bin index by floor((pos * nbins) / L)
    return range(int((s*nbins)//L), int(((e-1)*nbins)//L)+1)

def accumulate_to_bins(vec, L, nbins, blocks, mode="presence", min_frac=0.0):
    """
    Accumulate alignment blocks into bins.
    - mode="presence": if overlap/bin_len >= min_frac -> +1.0 (each read contributes <=1 per bin)
    - mode="depth":    add fractional coverage overlap/bin_len (0..1 per read per bin), multiple reads can sum >1
    """
    bin_len = L / nbins
    for s, e in blocks:
        s0 = max(0.0, float(s)); e0 = min(float(L), float(e))
        if e0 <= s0: continue
        for b in bins_for_interval(L, nbins, int(s0), int(e0)):
            b_start = b * bin_len
            b_end   = (b+1) * bin_len
            ov = max(0.0, min(e0, b_end) - max(s0, b_start))
            if ov <= 0: continue
            if mode == "presence":
                if (ov / bin_len) >= min_frac:
                    vec[b] += 1.0
            else:  # depth
                vec[b] += (ov / bin_len)

def presence_update(vec, L, nbins, blocks, min_frac=0.0):
    accumulate_to_bins(vec, L, nbins, blocks, mode="presence", min_frac=min_frac)

def write_gz_tsv(df, path):
    df.to_csv(path, sep="\t", index=False, compression="gzip")

def compute_forward_exclusions(forward_bam, tx_len, mapq_min, min_tx_len):
    """On forward-only BAM (primary only), count exclusions that prevent classification."""
    bam = pysam.AlignmentFile(forward_bam, "rb")
    excl_low = excl_sa = excl_ref = excl_short = 0
    total = 0
    for aln in bam.fetch(until_eof=True):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        total += 1
        tid = aln.reference_name
        L = tx_len.get(tid, 0)
        if aln.mapping_quality < mapq_min: excl_low += 1; continue
        if aln.has_tag("SA"): excl_sa += 1; continue
        if L <= 0: excl_ref += 1; continue
        if L < min_tx_len: excl_short += 1; continue
    bam.close()
    return {
        "forward_primary_for_qc": total,
        "excluded_low_mapq": excl_low,
        "excluded_split_alignment": excl_sa,
        "excluded_ref_not_in_fasta": excl_ref,
        "excluded_short_transcript": excl_short
    }


# ---------------- Visualization (presence/depth) ----------------

def qc_visualization(bam_path, ref_fa, out_prefix,
                     nbins, mapq_min, tol5, tol3, min_cov_frac,
                     min_reads, reads_cap, no_row_labels,
                     figwidth, row_height, max_figheight, dpi, clip_quantile,
                     venn_min_reads, bin_cover_min_frac, no_polya_relax, min_tx_len,
                     # metric switches
                     genebody_metric="presence", heatmap_metric="presence", read_heatmap_metric="presence",
                     depth_cap=None,
                     # read-level heatmap opts
                     heatmap_read_level=False, heatmap_max_reads=5000,
                     heatmap_read_sort="completeness", heatmap_read_subsample="stratified",
                     show_read_class_strip=False):
    tx_len = load_tx_len(ref_fa)
    bam = pysam.AlignmentFile(bam_path, "rb")

    per_read = []
    # global accumulators
    global_presence = np.zeros(nbins, dtype=np.float64)
    global_depth    = np.zeros(nbins, dtype=np.float64)
    n_reads_all = 0

    # per-transcript accumulators
    tx_presence = defaultdict(lambda: np.zeros(nbins, dtype=np.float64))
    tx_depth    = defaultdict(lambda: np.zeros(nbins, dtype=np.float64))
    tx_counts   = Counter()

    # read-level rows
    per_read_rows = []

    # --- pass & classify ---
    for aln in tqdm(bam.fetch(until_eof=True), desc="QC scan"):
        rec = classify_read(
            aln, tx_len, tol5, tol3, mapq_min, min_cov_frac,
            enable_polya_relax=(not no_polya_relax),
            min_tx_len=min_tx_len
        )
        if rec is None: continue
        per_read.append(rec)

        tid = rec["transcript_id"]; L = rec["tx_len"]
        blocks = aligned_blocks_from_cigar(aln)

        pres = np.zeros(nbins, dtype=np.float64)
        dep  = np.zeros(nbins, dtype=np.float64)

        # presence: thresholded touch
        accumulate_to_bins(pres, L, nbins, blocks, mode="presence", min_frac=bin_cover_min_frac)
        # depth: fractional base coverage
        accumulate_to_bins(dep,  L, nbins, blocks, mode="depth")

        # >>> 必修修补：单读向量截到 ≤1，防同一 read 对同一 bin 累加过量
        pres = np.minimum(pres, 1.0)
        dep  = np.minimum(dep,  1.0)

        # transcript aggregation
        tx_presence[tid] += pres
        tx_depth[tid]    += dep
        tx_counts[tid]   += 1

        # global (read-weighted)
        global_presence += pres
        global_depth    += dep
        n_reads_all += 1

        # read-level store
        per_read_rows.append({
            "read_id": rec["read_id"],
            "transcript_id": tid,
            "class": rec["class"],
            "presence": pres,
            "depth": dep,
            "comp_presence": float(pres.mean()),
            "comp_depth": float(dep.mean())
        })

    bam.close()
    classified_reads = len(per_read)
    if classified_reads == 0:
        print("[WARN] No passing reads for QC; skip plots.")
        return 0

    # --- per-read table ---
    df = pd.DataFrame(per_read)
    write_gz_tsv(df, out_prefix + ".read_class.tsv.gz")

    # --- per-transcript stats ---
    agg = df.groupby("transcript_id").agg(
        tx_len=('tx_len','first'),
        n_reads=('read_id','count'),
        n_FL=('class', lambda s: (s=='full_length').sum()),
        n_5tr=('class', lambda s: (s=='5prime_trunc').sum()),
        n_3tr=('class', lambda s: (s=='3prime_trunc').sum()),
        n_internal=('class', lambda s: (s=='internal').sum()),
        n_lowcov=('class', lambda s: (s=='low_cov').sum())
    ).reset_index()
    agg['FL_rate'] = (agg['n_FL'] / agg['n_reads']).fillna(0).round(6)
    write_gz_tsv(agg, out_prefix + ".tx_stats.tsv.gz")

    # ---------- 1) Classification bar ----------
    cls = df["class"].value_counts().reindex(
        ["full_length","5prime_trunc","3prime_trunc","internal","low_cov"], fill_value=0)
    plt.figure(figsize=(6,4))
    plt.bar(cls.index, cls.values)
    plt.ylabel("Reads"); plt.title("Read classification")
    plt.xticks(rotation=20); plt.tight_layout()
    plt.savefig(out_prefix + ".class_bar.png", dpi=dpi); plt.close()

    # ---------- 2) Gene body curves (presence/depth) ----------
    x = (np.arange(nbins)+0.5)/nbins * 100.0
    metric = genebody_metric  # "presence" or "depth"

    if metric == "presence":
        y_read = global_presence / max(1, n_reads_all)
        frac_per_tx = []
        for tid in tx_presence:
            r = max(1, tx_counts[tid])
            frac_per_tx.append(tx_presence[tid] / r)
        y_tx = np.mean(np.vstack(frac_per_tx), axis=0) if frac_per_tx else np.zeros_like(y_read)
        ylabel = "Presence fraction (0..1)"
    else:  # depth
        y_read = global_depth / max(1, n_reads_all)
        frac_per_tx = []
        for tid in tx_depth:
            r = max(1, tx_counts[tid])
            frac_per_tx.append(tx_depth[tid] / r)
        y_tx = np.mean(np.vstack(frac_per_tx), axis=0) if frac_per_tx else np.zeros_like(y_read)
        if depth_cap and depth_cap > 0:
            y_read = np.clip(y_read, None, depth_cap)
            y_tx   = np.clip(y_tx,   None, depth_cap)
        ylabel = "Mean depth per read"

    plt.figure(figsize=(8,4))
    plt.plot(x, y_read, label="Read-weighted")
    plt.plot(x, y_tx,   label="Transcript-weighted")
    plt.xlabel("Transcript body (5′ → 3′, %)")
    plt.ylabel(ylabel)
    plt.ylim(bottom=0)
    plt.title(f"Gene body ({metric}), NBINS={nbins}")
    plt.legend(); plt.tight_layout()
    plt.savefig(out_prefix + f".genebody_{metric}.png", dpi=dpi); plt.close()

    # ---------- 3) Transcript-level heatmap (presence/depth) ----------
    tids_all, comp_all, mat_rows, reads_list = [], [], [], []
    for tid in tx_counts:
        r = max(1, tx_counts[tid])
        if heatmap_metric == "presence":
            vec = tx_presence[tid] / r
            comp = float(vec.mean())
        else:
            vec = tx_depth[tid] / r
            comp = float(vec.mean())
            if depth_cap and depth_cap > 0:
                vec = np.clip(vec, None, depth_cap)
        tids_all.append(tid); comp_all.append(comp); mat_rows.append(vec); reads_list.append(r)

    if mat_rows:
        mat = np.vstack(mat_rows)
        df_heat = pd.DataFrame({
            "transcript_id": tids_all,
            "completeness": comp_all,
            "n_reads": reads_list
        })
        df_heat = df_heat.merge(agg[["transcript_id","tx_len","FL_rate"]], on="transcript_id", how="left")
        df_heat = df_heat[df_heat["n_reads"] >= min_reads].copy()

        if not df_heat.empty:
            # reorder and sort
            order_idx = [tids_all.index(t) for t in df_heat["transcript_id"]]
            mat = mat[order_idx, :]
            sidx = np.argsort(-df_heat["completeness"].to_numpy())
            df_heat = df_heat.iloc[sidx].reset_index(drop=True)
            mat = mat[sidx, :]

            # save matrix table
            mdf = pd.DataFrame(mat, columns=[f"bin_{i}" for i in range(nbins)])
            mdf.insert(0, "transcript_id", df_heat["transcript_id"])
            mdf.insert(1, "completeness", df_heat["completeness"])
            mdf.insert(2, "n_reads", df_heat["n_reads"])
            mdf.to_csv(out_prefix + f".tx_{heatmap_metric}_matrix.tsv.gz", sep="\t", index=False, compression="gzip")

            # robust optional clipping (secondary)
            if 0.0 < clip_quantile < 1.0:
                hi = np.quantile(mat, clip_quantile)
                mat = np.clip(mat, a_min=None, a_max=hi)

            # right-side reads strip (capped)
            cap = max(1, int(reads_cap))
            nr_disp = np.minimum(df_heat["n_reads"].to_numpy(), cap) / cap
            nr_strip = nr_disp.reshape(-1,1)

            # layout: heatmap | strip | colorbars
            fig_h = min(max_figheight, max(3.0, row_height * mat.shape[0]))
            fig_w = figwidth
            fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
            gs = GridSpec(1, 3, width_ratios=[20, 1, 0.9], wspace=0.05)

            # heatmap
            ax0 = fig.add_subplot(gs[0,0])
            im0 = ax0.imshow(mat, aspect='auto', interpolation='nearest')
            ax0.set_xticks([])
            if no_row_labels:
                ax0.set_yticks([]); ax0.set_ylabel("")
                ax0.tick_params(left=False); ax0.spines['left'].set_visible(False)
            else:
                ax0.set_yticks(np.arange(len(df_heat)))
                ax0.set_yticklabels(df_heat["transcript_id"], fontsize=6)
            title_metric = "Presence frac" if heatmap_metric=="presence" else "Depth/Read"
            ax0.set_title(f"Genebody heatmap ({title_metric})\nSorted by completeness", fontsize=10)

            # reads strip
            ax1 = fig.add_subplot(gs[0,1])
            im1 = ax1.imshow(nr_strip, aspect='auto', interpolation='nearest', cmap='Greys')
            ax1.set_xticks([]); ax1.set_yticks([])
            ax1.set_title(f"n_reads (cap={cap})", fontsize=9, pad=6)

            # far-right colorbars
            axc = fig.add_subplot(gs[0,2]); axc.axis('off')
            bbox = axc.get_position()
            # main colorbar
            cax0 = fig.add_axes([bbox.x0 + 0.10*bbox.width, bbox.y0 + 0.55*bbox.height,
                                 0.30*bbox.width, 0.35*bbox.height])
            cb0 = fig.colorbar(im0, cax=cax0, orientation='vertical')
            cb0.ax.tick_params(labelsize=8)
            cb0.set_label('Presence frac' if heatmap_metric=="presence" else 'Depth/Read', fontsize=8)
            # reads colorbar
            cax1 = fig.add_axes([bbox.x0 + 0.10*bbox.width, bbox.y0 + 0.08*bbox.height,
                                 0.30*bbox.width, 0.35*bbox.height])
            cb1 = fig.colorbar(im1, cax=cax1, orientation='vertical')
            cb1.ax.tick_params(labelsize=8); cb1.set_ticks([0.0, 1.0]); cb1.set_ticklabels(['0', f'{cap}'])

            plt.tight_layout()
            fig.savefig(out_prefix + f".heatmap_{heatmap_metric}.png", dpi=dpi)
            plt.close(fig)
        else:
            print(f"[WARN] No transcripts with n_reads >= {min_reads} to draw heatmap.")
    else:
        print("[WARN] No transcripts aggregated for heatmap.")

    # ---------- 3b) Read-level heatmap (presence/depth, optional) ----------
    if heatmap_read_level and per_read_rows:
        metric_r = read_heatmap_metric  # "presence" or "depth"
        df_r = pd.DataFrame([
            {"read_id":r["read_id"], "transcript_id":r["transcript_id"],
             "class":r["class"],
             "comp": r["comp_presence"] if metric_r=="presence" else r["comp_depth"]}
            for r in per_read_rows
        ])
        mat_r = np.vstack([
            (r["presence"] if metric_r=="presence" else r["depth"])
            for r in per_read_rows
        ])

        if metric_r=="depth" and depth_cap and depth_cap>0:
            mat_r = np.clip(mat_r, None, depth_cap)

        # subsample
        max_rows = max(1, int(heatmap_max_reads))
        if mat_r.shape[0] > max_rows:
            if heatmap_read_subsample == "head":
                sel_idx = np.arange(max_rows)
            elif heatmap_read_subsample == "random":
                sel_idx = np.random.choice(mat_r.shape[0], size=max_rows, replace=False)
            else:  # stratified by class
                sel_idx_list = []
                for _, sub in df_r.groupby("class"):
                    k = max(1, int(max_rows * len(sub) / len(df_r)))
                    sel_idx_list.append(sub.sample(n=min(k, len(sub)), random_state=1).index.to_numpy())
                sel_idx = np.concatenate(sel_idx_list) if sel_idx_list else np.arange(max_rows)
                if len(sel_idx) > max_rows:
                    sel_idx = np.random.choice(sel_idx, size=max_rows, replace=False)
            df_r = df_r.iloc[sel_idx].reset_index(drop=True)
            mat_r = mat_r[sel_idx, :]

        # sort
        if heatmap_read_sort == "class":
            order = np.argsort(df_r["class"].to_numpy())
        else:
            order = np.argsort(-df_r["comp"].to_numpy())
        df_r = df_r.iloc[order].reset_index(drop=True)
        mat_r = mat_r[order, :]

        # draw
        fig_h = min(max_figheight, max(3.0, row_height * mat_r.shape[0]))
        fig_w = figwidth
        if show_read_class_strip:
            gs = GridSpec(1, 3, width_ratios=[20, 1, 0.9], wspace=0.05)
        else:
            gs = GridSpec(1, 2, width_ratios=[20, 0.9], wspace=0.05)

        fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
        ax0 = fig.add_subplot(gs[0,0])
        im0 = ax0.imshow(mat_r, aspect='auto', interpolation='nearest')
        ax0.set_xticks([])
        ax0.set_yticks([]); ax0.set_ylabel("")
        ax0.tick_params(left=False); ax0.spines['left'].set_visible(False)
        ax0.set_title(f"Read-level genebody heatmap ({metric_r})", fontsize=10)

        if show_read_class_strip:
            class2idx = {"full_length":0,"5prime_trunc":1,"3prime_trunc":2,"internal":3,"low_cov":4}
            cmap_classes = plt.get_cmap("tab10")
            class_strip_idx = np.array([class2idx.get(c,4) for c in df_r["class"]], dtype=int).reshape(-1,1)
            ax1 = fig.add_subplot(gs[0,1])
            im1 = ax1.imshow(class_strip_idx, aspect='auto', interpolation='nearest',
                             cmap=cmap_classes, vmin=0, vmax=4)
            ax1.set_xticks([]); ax1.set_yticks([]); ax1.set_title("class", fontsize=9, pad=6)

            axc = fig.add_subplot(gs[0,2]); axc.axis('off')
            bbox = axc.get_position()
            cax0 = fig.add_axes([bbox.x0 + 0.10*bbox.width, bbox.y0 + 0.55*bbox.height,
                                 0.30*bbox.width, 0.35*bbox.height])
            cb0 = fig.colorbar(im0, cax=cax0, orientation='vertical')
            cb0.ax.tick_params(labelsize=8)
            cb0.set_label('Presence' if metric_r=="presence" else 'Depth/Read', fontsize=8)
            cax1 = fig.add_axes([bbox.x0 + 0.10*bbox.width, bbox.y0 + 0.08*bbox.height,
                                 0.30*bbox.width, 0.35*bbox.height])
            norm = BoundaryNorm(np.arange(-0.5,5.5,1), cmap_classes.N)
            cb1 = mcolorbar.ColorbarBase(cax1, cmap=cmap_classes, norm=norm, orientation='vertical')
            cb1.set_ticks([0,1,2,3,4]); cb1.set_ticklabels(["FL","5T","3T","INT","LOW"])
        else:
            axc = fig.add_subplot(gs[0,1]); axc.axis('off')
            bbox = axc.get_position()
            cax0 = fig.add_axes([bbox.x0 + 0.15*bbox.width, bbox.y0 + 0.30*bbox.height,
                                 0.35*bbox.width, 0.40*bbox.height])
            cb0 = fig.colorbar(im0, cax=cax0, orientation='vertical')
            cb0.ax.tick_params(labelsize=8)
            cb0.set_label('Presence' if metric_r=="presence" else 'Depth/Read', fontsize=8)

        plt.tight_layout()
        fig.savefig(out_prefix + f".heatmap_{metric_r}.readlevel.png", dpi=dpi)
        plt.close(fig)

    # ---------- 4) Venn ----------
    if HAS_VENN:
        detected = set(agg.loc[agg["n_reads"] >= venn_min_reads, "transcript_id"].tolist())
        db_set   = set(load_tx_len(ref_fa).keys())
        plt.figure(figsize=(4.5,4.0), dpi=dpi)
        venn2(subsets=(len(detected-db_set), len(db_set-detected), len(detected & db_set)),
              set_labels=('Detected', 'In DB'))
        plt.title(f"Detected vs DB (min_reads={venn_min_reads})")
        plt.tight_layout()
        plt.savefig(out_prefix + ".venn_detected_vs_db.png", dpi=dpi)
        plt.close()
    else:
        print("[INFO] matplotlib-venn not installed; skip Venn.")

    return classified_reads


# ---------------- Main ----------------

def main():
    ap = argparse.ArgumentParser(description="Alignment + counts + transcript/read-level QC with presence/depth metrics")
    ap.add_argument("--ref", required=True, help="Transcript FASTA")
    ap.add_argument("--fastq", required=True, help="FASTQ/FASTQ.GZ")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--use-forward-bam", action="store_true", help="Use forward-only BAM for QC")

    # QC params
    ap.add_argument("--nbins", type=int, default=100)
    ap.add_argument("--mapq-min", type=int, default=20)
    ap.add_argument("--tol5", type=int, default=20)
    ap.add_argument("--tol3", type=int, default=20)
    ap.add_argument("--min-cov-frac", type=float, default=0.30)

    # Robustness knobs
    ap.add_argument("--min-tx-len", type=int, default=200, help="Filter transcripts shorter than this length")
    ap.add_argument("--bin-cover-min-frac", type=float, default=0.0,
                    help="Minimal fraction overlap required to mark a bin as present (0..1)")
    ap.add_argument("--no-polya-relax", action="store_true",
                    help="Disable 3' polyA soft-clip relaxation")

    # Visualization sizing
    ap.add_argument("--min-reads", type=int, default=10)
    ap.add_argument("--reads-cap", type=int, default=100)
    ap.add_argument("--no-row-labels", action="store_true")
    ap.add_argument("--figwidth", type=float, default=24.0)
    ap.add_argument("--row-height", type=float, default=0.10)
    ap.add_argument("--max-figheight", type=float, default=48.0)
    ap.add_argument("--dpi", type=int, default=180)
    ap.add_argument("--clip-quantile", type=float, default=0.995)
    ap.add_argument("--venn-min-reads", type=int, default=1)

    # Metric switches
    ap.add_argument("--genebody-metric", choices=["presence","depth"], default="presence",
                    help="Metric for gene body curves")
    ap.add_argument("--heatmap-metric", choices=["presence","depth"], default="presence",
                    help="Metric for transcript-level heatmap")
    ap.add_argument("--read-heatmap-metric", choices=["presence","depth"], default="presence",
                    help="Metric for read-level heatmap")
    ap.add_argument("--depth-cap", type=float, default=None,
                    help="Cap for depth-based plotting (e.g., 5.0)")

    # Read-level heatmap options
    ap.add_argument("--heatmap-read-level", action="store_true",
                    help="Draw read-level genebody heatmap (each row = one read)")
    ap.add_argument("--heatmap-max-reads", type=int, default=5000,
                    help="Maximum number of reads to display in read-level heatmap")
    ap.add_argument("--heatmap-read-sort", choices=["completeness","class"], default="completeness",
                    help="Sort reads by completeness (default) or by class")
    ap.add_argument("--heatmap-read-subsample", choices=["head","random","stratified"], default="stratified",
                    help="Subsample strategy when reads > max: head/random/stratified-by-class")
    ap.add_argument("--show-read-class-strip", action="store_true",
                    help="Add a right-side class color strip for read-level heatmap")

    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ALL = os.path.join(args.outdir, "aln.tx.bam")
    FWD = os.path.join(args.outdir, "aln.tx.forward.bam")

    # 1) Alignment (ALL)
    run_pipe_alignment(args.minimap2, args.samtools, args.ref, args.fastq, args.threads, ALL)
    run_cmd([args.samtools, "index", ALL], "Index ALL")

    # 2) Forward-only BAM
    run_cmd([args.samtools, "view", "-@", str(args.threads), "-b", "-F", "16", ALL, "-o", FWD], "Filter forward-only")
    run_cmd([args.samtools, "index", FWD], "Index FWD")

    # 3) Counts & exclusions
    input_reads = count_fastq_reads(args.fastq)
    PRIMARY_MASK = 2308  # exclude: unmapped(4)+secondary(256)+supplementary(2048)
    mapped_all  = int(run_cmd([args.samtools,"view","-@",str(args.threads),"-c","-F",str(PRIMARY_MASK),ALL]).strip())
    mapped_fwd  = int(run_cmd([args.samtools,"view","-@",str(args.threads),"-c","-F",str(PRIMARY_MASK+16),ALL]).strip())
    mapped_rev  = int(run_cmd([args.samtools,"view","-@",str(args.threads),"-c","-f","16","-F",str(PRIMARY_MASK),ALL]).strip())
    mapped_fwd_bam = int(run_cmd([args.samtools,"view","-@",str(args.threads),"-c","-F",str(PRIMARY_MASK),FWD]).strip())

    tx_len_dict = load_tx_len(args.ref)
    excl = compute_forward_exclusions(FWD, tx_len_dict, args.mapq_min, args.min_tx_len)

    # 4) QC on chosen BAM
    qc_bam = FWD if args.use_forward_bam else ALL
    out_prefix = os.path.join(args.outdir, "sample_tx_qc")
    classified = qc_visualization(
        bam_path=qc_bam, ref_fa=args.ref, out_prefix=out_prefix,
        nbins=args.nbins, mapq_min=args.mapq_min, tol5=args.tol5, tol3=args.tol3,
        min_cov_frac=args.min_cov_frac,
        min_reads=args.min_reads, reads_cap=args.reads_cap,
        no_row_labels=args.no_row_labels,
        figwidth=args.figwidth, row_height=args.row_height, max_figheight=args.max_figheight,
        dpi=args.dpi, clip_quantile=args.clip_quantile,
        venn_min_reads=args.venn_min_reads,
        bin_cover_min_frac=args.bin_cover_min_frac,
        no_polya_relax=args.no_polya_relax,
        min_tx_len=args.min_tx_len,
        genebody_metric=args.genebody_metric,
        heatmap_metric=args.heatmap_metric,
        read_heatmap_metric=args.read_heatmap_metric,
        depth_cap=args.depth_cap,
        heatmap_read_level=args.heatmap_read_level,
        heatmap_max_reads=args.heatmap_max_reads,
        heatmap_read_sort=args.heatmap_read_sort,
        heatmap_read_subsample=args.heatmap_read_subsample,
        show_read_class_strip=args.show_read_class_strip
    )

    # 5) Summary TSV
    with open(os.path.join(args.outdir, "alignment_counts.tsv"), "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"input_fastq_reads\t{input_reads}\n")
        fh.write(f"mapped_primary_all\t{mapped_all}\n")
        fh.write(f"mapped_primary_forward\t{mapped_fwd}\n")
        fh.write(f"mapped_primary_reverse\t{mapped_rev}\n")
        fh.write(f"mapped_primary_forward_bam\t{mapped_fwd_bam}\n")
        for k,v in excl.items(): fh.write(f"{k}\t{v}\n")
        fh.write(f"classified_reads\t{classified}\n")

    print("[DONE] Outputs in", args.outdir)
    print("  - aln.tx.bam(.bai), aln.tx.forward.bam(.bai)")
    print("  - alignment_counts.tsv")
    print("  - sample_tx_qc.class_bar.png")
    print("  - sample_tx_qc.genebody_presence.png OR sample_tx_qc.genebody_depth.png")
    print("  - sample_tx_qc.heatmap_presence.png OR sample_tx_qc.heatmap_depth.png")
    print("  - sample_tx_qc.heatmap_presence.readlevel.png / heatmap_depth.readlevel.png (if --heatmap-read-level)")
    print("  - sample_tx_qc.venn_detected_vs_db.png (if matplotlib-venn installed)")
    print("  - sample_tx_qc.read_class.tsv.gz")
    print("  - sample_tx_qc.tx_stats.tsv.gz")
    print("  - sample_tx_qc.tx_presence_matrix.tsv.gz OR tx_depth_matrix.tsv.gz")

if __name__ == "__main__":
    main()
