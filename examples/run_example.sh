#!/bin/bash
# Example run for IsoformBody

python ../pipeline_tx_qc.py   --ref transcripts.fa   --fastq sample.fastq.gz   --outdir results_example   --threads 8   --use-forward-bam   --nbins 100   --mapq-min 20 --tol5 20 --tol3 20   --min-cov-frac 0.30 --min-tx-len 300 --bin-cover-min-frac 0.10   --min-reads 10 --reads-cap 100   --no-row-labels   --figwidth 24 --row-height 0.10 --max-figheight 48 --dpi 180   --clip-quantile 0.995   --venn-min-reads 1   --genebody-metric presence   --heatmap-metric presence   --read-heatmap-metric presence   --heatmap-read-level   --show-read-class-strip
