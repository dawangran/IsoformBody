# IsoformBody

**IsoformBody** 是一个面向 RNA 直测 (direct RNA sequencing) 数据的 **转录本覆盖度与完整性评估工具**。  
它以 isoform (转录本) 为分析单元，绘制 **genebody 覆盖曲线、热图、分类统计**，帮助快速发现 5′/3′ 截短偏倚与覆盖不均。

---

## ✨ 功能
- **自动对齐**：调用 minimap2 + samtools，比对 reads → transcript FASTA
- **分类**：基于起止位置和覆盖度，将 reads 分为 full-length / 5′-trunc / 3′-trunc / internal / low_cov
- **可视化**：
  - Reads 分类柱状图
  - Gene body 曲线（presence / depth 双口径）
  - Transcript-level 热图（可选 reads strip）
  - Read-level 热图（每条 read 可视化）
  - Venn 图（检测到 vs 数据库转录本）
- **统计输出**：
  - 每条 read 分类表
  - 每个转录本的 reads 数和 FL 比例
  - 热图矩阵 (presence / depth)
  - alignment counts 汇总

---

## 📦 安装
```bash
conda create -n isoformbody python=3.10
conda activate isoformbody
pip install -r requirements.txt
```

外部依赖需提前安装：
- minimap2
- samtools

---

## 🚀 使用示例

### Presence 模式（默认，推荐先跑）
```bash
python pipeline_tx_qc.py --ref transcripts.fa --fastq sample.fastq.gz --outdir results_presence --threads 16 --use-forward-bam --nbins 100 --mapq-min 20 --tol5 20 --tol3 20 --min-cov-frac 0.30 --min-tx-len 300 --bin-cover-min-frac 0.10 --min-reads 10 --reads-cap 100 --no-row-labels --figwidth 24 --row-height 0.10 --max-figheight 48 --dpi 180 --clip-quantile 0.995 --venn-min-reads 1 --genebody-metric presence --heatmap-metric presence --read-heatmap-metric presence --heatmap-read-level --show-read-class-strip
```

### Depth 模式
```bash
python pipeline_tx_qc.py --ref transcripts.fa --fastq sample.fastq.gz --outdir results_depth --threads 16 --use-forward-bam --nbins 100 --mapq-min 20 --tol5 20 --tol3 20 --min-cov-frac 0.30 --min-tx-len 300 --bin-cover-min-frac 0.10 --min-reads 10 --reads-cap 100 --no-row-labels --figwidth 24 --row-height 0.10 --max-figheight 48 --dpi 180 --clip-quantile 0.995 --venn-min-reads 1 --genebody-metric depth --heatmap-metric depth --read-heatmap-metric depth --depth-cap 5 --heatmap-read-level --show-read-class-strip
```

---

## 📊 输出文件
- 比对文件：aln.tx.bam(.bai), aln.tx.forward.bam(.bai)
- 表格：alignment_counts.tsv, read_class.tsv.gz, tx_stats.tsv.gz, tx_matrix.tsv.gz
- 图形：class_bar.png, genebody_presence.png, heatmap_presence.png, venn_detected_vs_db.png

---

## 📄 License
MIT License
