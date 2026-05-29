# Figure Generation Code for NHPABC

Code for generating figures in **"Multimodal brain cell atlas across the adult macaque lifespan"**.

## About

This repository contains R scripts used to produce all publication figures for the NHPABC (Non-Human Primate Aging Brain Cell atlas) project — a multimodal single-cell atlas of 2,955,873 nuclei from eight brain regions of 23 female cynomolgus macaques across the adult lifespan.

## Repository Structure

```
figures/              → Main figure code (Figure 1–7)
supplemental_figures/ → Supplemental figure code (Figure S1–S7)
```

## Code Organization

Scripts are organized by **visualization type**. Each script generates one type of plot and may cover multiple figure panels across different figures.

### Main Figures

| # | Script | Type | Covers |
|---|--------|------|--------|
| 1 | `plot_schematic.R` | schematic | Fig1A |
| 2 | `plot_umap.R` | UMAP embedding | Fig1B-C, Fig2A, Fig7A |
| 3 | `plot_scatter.R` | scatter plot | Fig1D |
| 4 | `plot_barplot.R` | bar plot | Fig1E-F, Fig2G-H, Fig3A,C, Fig5B,H, Fig6A, Fig7E |
| 5 | `plot_heatmap.R` | heatmap | Fig1G-H, Fig2M, Fig3F, Fig4A,C,E, Fig5C,G, Fig6D, Fig7B,D |
| 6 | `plot_dotplot.R` | dot plot | Fig2B, Fig3D-E,G-H, Fig4B,D,F, Fig5D,I, Fig6B, Fig7C |
| 7 | `plot_composite.R` | multi-panel | Fig2C-F, Fig4G,L-M, Fig6F |
| 8 | `plot_lineplot.R` | line/trend | Fig2I-L, Fig3I-J, Fig4H-K, Fig5E,J, Fig6E |
| 9 | `plot_boxplot.R` | box plot | Fig2N |
| 10 | `plot_histogram.R` | histogram | Fig3B |
| 11 | `plot_network.R` | network | Fig5A,F |
| 12 | `plot_genome_track.R` | genome track | Fig6C |

### Supplemental Figures

| # | Script | Type | Covers |
|---|--------|------|--------|
| 1 | `plot_violin.R` | violin plot | FigS1A, FigS6L |
| 2 | `plot_barplot.R` | bar plot | FigS1B, FigS2C, FigS3A,C,D,G,K, FigS5C, FigS6B,K, FigS7E,G |
| 3 | `plot_umap.R` | UMAP embedding | FigS1C-D, FigS2A-B, FigS7A |
| 4 | `plot_boxplot.R` | box plot | FigS2D-E |
| 5 | `plot_histogram.R` | histogram | FigS3B |
| 6 | `plot_dotplot.R` | dot plot | FigS3E-F,L-N, FigS4B,D,F,H, FigS5D, FigS6I, FigS7F |
| 7 | `plot_heatmap.R` | heatmap | FigS3H-J, FigS4A,C,E,G,K-L, FigS5B, FigS6E,G, FigS7B |
| 8 | `plot_lineplot.R` | line/trend | FigS3I,M,O-P, FigS5E, FigS6D,F, FigS7D |
| 9 | `plot_composite.R` | multi-panel | FigS4I-J, FigS6J |
| 10 | `plot_network.R` | network | FigS5A,F |
| 11 | `plot_genome_track.R` | genome track | FigS6A,H |
| 12 | `plot_scatter.R` | scatter plot | FigS6C, FigS7C |

## Dependencies

- R >= 4.0
- Key packages: `Seurat`, `Signac`, `ggplot2`, `ComplexHeatmap`, `dplyr`, `GenomicRanges`

## Data

Processed data are available through the [NHPABC portal](https://db.cngb.org/stomics/nhpabc/). See the paper's **Data Availability** section for details.

## Citation

> Xiao Zhang, Guangyao Lai, Xiangyu Guo, et al. *Multimodal brain cell atlas across the adult macaque lifespan*. **Cell**, 2025.
