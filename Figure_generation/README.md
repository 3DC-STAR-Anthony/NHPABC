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
| 1 | `plot_umap.R` | UMAP embedding | Fig1B, Fig2A left, Fig2C left, Fig2E left, Fig7C |
| 2 | `plot_scatter.R` | scatter plot | Fig1C, Fig3C, Fig4A, Fig5A, Fig6A |
| 3 | `plot_barplot.R` | bar plot | Fig2A, Fig2C, Fig2E, Fig3C-D, Fig3G-H, Fig4A, Fig5B, Fig6A, Fig7B |
| 4 | `plot_heatmap.R` | heatmap | Fig1D, Fig3F, Fig4B-C, Fig4E-F, Fig4G, Fig4I-J, Fig4L-M, Fig6D, Fig7A, Fig7D-E |
| 5 | `plot_dotplot.R` | dot plot | Fig2B, Fig2D, Fig2F-G, Fig5D, Fig5F, Fig5H, Fig5J |
| 6 | `plot_line&dotplot.R` | line/trend | Fig3E, Fig3I, Fig4D, Fig4H, Fig4K, Fig4N, Fig5E, Fig5G, Fig5I, Fig5K |
| 7 | `plot_vlnplot.R` | violinplot | Fig3J, Fig4D, Fig4H, Fig4K, Fig4H, Fig5E, Fig5G, Fig5I, Fig5K |
| 8 | `plot_boxplot.R` | box plot | Fig3E, Fig5C |
| 9 | `plot_histogram.R` | histogram | Fig3B |
| 10 | `plot_volcanoplot.R` | volcano | Fig2I, Fig2K |
| 11 | `plot_footprinting.R` | footprinting | Fig6E |
| 12 | `plot_genome_track.R` | genome track | Fig6B, Fig6C |


### Supplemental Figures

| # | Script | Type | Covers |
|---|--------|------|--------|
| 1 | `plot_umap.R` | UMAP embedding | FigS2A-C |
| 2 | `plot_scatter.R` | scatter plot | FigS7A |
| 3 | `plot_barplot.R` | bar plot | FigS2B, FigS2C, FigS3F, FigS3L, FigS3N, FigS5B, FigS6D, FigS6K, FigS6M, FigS7C, FigS7E |
| 4 | `plot_heatmap.R` | heatmap | FigS3I, FigS4F, FigS4I-N, FigS5C, FigS6C, FigS6G, FigS6O |
| 5 | `plot_dotplot.R` | dot plot | FigS1C |
| 6 | `plot_line&dotplot.R` | line/trend | FigS3M, FigS3O-P, FigS4G-H, FigS6B |
| 7 | `plot_vlnplot.R` | violinplot | FigS3M, FigS3O-P, FigS4G-H |
| 8 | `plot_boxplot.R` | box plot | FigS1A, FigS1B, FigS2D, FigS2E, FigS3A-D, FigS3J, FigS4A-D, FigS5A, FigS6I-J, FigS6L, FigS7D |
| 9 | `plot_histogram.R` | histogram | FigS3G, FigS4E |
| 10 | `plot_volcanoplot.R` | Volcano | FigS2E, FigS6A |
| 11 | `plot_vennplot.R` | venn | FigS3K, FigS3H, FigS6E |
| 12 | `plot_composite.R` | multi-panel | FigS6F |
| 13 | `plot_pieplot.R` | pie plot | FigS6H |
| 14 | `plot_genome_track.R` | genome track | FigS6N, FigS7B |

## Dependencies

- R >= 4.0
- Key packages: `Seurat`, `Signac`, `ggplot2`, `ComplexHeatmap`, `dplyr`, `GenomicRanges`

## Data

Processed data are available through the [NHPABC portal](https://db.cngb.org/stomics/nhpabc/). See the paper's **Data Availability** section for details.

