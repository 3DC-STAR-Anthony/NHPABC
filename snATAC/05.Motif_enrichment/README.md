Motif Enrichment Analysis of NHPABC

- Motif Database: JASPAR2022 (JB14.1).
- TF Filtering: Excluded TFs expressed in <5% of cells (based on matched snRNA-seq clusters).
- Peak Annotation: Annotated peaks with known motifs using ArchR's addMotifAnnotations(default settings).
- Background Selection: Generated GC-matched background peaks for each peak using ArchR's getBgdPeak(default settings).
- Enrichment Calculation: Computed motif enrichment in target peak sets versus background using ArchR's computeEnrichment.

# Step 1: Preparations of Input File
