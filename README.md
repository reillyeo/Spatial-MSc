# Spatial-MSc

This workflow is designed for primary analysis of imaging-based spatial transcriptomic data.
The input data in the src directory includes DAPI images of each tissue section in .tif format and coordinates of transcript reads in 3-column .csv files ("gene", "x", "y").

Images and reads originating from the same tissue section share the same file prefix, referenced in the Snakefile as 'SECTION_PREFIXES'.

The src directory here also contains the file "scRNA_refdata". This is a seurat object containing single cell RNA-seq data that is used as a reference when assigning cell-type labels to spatial data (loaded in R using qread).

The pipeline first segments the microscopy images using Cellpose, followed by Baysor, and assigns transcript reads to cells based on the generated segmentation. Baysor configuration parameters can be found in src/baysor_configs.toml

Following the assignment of transcripts to cells, typical single cell analysis steps can be performed (dimension reduction, cell type assignment)
