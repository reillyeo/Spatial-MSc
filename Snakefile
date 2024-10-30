SECTION_PREFIXES = ["B8_6B_1A", "B8_6B_1B", "B8_6B_1C", "B1_11B_1A", "B3_9B_1B", "B3_9B_1D", "B5_5C_1B", "B5_5C_1A", "B4_6B_1C", "B6_6B_1C", "B7_9B_1B", "B4_8C_1A", "B4_8C_1C", "A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2"]



rule all:
    input:
        expand("out/cellpose_output/{section}_cp.csv", section = SECTION_PREFIXES),
        expand("out/baysor_results/{section}.csv", section = SECTION_PREFIXES),
        expand("src/data_post_seg/{section}.csv, section = SECTION_PREFIXES"),
        expand("out/seurat_objects/{section}_seurat", section = SECTION_PREFIXES),
        expand("out/labelled_spatial/{section}_labelled_corral", section = SECTION_PREFIXES),
        expand("out/labelled_spatial/{section}_labelled_sct", section = SECTION_PREFIXES),
        expand("out/labelled_fine/{section}_labelled_sct", section = SECTION_PREFIXES),
        expand("out/labelled_fine/{section}_labelled_corral", section = SECTION_PREFIXES),
        expand("out/sctype_labelled/{section}_sctyped", section = SECTION_PREFIXES)


# run cellpose segmentation on all dapi images and add segmentation labels to reads dataframes
rule run_cellpose:
    input:
        dapi = "src/raw_data/dapi_images/{section}.tif",
        reads = "src/raw_data/transcript_coords/{section}_reads.csv"
    output:
        "out/cellpose_output/{section}_cp.csv"
    shell:
        "python src/scripts/cellpose_seg.py {input.dapi} {input.reads} {output}"


# run baysor segmentation on all sections using prior cellpose segmentation
rule run_baysor:
    input:
        df = "out/cellpose_output/{section}_cp.csv",
        config = "src/baysor_configs.toml"
    output:
        "out/baysor_results/{section}"
    shell:
        "./baysor/bin/baysor run -c {input.config} {input.df} :cp_label -o {output}"


# read baysor results and convert to single-cell format csv that can be read into seurat
rule run_data_gen:
    input:
        "out/baysor_results/{section}"
    output:
        "src/data_post_seg/{section}.csv"
    shell:
        "python src/scripts/data_gen.py {input} {output}"


# run various dimension reductions on seurat objects
rule run_dimreduc:
    input:
        df = "src/data_post_seg/{section}.csv"
    output:
        "out/seurat_objects/{section}_seurat"
    shell:
        "Rscript src/scripts/dimreducs.R {input} {output}"


# following 5 rules all perform cell type labelling in different manners or based on different dimreducs

rule run_labeltransfer_corral:
    input:
        st = "out/seurat_objects/{section}_seurat",
        sc = "src/data/scRNA_refdata"
    output:
        "out/labelled_spatial/{section}_labelled_corral"
    shell:
        "Rscript src/scripts/label_transfer_corral.R {input.st} {input.sc} {output}"



rule run_labeltransfer_sct:
    input:
        st = "out/seurat_objects/{section}_seurat",
        sc = "src/data/scRNA_refdata"
    output:
        "out/labelled_spatial/{section}_labelled_sct"
    shell:
        "Rscript src/scripts/label_transfer_sct.R {input.st} {input.sc} {output}"



rule run_labeltransfer_corral_fine:
    input:
        st = "out/seurat_objects/{section}_seurat",
        sc = "src/data/sc_final_labels"
    output:
        "out/labelled_fine/{section}_labelled_corral"
    shell:
        "Rscript src/scripts/label_transfer_corral_fine.R {input.st} {input.sc} {output}"



rule run_labeltransfer_sct_fine:
    input:
        st = "out/seurat_objects/{section}_seurat",
        sc = "src/data/sc_final_labels"
    output:
        "out/labelled_fine/{section}_labelled_sct"
    shell:
        "Rscript src/scripts/label_transfer_sct_fine.R {input.st} {input.sc} {output}"



rule run_sctype:
    input:
        st = "out/seurat_objects/{section}_seurat",
        sc = "src/data/sc_final_labels"
    output:
        "out/sctype_labelled/{section}_sctyped"
    shell:
        "Rscript src/scripts/sctype.R {input.st} {input.sc} {output}"



rule dimred_qc:
    input:
        "out/{section}_sc_df.csv"
    output:
        "out/{section}_dr_qc3.csv"
    shell:
        "Rscript src/scripts/dimred_qc.R {input} {output}"
