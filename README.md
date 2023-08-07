Snakemake workflow for segmentation and primary analysis of imaging-based spatial transcriptomic data. 

Raw data used in this project (src directory), and relevant outputs (out directory) can be accessed here: https://drive.google.com/drive/folders/1Tqv58Bto6mCE_oyBI9vyJodWHe_qBwh1?usp=sharing

To run workflow, Baysor binary needs to be downloaded into main directory. Details here: https://github.com/kharchenkolab/Baysor#binary-download

Segmentation of image data is performed using Cellpose followed by Baysor. Baysor parameter configuration is in Drive folder linked above as baysor_configs.toml

After segmentation, various dimension reduction and cell type labelling methods are tested. 
