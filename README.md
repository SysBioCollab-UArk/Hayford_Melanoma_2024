## Data repository for &quot;A heterogeneous drug-tolerant persister state is characterized by altered ion channel activity and susceptibility to ferroptosis,&quot; Hayford, Stauffer, Baleami et al. (2024), *[journal] [vol] : [pages]; [DOI: ...](...)*

---

### ***&ast;Instructions for creating panels in all main and supplementary figures based on experimental data in this repository***

- #### <ins>MAIN FIGURES</ins>

  - #### <ins>FIGURE 1</ins>: 

    **Panel A**: In the [scRNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/code) directory, run ``Seurat_v5_SKMEL5_combined_hg38.R``, which pulls data from [scRNA/data/Idling](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/Idling) and [scRNA/data/Untreated](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/Untreated) for the ``barcodes.tsv.gz``, ``features.tsv.gz``, and ``matrix.mtx.gz`` files for idling and untreated conditions, respectively, and lineage information from ``Treated_LineageBC_cellBC.csv`` and``Untreated_LineageBC_cellBC.csv`` in [scRNA/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data).
    
    **Panel B**: In the [RNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/RNA/code) directory, run ``analysis_subclones.R``, which pulls data from ``featureCounts_matrix_all.csv`` in [RNA/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/RNA/data).
    
    **Panel C**: In the [scRNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/code) directory, run ``VISION_SKMEL5.R``, which pulls data from ``combined_includingState.RData`` on a local drive and VISION hallmark data from [scRNA/data/VISION_hallmark](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/VISION_hallmark).
    
    **Panel D**: In the [scRNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/code) directory, run ``Seurat_v5_SKMEL5_combined_hg38.R``, which pulls data from [scRNA/data/Idling](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/Idling) and [scRNA/data/Untreated](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/Untreated) for the ``barcodes.tsv.gz``, ``features.tsv.gz``, and ``matrix.mtx.gz`` files for idling and untreated conditions, respectively, and lineage information from ``Treated_LineageBC_cellBC.csv`` and``Untreated_LineageBC_cellBC.csv`` in [scRNA/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data).
    
    **Panel E**: In the [scRNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/code) directory, run ``Seurat_v5_SKMEL5_combined_hg38.R``, which pulls data from [scRNA/data/Idling](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/Idling) and [scRNA/data/Untreated](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/Untreated) for the ``barcodes.tsv.gz``, ``features.tsv.gz``, and ``matrix.mtx.gz`` files for idling and untreated conditions, respectively, and lineage information from ``Treated_LineageBC_cellBC.csv`` and``Untreated_LineageBC_cellBC.csv`` in [scRNA/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data).

  - #### <ins>FIGURE 2</ins>: 

    **Panel A**: In the [barcoding/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/barcoding/code) directory, run ``bcSamplingPlot.R``, which pulls data from ``bcPlot.csv`` and ``bcCount_allTechReps.csv`` files in [barcoding/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/barcoding/data).
    
    **Panel B**: In the [barcoding/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/barcoding/code) directory, run ``bcSamplingPlot.R``, which pulls data from ``bcPlot.csv`` and ``bcCount_allTechReps.csv`` files in [barcoding/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/barcoding/data).
    
    

  - #### <ins>FIGURE 3</ins>: 

    xxx

    **Panel E**: In the [ion_flux/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ion_flux/code) directory, run ``analysisforpaper.Rmd``, which pulls data from ``VQU_008_CPA_3add_1tip_FLUO-8_AM_Feb_12_21_1023.txt`` in [ion_flux/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ion_flux/data) and well conditions from ``platemap_20210211.tsv`` in [ion_flux/platemaps](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ion_flux/platemaps). 

  - #### <ins>FIGURE 4</ins>: 

    xxx

  - #### <ins>FIGURE 5</ins>: 

    xxx

- #### <ins>SUPPLEMENTARY FIGURES</ins>

  - #### <ins>FIGURE S1</ins>:

    xxx

  - #### <ins>FIGURE S2</ins>:

    xxx

  - #### <ins>FIGURE S3</ins>:

    xxx

  - #### <ins>FIGURE S4</ins>:

    xxx
