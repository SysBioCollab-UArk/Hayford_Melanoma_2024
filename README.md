## Data repository for &quot;A heterogeneous drug-tolerant persister state is characterized by altered ion channel activity and susceptibility to ferroptosis,&quot; Hayford, Stauffer, Baleami et al. (2024), *[journal] [vol] : [pages]; [DOI: ...](...)*

---

### ***&ast;Instructions for creating panels in all main and supplementary figures based on experimental data in this repository***

- #### Before running any code, open the `Hayford_Melanoma_2024.Rproj` file and activate the virtual environment by running the following code in the R console:
  `library(renv)`
  
  `renv::restore()`
  
  This will install all of the packages in the virtual environment.

  If you get an error about `renv` not being installed, run `install.packages("renv")` to install `renv`. Then run the above commands.

#### To create all of the main and supplementary figures in one script, run ``createPlotsForPaper2024.R`` in the [Hayford_Melanoma_2024](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main) directory. 
  - When the ``Seurat_v5_SKMEL5_combined_hg38.R`` script runs in the [scRNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/code) directory, a prompt for assigning the Seurat clusters will ask you to enter which Seurat clusters are associated with the small idling cluster (IS), large idling cluster (IL), small untreated cluster (UTS), and large untreated cluster (UTL). 
    - Open the ``UMAP_combined_SKMEL5_hg38_qcCCReg_allPlots.pdf`` file in the [scRNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/code) directory. 
    - The small, fast-dividing idling cluster (IS) includes the Seurat cluster with the largest proportion of fast-dividing cells, i.e., the cluster that includes the blue points in the second plot and mostly green points in the third plot. In the example below, this corresponds to Seurat cluster number 5.
    - The large, slow-dividing idling cluster (IL) includes the remaining Seurat clusters in the idling cluster, i.e., the predominantly blue cluster(s) in the second plot. In the example below, this corresponds to Seurat cluster numbers 0 3 6 7.
    - The small, fast-dividing untreated cluster (UTS) includes the Seurat cluster that corresponds to the smaller untreated cluster, i.e., the small red cluster in the second plot. In the example below, this corresponds to Seurat cluster number 8.
    - The large, slow-dividing untreated cluster (UTL) includes the remaining Seurat clusters in the largest untreated cluster, i.e., the large red cluster in the second plot. In the example below, this corresponds to Seurat cluster numbers 1 2 4 9.
    - **Example for defining Seurat clusters**:
      
        - Initial prompt in the R console


          <img width="874" alt="Hayford_Melanoma_2024_SeuratClusterPrompt_Initial" src="https://github.com/user-attachments/assets/3852d939-3901-4ad2-ab77-22a8e576190e">
          

        - `UMAP_combined_SKMEL5_hg38_qcCCReg_allPlots.pdf` for this example
          
        <img width="404" alt="Screenshot_UMAP_combined_SKMEL5_hg38_qcCCReg_allPlots" src="https://github.com/user-attachments/assets/b5fe33f2-6803-47ec-ad3e-dd2461310389">

  
        - Final image of prompt in the R console
          <img width="872" alt="Hayford_Melanoma_2024_SeuratClusterPrompt_Final" src="https://github.com/user-attachments/assets/22beae9f-9e97-41cd-b2e0-011e8fdf372a">

      
  - All panels will be saved to a ``_HAYFORD_2024_PAPER_PLOTS`` directory created by the ``createPlotsForPaper2024.R`` script.
    
#### To create the figures individually, use the information below.

- #### <ins>MAIN FIGURES</ins>

  - #### <ins>FIGURE 1</ins>: 

    **Panels A, D, and E**: In the [scRNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/code) directory, run ``Seurat_v5_SKMEL5_combined_hg38.R``, which pulls data from [scRNA/data/Idling](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/Idling) and [scRNA/data/Untreated](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/Untreated) for the ``barcodes.tsv.gz``, ``features.tsv.gz``, and ``matrix.mtx.gz`` files for idling and untreated conditions, respectively, and lineage information from ``Treated_LineageBC_cellBC.csv`` and``Untreated_LineageBC_cellBC.csv`` in [scRNA/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data).
    
    **Panel B**: In the [RNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/RNA/code) directory, run ``analysis_subclones.R``, which pulls data from ``featureCounts_matrix_all.csv`` in [RNA/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/RNA/data).
    
    **Panel C**: In the [scRNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/code) directory, run ``VISION_SKMEL5.R``, which pulls data from ``combined_includingState.RData`` created from the ``Seurat_v5_SKMEL5_combined_hg38.R`` script in [scRNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/code) and VISION hallmark data from [scRNA/data/VISION_hallmark](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/VISION_hallmark).



  - #### <ins>FIGURE 2</ins>: 

    **Panels A and B**: In the [barcoding/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/barcoding/code) directory, run ``bcSamplingPlot.R``, which pulls data from ``bcPlot.csv`` and ``bcCount_allTechReps.csv`` files in [barcoding/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/barcoding/data).
    
    **Panels C and D**: In the [scRNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/code) directory, run ``Seurat_v5_SKMEL5_combined_hg38.R``, which pulls data from [scRNA/data/Idling](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/Idling) and [scRNA/data/Untreated](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/Untreated) for the ``barcodes.tsv.gz``, ``features.tsv.gz``, and ``matrix.mtx.gz`` files for idling and untreated conditions, respectively, lineage information from ``Treated_LineageBC_cellBC.csv`` and``Untreated_LineageBC_cellBC.csv``, and data for most abundant barcodes from ``top25BCs.RData`` in [scRNA/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data).
    
    **Panel E**: In the [scRNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/code) directory, run ``MetaCluster_bcCycling.R``, which pulls data from ``combined_includingState.RData`` created from the ``Seurat_v5_SKMEL5_combined_hg38.R`` script in [scRNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/code), data for most abundant barcodes from ``top25BCs.RData`` in [scRNA/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data), and ``barcodeFC_forCorPlot.RData`` in [scRNA/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data) for the fold change information of the 25 most abundant barcodes. 
    
    

  - #### <ins>FIGURE 3</ins>: 

    **Panel A**: In the [RNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/RNA/code) directory, run ``analysis_subclones.R``, which pulls data from ``featureCounts_matrix_all.csv`` in [RNA/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/RNA/data).
    
    **Panel B**: In the [RNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/RNA/code) directory, run ``idling_RNAseq_analysis.Rmd``, which pulls data from ``featureCounts_matrix_all.csv`` in [RNA/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/RNA/data).
    
    **Panel C**: In the [ATAC/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ATAC/code) directory, run ``ATAC_analysis_subsample.R``, which pulls data from ``ENCFF356LFX.bed.gz``, ``untreated_all_peaks.narrowPeak_peaks.narrowPeak``, ``trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique_fullClean.bam``, ``idling_peaks.narrowPeak``, and ``idling_subsample_25p_dedup_unique_fullClean.bam`` in [ATAC/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ATAC/data).
    
    **Panel D**: In the [GO_correlations/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/GO_correlations/code) directory, run ``SKMEL5_RNA-ATAC_logq-q.R``, which pulls data from ``untreatedIdling_DEA.RData`` in [GO_correlations/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/GO_correlations/data) and ``SKMEL5_ATACsub25_annotatedPeaks_uniqueShared.RData`` in [ATAC/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ATAC/data).

    **Panel E**: In the [ion_flux/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ion_flux/code) directory, run ``analysisforpaper.Rmd``, which pulls data from ``VQU_008_CPA_3add_1tip_FLUO-8_AM_Feb_12_21_1023.txt`` in [ion_flux/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ion_flux/data) and well conditions from ``platemap_20210211.tsv`` in [ion_flux/platemaps](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ion_flux/platemaps). 

  - #### <ins>FIGURE 4</ins>: 

    **Panel A**: In the [RNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/RNA/code) directory, run ``Ferroptosis_heatmaps.Rmd``, which pulls data from ``merged_countdata.RData``, ``RLD_SC-1,7,10_0,3,8d_20180701.RData``, and ``WP_FERROPTOSIS.v2023.2.Hs.ts`` in [RNA/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/RNA/data).
    
    **Panel B**: In the [dose_response/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/dose_response/code) directory, run ``RSL3_DRC.Rmd``, which pulls data from ``SKMEL5-H2BGFP-postidle_ErastinDR-40uM-2fold_Fer1-0-1uM.txt``, ``SKMEL5-H2BGFP-postidle_RSL3DR-40uM-2fold_Fer1-0-1uM.txt``, ``SKMEL5-H2BGFP-untreated_ErastinDR-40uM-2fold_Fer1-0-1uM.txt``, and ``SKMEL5-H2BGFP-untreated_RSL3DR-40uM-2fold_Fer1-0-1uM.txt`` in [dose_response/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/dose_response/data) and ``ErastinDR plate map.txt`` and ``RSL3DR plate map.txt`` in [dose_response/Platemaps](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/dose_response/Platemaps). 

  - #### <ins>FIGURE 5</ins>: 

    **Panels A and B**: These panels were made using BioRender.
    

- #### <ins>SUPPLEMENTARY FIGURES</ins>

  - #### <ins>FIGURE S1</ins>:

    **Panel A and B**: In the [scRNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/code) directory, run ``Seurat_v5_SKMEL5_combined_hg38.R``, which pulls data from [scRNA/data/Idling](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/Idling) and [scRNA/data/Untreated](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/Untreated) for the ``barcodes.tsv.gz``, ``features.tsv.gz``, and ``matrix.mtx.gz`` files for idling and untreated conditions, respectively, and lineage information from ``Treated_LineageBC_cellBC.csv`` and``Untreated_LineageBC_cellBC.csv`` in [scRNA/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data).

  - #### <ins>FIGURE S2</ins>:

    **Panels A, B, and C**: In the [barcoding/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/barcoding/code) directory, run ``bcSamplingPlot.R``, which pulls data from ``bcPlot.csv`` and ``bcCount_allTechReps.csv`` files in [barcoding/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/barcoding/data).
    
    **Panel D**: In the [scRNA/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/code) directory, run ``Seurat_v5_SKMEL5_combined_hg38.R``, which pulls data from [scRNA/data/Idling](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/Idling) and [scRNA/data/Untreated](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data/Untreated) for the ``barcodes.tsv.gz``, ``features.tsv.gz``, and ``matrix.mtx.gz`` files for idling and untreated conditions, respectively, lineage information from ``Treated_LineageBC_cellBC.csv`` and``Untreated_LineageBC_cellBC.csv``, and data for most abundant barcodes from ``top25BCs.RData`` in [scRNA/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/scRNA/data).

  - #### <ins>FIGURE S3</ins>:

    **Panels A and B**: In the [ATAC/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ATAC/code) directory, run ``preseq_subsample.R``, which pulls data from ``trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_cCurveResults.txt``, ``trimmed_3334-CH-2-CGTACTAG-CTCTCTAT_S122_aligned_sorted_cCurveResults.txt``, ``idling_subsample_25p_dedup_cCurveResults.txt``, ``idling_subsample_33p_dedup_cCurveResults.txt``,  ``idling_subsample_50p_dedup_cCurveResults.txt``, ``untreated_dedup_ISM.txt``,  ``idling_dedup_ISM.txt``, ``idling_subsample_25p_dedup_ISM.txt``, ``idling_subsample_33p_dedup_ISM.txt``, and ``idling_subsample_50p_dedup_ISM.txt`` in [ATAC/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ATAC/data).
    
    **Panels C, D, and E**: In the [ATAC/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ATAC/code) directory, run ``ATAC_analysis_subsample.R``, which pulls data from ``ENCFF356LFX.bed.gz``, ``untreated_all_peaks.narrowPeak_peaks.narrowPeak``, ``trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique_fullClean.bam``, ``idling_peaks.narrowPeak``, and ``idling_subsample_25p_dedup_unique_fullClean.bam`` in [ATAC/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ATAC/data) and ``idling_peaks.narrowPeak`` and ``untreated_all.narrowPeak`` in [ATAC/data/venn](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ATAC/data/venn).
    
    **Panels F and G**: In the [GO_correlations/code](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/GO_correlations/code) directory, run ``SKMEL5_RNA-ATAC_logq-q.R``, which pulls data from ``untreatedIdling_DEA.RData`` in [GO_correlations/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/GO_correlations/data) and ``SKMEL5_ATACsub25_annotatedPeaks_uniqueShared.RData`` in [ATAC/data](https://github.com/SysBioCollab-UArk/Hayford_Melanoma_2024/tree/main/ATAC/data).

  - #### <ins>FIGURE S4</ins>:

    Figure created using WikiPathways. 
