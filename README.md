# Adaptation of <i>Ustilago maydis</i> SG200 to hydrogen peroxide (H<sub>2</sub>O<sub>2</sub>)

## Repository created by Jorge Cuamatzi-Flores

Repository with scripts to reproduce the analysis and figures presented in: [insert link to paper here]


## Bioinformatic analysis of sequenced genomes

In this experiment we sequenced by short-reads the genome of <i>U. maydis</i> SG200 and UmH<sub>2</sub>O<sub>2</sub>-R

Pipeline to process the `fastq` files:

The raw fastq of the sequenced genomes of U. maydis SG200 and the adapted H<sub>2</sub>O<sub>2</sub> strain are accessible on NCBI with accession numbers SRR25650020 and SRR25650011, respectively.

The document `Data_Sheet_to_Download_Fastq_Files.csv` in the main directory `~/Umaydis_experimental_Evolution` contains the accession number assigned to each sample.


### Cleaning

See the `01_Cleaning` subdirectory

### Mapping and coverage plotting

See the `02_MappingAndCoveragePlotting` subdirectory

The <b>Figure 3</b> of the main manuscript can be found in `02_MappingAndCoveragePlotting/Figures/Plot_USMA_521_v2_9.Log2Ratio.png`

### SNP Calling

See the `03_SNP_Calling` subdirectory

## Analysis of collected data

Analysis of data collected from experimental evolution, analysis of gene expression of UMAG_11067 and phenotypic analysis of oexUMAG_11067 (strain that overexpress UMAG_11067 gene)

### Manuscript Figures

See the subdirectory `04_USMA_Data_Growth_qPCR_Infection`

### Supplementary Material

See the `03_SNP_Calling` subdirectory




