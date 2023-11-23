# Adaptation of <i>Ustilago maydis</i> SG200 to hydrogen peroxide (H<sub>2</sub>O<sub>2</sub>)

## Repository created by Jorge Cuamatzi-Flores

Repository with scripts to reproduce the analysis and figures presented in: [Enhanced Oxidative Stress Resistance in Ustilago maydis and Its Implications on the Virulence](https://www.researchgate.net/publication/375481495_Enhanced_Oxidative_Stress_Resistance_in_Ustilago_maydis_and_Its_Implications_on_the_Virulence#fullTextFileContent)


## Bioinformatic analysis of sequenced genomes

In this experiment we sequenced by short-reads the genome of <i>U. maydis</i> SG200 and UmH<sub>2</sub>O<sub>2</sub>-R

Pipeline to process the `fastq` files:

The raw fastq of the sequenced genomes of U. maydis SG200 and the adapted H<sub>2</sub>O<sub>2</sub> strain are accessible on NCBI with accession numbers SRR25650020 and SRR25650011, respectively.

The document `Data_Sheet_to_Download_Fastq_Files.csv` in the main directory `~/Umaydis_experimental_Evolution` contains the accession number assigned to each sample.


### Cleaning

See the `01_Cleaning` subdirectory

### Mapping and coverage plotting

See the `02_MappingAndCoveragePlotting` subdirectory

### SNP Calling

See the `03_SNP_Calling` subdirectory

## Analysis of collected data

Analysis of data collected from experimental evolution, analysis of gene expression of UMAG_11067 and phenotypic analysis of oexUMAG_11067 (strain that overexpress UMAG_11067 gene)

### Manuscript Figures

See the subdirectory `04_USMA_Data_Growth_qPCR_Infection`

### Supplementary Material

See the `03_SNP_Calling` subdirectory




