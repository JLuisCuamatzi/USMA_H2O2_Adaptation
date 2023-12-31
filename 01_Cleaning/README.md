## 01.- Cleaning of genomes from sequenced colonies


The cleaning of raw `fastq` files was done with [<b>fastp/0.20.0</b>](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234) using the next code:


```
## Example for sample: 2021EE01

# Inputs:
input_fastqR1="2021EE01_R1.fq.gz"
input_fastqR2="2021EE01_R2.fq.gz"

# Outputs:
output_fastqR1="2021EE01_R1_clean.fastq.gz"
output_fastqR2="2021EE01_R2_clean.fastq.gz"

# Unpaired fastq files
unpaired1="2021EE01_unpaired_R1_clean.fastq.gz"
unpaired2="2021EE01_unpaired_R2_clean.fastq.gz"

# json file
json_file="2021EE01_fastq.json"

# html file (Quality Control Report)
html_file="2021EE01_fastp.html"

# Running fastp
fastp -i ${input_fastqR1} -o ${output_fastqR1}           \
      -I ${input_fastqR2} -O ${output_fastqR2}           \
      --unpaired1 ${unpaired1} --unpaired2 ${unpaired2}  \
      -w 2 -y -x -z 9 -j ${json_file} -h ${html_file}

```

To optimize the execution of this process, I wrote the `01_Cleaning.py` python script that generates sh files with the aforementioned commands.

This script needs a sample sheet file in `csv` format where the ID of each sample is listed. 
 - The sample sheet is in the main directory: `~/Umaydis_experimental_Evolution/USMA_EE_Colonies_SampleSheet.csv`.
 - `python3` is required to execute this script.

Execute `01_Cleaning.py` script
```
python3 01_Cleaning.py -c ../USMA_EE_Colonies_SampleSheet.csv

```

Once this script is executed, a directory named `shFiles` will be created. In this directory will be the `.sh` files to be executed.

You can run these `.sh` files with the next command:

```
# 

sh shFiles/2021EE01_Cleaning.sh

sh shFiles/2021EE18_Cleaning.sh

```
