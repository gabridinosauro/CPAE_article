#!/bin/bash -l

# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
### Optional. Set the job name
#SBATCH --job-name=coverm_genomes1
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=%x-%j.out
### REQUIRED. Specify the PI group for this job
#SBATCH --account=laubitz
### Optional. Request email when job begins and ends
### SBATCH --mail-type=ALL
### Optional. Specify email address to use for notification
### SBATCH --mail-user=schiro@email.arizona.edu
### REQUIRED. Set the partition for your job.
#SBATCH --partition=standard
### REQUIRED. Set the number of cores that will be used for this job. 
#SBATCH --ntasks=48
### REQUIRED. Set the number of nodes
#SBATCH --nodes=1
### REQUIRED. Set the memory required for this job.
#SBATCH --mem=240gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=100:00:00


# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------
conda activate coverm

cd /xdisk/pkiela/schiro/CPAE/clean

for file in *_clean_1.fastq; do
	name=${file%%_clean_1.fastq}
	coverm genome --coupled ${name}_clean_1.fastq ${name}_clean_2.fastq --genome-fasta-directory /xdisk/pkiela/schiro/CPAE/bins/ -o ${name}_reads_base.tsv -t 48 --min-read-percent-identity 0.95 --min-read-aligned-length 45 -m reads_per_base --min-covered-fraction 0
done

```
now I create a table with them, in R.
```
### Quantification Bins_cover_M
setwd("/Users/gabri/Library/CloudStorage/Box-Box/CPAE_manuscript/new_analyses/quantification_bins/")
file_list = as.list(list.files())

abundance = data.frame(matrix(nrow = 159, ncol = 0))
for (files in file_list)
{
  temp = read.delim(files,row.names = 1)
  abundance = cbind(abundance, temp)
  stopifnot(rownames(abundance) == rownames(temp))
}

saveRDS(abundance,"coverm_rel_ab.RDS")



```



