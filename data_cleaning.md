### First we quality trimmed reads
```
#!/bin/bash -l

# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
### Optional. Set the job name
#SBATCH --job-name=trimming
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
#SBATCH --ntasks=20
### REQUIRED. Set the number of nodes
#SBATCH --nodes=1
### REQUIRED. Set the memory required for this job.
#SBATCH --mem=80gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=05:00:00


# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------
conda activate quality_check

cd /xdisk/laubitz/schiro/CPAE/concat

mkdir trimmed
mkdir unpaired 
mkdir trimmed/logs

for file in *1.f*; do
	name=${file%%_R1.*}

trimmomatic PE ${name}_R1.fastq.gz ${name}_R2.fastq.gz \
trimmed/${name}_R1_trm.fastq.gz unpaired/${name}_R1_unpaired.fastq.gz \
trimmed/${name}_R2_trm.fastq.gz unpaired/${name}_R2_unpaired.fastq.gz \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:50  -threads 20  &>> trimmed/logs/${name}_log.txt

done



```
### Then we remove host contamination
```
#!/bin/bash -l

# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
### Optional. Set the job name
#SBATCH --job-name=Human_rem
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=%x-%j.out
### REQUIRED. Specify the PI group for this job
#SBATCH --account=laubitz
### Optional. Request email when job begins and ends
### SBATCH --mail-type=ALL
### Optional. Specify email address to use for notification
### SBATCH --mail-user=schiro@arizona.edu
### REQUIRED. Set the partition for your job.
#SBATCH --partition=standard
### REQUIRED. Set the number of cores that will be used for this job. 
#SBATCH --ntasks=48
### REQUIRED. Set the number of nodes
#SBATCH --nodes=1
### REQUIRED. Set the memory required for this job.
#SBATCH --mem=200gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=50:00:00


# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------
conda activate bowtie

cd /xdisk/laubitz/schiro/CPAE/concat/trimmed
mkdir clean

for file in *R1_trm.fastq.gz*; do
	name=${file%%_R1_trm.fastq.gz}
	bowtie2 -p 94 -x /groups/laubitz/databases/reference_genomes/GRCh38_noalt_as/GRCh38_noalt_as -1 ${name}_R1_trm.fastq.gz -2 ${name}_R2_trm.fastq.gz 	--un-conc-gz clean/${name}_clean  > ${name}_mapped_and_unmapped.sam
done

cd /xdisk/laubitz/schiro/CPAE/concat/trimmed/clean
conda deactivate
conda activate multiqc

mkdir fastqc
fastqc *.gz* -t 20 -o fastqc

cd fastqc
multiqc .

done
```




