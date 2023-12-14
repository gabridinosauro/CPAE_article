First part of binning was running a co-assembly
```
conda activate metawrap-env
cd /home/u1/schiro/metaWRAP
export PATH=$PWD/bin:$PATH
cd /xdisk/laubitz/schiro/CPAE/clean/concat_clean_files
metawrap assembly -1 ALL_READS_1.fastq.gz -2 ALL_READS_2.fastq.gz -m 192 -t 48  -o ASSEMBLY

```


The we need to run a binning 
```
cd /xdisk/laubitz/schiro/CPAE/clean/
metawrap binning -o INITIAL_BINNING -t 48 -a concat_clean_files/ASSEMBLY/final_assembly.fasta  --concoct *fastq
```


and refinement
```
cd /xdisk/laubitz/schiro/CPAE/clean/
metawrap bin_refinement -o BIN_REFINEMENT -t 48 -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ 
```
and reassembly
```
metawrap reassemble_bins -o BIN_REASSEMBLY -1 concat_clean_files/ALL_READS_1.fastq.gz -2 concat_clean_files/ALL_READS_2.fastq.gz -t 48 -m 192 -c 70 -x 10 -b BIN_REFINEMENT/metawrap_70_10_bins
```

Now we can assign taxonomy with GTDBTk
```
conda activate /contrib/laubitz/GTDB
gtdbtk classify_wf --genome_dir /xdisk/laubitz/schiro/CPAE/clean/BIN_REASSEMBLY/reassembled_bins --out_dir /xdisk/laubitz/schiro/CPAE/clean/BIN_REASSEMBLY/reassembled_bins/classification --cpus 48 --extension fa
```


Do quantification with coverM

```
conda activate coverm

cd /xdisk/pkiela/schiro/CPAE/clean

for file in *_clean_1.fastq; do
	name=${file%%_clean_1.fastq}
	coverm genome --coupled ${name}_clean_1.fastq ${name}_clean_2.fastq --genome-fasta-directory /xdisk/pkiela/schiro/CPAE/bins/ -o ${name}_reads_base.tsv -t 48 --min-read-percent-identity 0.95 --min-read-aligned-length 45 -m reads_per_base --min-covered-fraction 0
done
```

Following up, we need to functionally annotate these genomes

```
conda activate prokka

cd /xdisk/pkiela/schiro/CPAE/bins
mkdir prodigal
for file in *.f*; do
	name=${file%%.fna}
	prodigal -i ${file} -a prodigal/${name}_pro.fasta -d prodigal/${name}_gene.fasta \
-f gff \
-p meta \
-q \
-m

done

```
And now run KOfam scan

```
conda activate KOfam


cd /xdisk/pkiela/schiro/CPAE/bins/prodigal/amino
for file in *_pro.fasta; do
	name=${file%%_pro.fasta}
	exec_annotation -f mapper -o ${name}_annotations.txt --cpu 20 $file
done



```









