1st step assembly reads
```
conda activate mega_hit

cd /xdisk/barberan/schiro/Metagen_Daniel/trimmed/human_clean


mkdir assemblies
for file in *_R1.fastq.gz; do
	name=${file%%_R1.fastq.gz}
	megahit -1 ${name}_R1.fastq.gz -2 ${name}_R2.fastq.gz \
        -o assemblies/${name} --out-prefix ${name} --k-min 27 -t 94 

rm -r assemblies/${name}/intermediate_contigs
done
```
Second step, run prodigal to identify ORF

```
conda activate quast

cd  /xdisk/barberan/schiro/Metagen_Daniel/contigs

mkdir filt_500
for file in *.contigs.fa*; do
	name=${file%%.contigs.fa*}
	reformat.sh in=${file} out=filt_500/${name}_con500.fa minlength=500
done

conda deactivate
conda activate prokka


cd /xdisk/barberan/schiro/Metagen_Daniel/contigs/filt_500


mkdir prodigal
for file in *.f*; do
	name=${file%%_con500.fa}
	prodigal -i ${file} -a prodigal/${name}_pro.fasta -d prodigal/${name}_gene.fasta \
-f gff \
-p meta \
-q \
-m

done
```
3rd step, clustering

```
conda activate quast

cd  /xdisk/barberan/schiro/Metagen_Daniel/contigs/filt_500/prodigal
## cat gene.fasta > all_genes.fasta

reformat.sh in=all_genes.fasta out=all_gene_rm.fasta minlength=100


conda deactivate
conda activate MMseq

mmseqs createdb all_gene_rm.fasta DB
mmseqs cluster DB all_gene_cluster all_gene_cluster_tmp --min-seq-id 0.95 -c 0.90 --cov-mode 1 --cluster-mode 2
mmseqs result2repseq DB all_gene_cluster all_gene_cluster_rep
mmseqs result2flat DB DB all_gene_cluster_rep all_gene_cluster_rep.fasta --use-fasta-header


```
Now convert the genes into AA into protein and annotate them with KOfam scan

```
conda activate KOfam

cd /xdisk/laubitz/schiro/CPAE/gene_catalog
exec_annotation -f mapper -o annotations.txt protein_catalog.fasta 

done
```

Quantification

```
conda activate  mapping



cd /xdisk/barberan/schiro/Metagen_Daniel/contigs/filt_500/prodigal/renamed
mkdir clustered
cat *gene.fasta > clustered/all_gene_cluster_rep.fasta
cd clustered
bwa index all_gene_cluster_rep.fasta -p all_gene_cluster_bwa.gene


cd /xdisk/barberan/schiro/Metagen_Daniel/trimmed/human_clean
mkdir sam_files_genes

for file in *_R1*; do
	name=${file%%_R1*}
	bwa mem /xdisk/barberan/schiro/Metagen_Daniel/contigs/filt_500/prodigal/renamed/clustered/all_gene_cluster_bwa.gene ${name}_R1.fastq.gz 	${name}_R2.fastq.gz -t 48 > sam_files_genes/${name}_gene.sam
done

cd sam_files_genes
for file in *.sam; do
	name=${file%%.sam}
	samtools view -@ 48 -bS ${file} > ${name}.bam
	samtools sort -m 4G -@ 48 -o ${name}.bam.sorted ${name}.bam
	rm ${file}
done



cd /xdisk/barberan/schiro/Metagen_Daniel/trimmed/human_clean/sam_files_genes

mkdir final_counts
for file in *gene.bam.sorted; do
	name=${file%%_gene.bam.sorted}
	bamToBed -i ${file} > ${name}.bed
	perl /xdisk/barberan/schiro/Metagen_Daniel/Commands/compute_abundance_from_bed.pl ${name}.bed /xdisk/barberan/schiro/Metagen_Daniel/contigs/filt_500/prodigal/renamed/clustered/all_gene_cluster_rep.fasta final_counts/${name}_gene_bwa_count.txt
done



conda activate coverm

cd /xdisk/pkiela/schiro/Kathi/clean/sam_files_genes  

for file in *_gene.bam.sorted; do
	name=${file%%_gene.bam.sorted}
	coverm contig -b $file -t 24 -m count --min-read-percent-identity 0.95 --min-read-aligned-length 45 -o ${name}_gene_counts.txt
	coverm contig -b $file -t 24 -m mean --min-read-percent-identity 0.95 --min-read-aligned-length 45 -o ${name}_gene_mean.txt
	coverm contig -b $file -t 24 -m trimmed_mean --min-read-percent-identity 0.95 --min-read-aligned-length 45 -o ${name}_gene_trim_mean.txt
	coverm contig -b $file -t 24 -m rpkm --min-read-percent-identity 0.95 --min-read-aligned-length 45 -o ${name}_gene_rpkm.txt
done


```




