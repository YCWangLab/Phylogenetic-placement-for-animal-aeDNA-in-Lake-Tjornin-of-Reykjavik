# 1. Call consensus genomes
multi_fasta= # multi-fasta file of population mitogenomes
prefix= # job prefix
thread=
# MAFFT v7.490
mafft --auto --adjustdirection --thread $thread $multi_fasta > $prefix.aligned
# EMBOSS v6.6.0.0
cons -sequence $prefix.aligned -outseq $prefix.consensus.temp -identity 0 -plurality 0 -name ${prefix}_consensus -auto
sed '/^>/!s/[a-z]/\U&/g' $prefix.consensus.temp | sed 's/-/N/g' > $prefix.consensus.fasta
rm $prefix.consensus.temp
bwa index $prefix.consensus.fasta


# 2. Get error-free pair-end reads from haploid genomes
multi_fasta= # multi-fasta file of population mitogenomes
out_dir= # output folder

temp_dir=$out_dir/temp
mkdir -p $temp_dir
# seqkit v2.10.0
seqkit split -i -O $temp_dir $multi_fasta
remove_base="$(basename "$multi_fasta" .fasta).part_"
for fasta in $temp_dir/*fasta; do
    header=$(basename "$fasta" .fasta)
    new_header=$(echo $header | sed s/^${remove_base}//)
    mv $fasta $temp_dir/$new_header.fasta
done
# ART v2.5.8
for fasta in $temp_dir/*fasta; do
    header=$(basename "$fasta" .fasta)
    art_illumina -ss HS25 -ef -p -l 150 -f 50 -m 500 -s 10  -i $fasta -o ${out_dir}/${header}_
done

rm -rf $temp_dir


# 3. Map reads to consensus genome using BWA MEM 
fq_path= # out_dir of Step. 2
bam_path= # output folder for bam files
mkdir -p $bam_path
ref= # consensus genome, i.e., $prefix.consensus.fasta
# bwa v0.7.19-r1273
bwa index $ref

for file in $fq_path/*_1.fq
    do
    sample=$(basename $file _1.fq)
    echo "Processing $sample"
    # samtools v1.22.1
	bwa mem -t 16 -M -R "@RG\tID:${sample}\tSM:${sample}" $ref $fq_path/${sample}_1.fq $fq_path/${sample}_2.fq  | samtools view -F 4 -q 30 -bS - | samtools sort -O bam -o $out_dir/$sample.sorted.bam -T $out_dir/$sample.sorted.temp
	# sambamba v1.0.1
    sambamba markdup -r $out_dir/$sample.sorted.bam $out_dir/$sample.dd.bam
	samtools index $out_dir/$sample.dd.bam
	rm $out_dir/$sample.sorted.bam $out_dir/$sample.dd.metrics
done


# 4. SNP calling with short reads using GATK HaplotypeCaller
ref= # consensus genome, i.e., $prefix.consensus.fasta
bam_path= # folder for bam files
suffix=.dd.bam
gvcf_path= # output folder for g.vcf files
mkdir -p $gvcf_path
samtools faidx $ref 
# gatk v4.3.0.0
gatk CreateSequenceDictionary -R $ref -O ${ref%.*}.dict

for bam in $bam_path/*$suffix
    do
    sample=$(basename $bam $suffix)
    echo "Processing $sample"
	gatk HaplotypeCaller -R $ref -I $bam -O $gvcf_path/$sample.g.vcf --emit-ref-confidence GVCF --sample-ploidy 1 --kmer-size 15 --kmer-size 25
done

prefix=
out_dir=
gvcf_files=$(ls $gvcf_path/*.g.vcf | sed 's/^/-V /' | tr '\n' ' ') 
gatk CombineGVCFs -R $ref $gvcf_files -O $out_dir/$prefix.combined.g.vcf
gatk GenotypeGVCFs -R $ref -V $out_dir/$prefix.combined.g.vcf -O $out_dir/$prefix.raw.vcf 
# bcftools v1.22
bcftools view -v snps $out_dir/$prefix.raw.vcf -Oz -o $out_dir/$prefix.raw.vcf.gz
bcftools index $out_dir/$prefix.raw.vcf.gz
bcftools view -m2 -M2 -v snps -i 'INFO/DP>10 & QUAL>30' $out_dir/$prefix.raw.vcf -o $out_dir/$prefix.biallelic.vcf
bcftools view -v snps $out_dir/$prefix.biallelic.vcf -Oz -o $out_dir/$prefix.biallelic.vcf.gz
bcftools index $out_dir/$prefix.biallelic.vcf.gz



