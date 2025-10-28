# 1 Prepare files for phylogenetic placement
vcf=
tree=
out_dir=
prepare=$out_dir/tree_data

# pathPhynder v1.2.3
phynder -B -o $out_dir/$prefix.snp $tree $vcf
cd $out_dir
pathPhynder -s prepare -i $tree -p $prefix -f $out_dir/$prefix.snp
cd -


# 2 Phylogenetic placement
fq=
sample=
ref= # $prefix.consensus.fasta
out_dir=
mkdir -p $out_dir
tree=
prepare=

bwa index $ref
bwa aln -n 0.05 -t 1 -l 1000 $ref $fq > $out_dir/$sample.sai
bwa samse $ref $out_dir/$sample.sai $fq | samtools view -bS -F4 | samtools sort -O bam - > $out_dir/$sample.sorted.bam
sambamba markdup -r $out_dir/$sample.sorted.bam $out_dir/$sample.dd.bam
pathPhynder -s all -i $tree -p $prepare -r $ref -b $out_dir/$sample.dd.bam -o $sample 
rm $out_dir/$sample.sai $out_dir/$sample.bam 

