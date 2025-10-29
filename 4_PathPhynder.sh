# 1 Prepare files for phylogenetic placement
vcf= # filtered bi-allelic vcf.gz file 
tree= # ultrametric tree, nwk format
phynder_dir= # output folder for phynder


# pathPhynder v1.2.3
phynder -B -o $phynder_dir/$prefix.snp $tree $vcf
cd $phynder_dir
pathPhynder -s prepare -i $tree -p $prefix -f $phynder_dir/$prefix.snp
cd -


# 2 Phylogenetic placement
fq= # $sample.enrich.fq
sample= # sample ID
ref= # $prefix.consensus.fasta
out_dir= # output folder
mkdir -p $out_dir
tree= # ultrametric tree, nwk format
prepare=$phynder_dir/tree_data

bwa index $ref
bwa aln -n 0.05 -t 1 -l 1000 $ref $fq > $out_dir/$sample.sai
bwa samse $ref $out_dir/$sample.sai $fq | samtools view -bS -F4 | samtools sort -O bam - > $out_dir/$sample.sorted.bam
sambamba markdup -r $out_dir/$sample.sorted.bam $out_dir/$sample.dd.bam
pathPhynder -s all -i $tree -p $prepare -r $ref -b $out_dir/$sample.dd.bam -o $sample 
rm $out_dir/$sample.sai $out_dir/$sample.bam 

