# 1. Enrich mitochondrial reads of target species
fq_path= # aeDNA fastq from different ages
ref= # $prefix.consensus.fasta for each species
suffix=
thread=
out_dir=

for fq in $fq_path/*$suffix
do
    sample=$(basename $fq | sed s/$suffix\$//)
    bwa aln -n 0.1 -t $thread -l 1000 $ref $fq > $out_dir/$sample.sai && \
    bwa samse $ref $out_dir/$sample.sai $fq | samtools view -F 4 -bS - | \
    samtools sort -O bam -o $out_dir/$sample.enrich.bam -T $out_dir/$sample.sorted.temp && \
    rm $out_dir/$sample.sai && \
    samtools fastq $out_dir/$sample.enrich.bam > $out_dir/$sample.enrich.fastq && \
    seqtk trimfq -b 3 -e 3 $out_dir/$sample.enrich.fastq > $out_dir/$sample.trimmed.fastq && \
    rm $out_dir/$sample.enrich.bam
done

# 2. Remove noise
# # 2.1 Competitive mapping
fq=
out_dir=
thread=
DB=
names_dmp=
nodes_dmp=
acc2tax=
taxa=
output_base=
ngsLCA_py=

bam_path=$out_dir/${output_base}_bam
mkdir -p $bam_path

cat $DB | while read line
do
    bDB=$(basename "$line" | sed 's/\//_/g')
    echo "Mapping $fq against $bDB"
    echo "bowtie2 output: $bam_path/$bDB.bam"
    # bowtie2 v2.5.4
    bowtie2 -k 1000 -x $line -U $fq --threads $thread --no-unal | samtools view -@ $thread -b > $bam_path/$bDB.bam
done

# ngsLCA v0.9
python $ngsLCA_py -b $bam_path --thread $thread --names $names_dmp --nodes $nodes_dmp -a $acc2tax --minedit 0 --maxedit 0 -o $out_dir/$output_base

grep -v $taxa $out_dir/$output_base.lca | cut -d $'\t' -f1 | cut -d ":" -f1- | rev | cut -d ":" -f4- | rev > $out_dir/$output_base.remove_id.txt

rm -rf $out_dir/$output_base.lca $bam_path

# # 2.2 Get unclassified or reads classified to target taxa
out_dir=
for file in $fq_path/*fastq
do
sample=`basename $file | cut -d. -f1 `
seqkit grep -v -f $output_base.remove_id.txt $file > $out_dir/$sample.cleaned.fq
done