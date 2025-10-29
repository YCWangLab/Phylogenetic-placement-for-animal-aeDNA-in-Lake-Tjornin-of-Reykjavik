# 1 Call assembly from VCF 
# # i.e., replace the reference with alternative SNPs and InDels for each samples
vcf= # the *raw.vcf.gz from 1_Reference_panel_construction.sh - Step. 4
ref= # consensus genome 
out_dir= # output folder
prefix= # job prefix
mkdir -p $out_dir
temp=$out_dir/temp
mkdir -p $temp

# bcftools v1.22
bcftools query -l $vcf > $temp/$prefix.sample.txt
while read sample
do
bcftools consensus -s $sample -f $ref $vcf | seqtk seq -A - | sed "s/^>.*$/>$sample/" > $temp/${sample}.fasta
done < $temp/$prefix.sample.txt
cat $temp/*.fasta > $out_dir/$prefix.vcf_consensus.fasta
rm -rf $temp


# 2 Construct ultramatric tree with BEAST and get the consensus tree

# # 2.1 alignment with mafft
input_fasta=$prefix.vcf_consensus.fasta
thread=
mafft --auto --adjustdirection --thread $thread $input_fasta > $prefix.aligned
trimal -in $prefix.aligned -out $prefix.aligned.trimmed.aln -automated1 $trimal_opt
rm $prefix.aligned

# # 2.2 make beast template xml
template=ultrametric_tree_template.xml 
chain_length=1000000000 # make sure the ESS for each parameter higher than 200
store_every=5000
prefix=
python make_beast_xml.py --xml_template $template --fasta_file $prefix.aligned.trimmed.aln --output_file $prefix.ultrametricTree.xml --prefix $prefix --chain_length $chain_length --store_every $store_every

# # 2.3 run beast
# beast v2.7.7
out_dir=
prefix=
xml= # xml file from step 2.2
beast -beagle -threads $thread  -prefix $out_dir/$prefix -overwrite $xml
# # Get consensus tree using logcombiner and treeannotator if multiple beast were run
logcombiner -b 10 -log run1.log -log run2.log -log run3.log -o $prefix.combined.log
logcombiner -b 10 -log run1.trees -log run2.trees -log run3.trees -o $prefix.combined.trees
treeannotator -burnin 10 -heights mean $prefix.combined.trees $prefix.combined.nex


# 2.4 Convert nexus format to newick format and remove posterior
nex=$prefix.combined.nex
nwk=$prefix.combined.tree
python nex_to_nwk.py -i $nex -o $nwk 
