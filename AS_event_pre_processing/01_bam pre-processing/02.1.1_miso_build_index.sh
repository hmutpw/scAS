#!/bin/bash
# preparing for MISO gff3 file.
# Stpe 1 Creating custom GFF annotations for MISO.
# first of all you should supply your gtf file for mapping before, 
# eg. 'gencode.vM22.annotation.gff' is what I used before for Tophat mapping.
# Please creating a 'ensGene.txt' file with 'gtfToGenePred' (UCSC Tools from
# http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/) from your '.gtf' file,
# please follow the instruction described in
# 'http://rnaseqlib.readthedocs.io/en/clip/#creating-custom-gff-annotations-for-miso'
conda activate miso-env
ref_dir="/mnt/odisk/raw4T/HSC_Autophagy/ref_genome"
out_index="/mnt/odisk/raw4T/HSC_Autophagy/STAR_res/MISO/preFileForMISO"
mkdir -p $out_index/gff;
echo "Generating ensGene.txt with gtfToGenePred";
/home/hmutpw/Applications/UCSCTools/gtfToGenePred -genePredExt -includeVersion -ignoreGroupsWithoutExons $ref_dir/gencode.vM22.primary_assembly.annotation.gtf $out_index/gff/gencGene.txt
sed "s/^/1\t/g" $out_index/gff/gencGene.txt > $out_index/gff/ensGene.txt
rm $out_index/gff/gencGene.txt
python /home/hmutpw/Applications/rnaseqlib-clip/rnaseqlib/gff/gff_make_annotation.py $out_index/gff/ $out_index/gff/ --flanking-rule commonshortest --genome-label vM22

# Stpe 2 Creating GFF Index for MISO.
mkdir -p $out_index/gffIndex;
echo "Generating gffIndex with index_gff";
index_gff --index $out_index/gff/commonshortest/A3SS.vM22.gff3 $out_index/gffIndex/A3SS
index_gff --index $out_index/gff/commonshortest/A5SS.vM22.gff3 $out_index/gffIndex/A5SS
index_gff --index $out_index/gff/commonshortest/RI.vM22.gff3 $out_index/gffIndex/RI
index_gff --index $out_index/gff/commonshortest/MXE.vM22.gff3 $out_index/gffIndex/MXE
index_gff --index $out_index/gff/commonshortest/SE.vM22.gff3 $out_index/gffIndex/SE

# step3 constitutive exons
mkdir -p $out_index/estimateInsertSize
exon_utils --get-const-exons $ref_dir/gencode.vM22.primary_assembly.annotation.gff3 --min-exon-size 1000 --output-dir $out_index/estimateInsertSize/
