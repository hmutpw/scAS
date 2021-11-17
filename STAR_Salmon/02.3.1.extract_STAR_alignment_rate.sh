#!/bin/bash
main_path="/mnt/odisk/raw4T/HSC_Autophagy"
for sample in `cat $main_path/sample_infor/total_sample_list.txt`
do
in_dir="$main_path/STAR_res/alignment_res/$sample"
out_dir="$main_path/STAR_res/alignment_summary/$sample"
mkdir -p $out_dir;
echo '[extract read alignment rate from STAR] '
if [ ! -f "$in_dir/Log.final.out" ];then
	echo "Log.final.out for sample $sample does not exists!";
	else
		echo "summarizing alignment result for sample $sample";
		header="item\tcount";
		totalpaired=`cat $in_dir/Log.final.out | grep "Number of input reads |" | sed 's/^[ \t]\+Number of input reads |[\t]\+\([0-9]\+\)$/total_reads\t\1/g'`;
		uniquemap=`cat $in_dir/Log.final.out | grep "Uniquely mapped reads number |" | sed 's/^[ \t]\+Uniquely mapped reads number |[\t]\+\([0-9]\+\)$/unique_mapped\t\1/g'`;
		multimap=`cat $in_dir/Log.final.out | grep 'Number of reads mapped to multiple loci |' | sed 's/^[ \t]\+Number of reads mapped to multiple loci |[\t]\+\([0-9]\+\)$/multi_mapped\t\1/g'`;
		toomanymap=`cat $in_dir/Log.final.out | grep 'Number of reads mapped to too many loci |' | sed 's/^[ \t]\+Number of reads mapped to too many loci |[\t]\+\([0-9]\+\)$/toomany_mapped\t\1/g'`;
		echo -e "$header\n$totalpaired\n$uniquemap\n$multimap\n$toomanymap\n" > $out_dir/$sample"_alignment_summary.tab"
	fi
done
