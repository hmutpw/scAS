#!/usr/bin/perl

#commond: perl 03.1.1.get_ERCC_transcript_fasta.pl ERCC92.gtf ERCC92.fa ERCC92.transcript.fa
open GTF_IN,"$ARGV[0]" or die $!;
open OUTPUT, "> $ARGV[1]" or die $!;

# get transcript gene id relation.
print OUTPUT "transcript_id".","."gene_id\n";
%Trans2gene;
while(<GTF_IN>){
        chomp;
        my @line=split(/\t/,$_);
		if($line[2]=~/exon/){
		my @infor=split(/;/,$line[8]);
		$infor[0]=~/^gene_id \"(.+)\"/;
		my $gene_id=$1;
		$infor[1]=~/transcript_id \"(.+)\"/;
		my $transcript_id=$1;
		my $newline="$transcript_id".","."$gene_id";
		$Trans2gene{$newline}++;
		}
}

foreach my $key (sort keys %Trans2gene){
print OUTPUT "$key\n";
}

close OUTPUT;
close GTF_IN;