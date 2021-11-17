#!/usr/bin/perl
open FASTA_IN,"$ARGV[0]" or die $!;
open SAM_IN,"$ARGV[1]" or die $!;
open SAM_OUT, "> $ARGV[2]" or die $!;

#get transcript ids from fasta header.
my %Transcrpt_id;
while(<FASTA_IN>){
        chomp;
        if($_=~/^>/){
        my @line=split(/\|/,$_);
        my $trans_id=$line[0];
        $trans_id=~/^>(.+)/;
        my $transcript_id=$1;
        $Transcrpt_id{$transcript_id}++;
        }
}

#filter sam files with transcript id only in fasta files.
while(<SAM_IN>){
        chomp;
        my @line=split(/\t/,$_);
        if($line[0]=~/^@/){
		$line[1]=~/^SN:(.+)/;
		if( exists $Transcrpt_id{$1}){
			print SAM_OUT "$_\n";
		}
        }
        elsif(exists $Transcrpt_id{$line[2]} ){
			print SAM_OUT "$_\n";
        }
}

close OUTPUT;
close SAM_IN;
close SAM_OUT;
