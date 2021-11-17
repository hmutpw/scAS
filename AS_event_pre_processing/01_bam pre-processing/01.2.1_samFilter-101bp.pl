#!/usr/bin/perl
open INPUT,"$ARGV[0]" or die $!;
open OUTPUT, "> $ARGV[1]" or die $!;
while(<INPUT>){
	chomp;
	my @line=split(/\t/,$_);
	if($line[0]=~/^@/){
		print OUTPUT "$_\n";
	}
	elsif(length($line[9]) eq 101){
		print OUTPUT "$_\n";
	}
}
close OUTPUT;
close INPUT;
