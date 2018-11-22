use strict;
use warnings;

my %hash=();

open(RF,"hg38.genenames") or die $!;
while(my $line=<RF>)
{
	next if($.==1);
	chomp($line);
	my @arr=split(/\t/,$line);
	$hash{$arr[0]}=$arr[1];
}
close(RF);

open(RF,"hg18_0.gtf") or die $!;
open(WF,">hg38_0.gtf.cp") or die $!;
while(my $line=<RF>)
{
	$line=~s/gene_id \"(.+?)\"\;/gene_id \"$hash{$1}\"\;/g;
	print WF $line;
}
close(WF);
close(RF);

