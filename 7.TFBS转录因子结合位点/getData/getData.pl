######Video source: https://shop119322454.taobao.com
use strict;
use warnings;

######perl getData.pl gene.txt hg38_names.txt hg38_upstream.fa
my %geneHash=();

open(RF,"$ARGV[0]") or die $!;
while(my $line=<RF>)
{
	chomp($line);
	$line=~s/\s+//g;
	$geneHash{$line}=1;
}
close(RF);

my %hash=();
open(RF,"$ARGV[1]") or die $!;
while(my $line=<RF>)
{
	chomp($line);
	my @arr=split(/\t/,$line);
	if(exists $geneHash{$arr[1]})
	{
		$hash{$arr[0]}=$arr[1] . '|' . $arr[0];
	}
}
close(RF);

open(RF,"$ARGV[2]") or die $!;
open(WF,">gene.fa") or die $!;
my $flag=0;
while(my $line=<RF>)
{
	chomp($line);
	if($line=~/^>(.+?)\_(.+?)\_(.+?)\s+/)
	{
		$flag=0;
		my $trans=$3;
		if(exists $hash{$trans})
		{
			$flag=1;
			print WF ">$hash{$trans}\n";
			next;
		}
	}
	if($flag==1)
	{
		print WF $line . "\n";
	}
}
close(WF);
close(RF);
