use strict;
use warnings;

my %hash=();

open(RF,"input.txt") or die $!;
open(WF,">PicTar.txt") or die $!;
while(my $line=<RF>)
{
	next if($line=~/^\n/);	
	chomp($line);
	my @arr=split(/\t/,$line);
	if($arr[$#arr]=~/.+\((.+?)\)\,/)
	{
		unless(exists $hash{$1})
		{
		  print WF $1 . "\n";
		  $hash{$1}=1;
		}
	}
	else
	{
		print $line . "\n";
	}
}
close(WF);
close(RF);
