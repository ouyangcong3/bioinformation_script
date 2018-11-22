use strict;
use warnings;

#video source: https://shop119322454.taobao.com

my $file=$ARGV[0];
my $mirna=$ARGV[1];

open(RF,"$file") or die $!;
open(WF,">miRanda.xls") or die $!;
while(my $line=<RF>)
{
	my @arr=split(/\t/,$line);
	if(($.==1)||($arr[1] eq $mirna))
	{
		print WF $line;
	}
}
close(WF);
close(RF);
