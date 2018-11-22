use strict;
use warnings;

my %hash=();
my %sampleHash=();
my @sampleIndex=();
my @sampleArr=();

open(DOSE,"$ARGV[0]") or die $!;
while(my $dose=<DOSE>)
{
	next if($.==1);
	chomp($dose);
	my @doseArr=split(/\t/,$dose);
	${$hash{$doseArr[3]}}{$doseArr[2]}=$doseArr[5];
	$sampleHash{$doseArr[2]}=0;
}
close(DOSE);

open(EXP,"$ARGV[1]") or die $!;
open(WF,">expression.txt") or die $!;
while(my $exp=<EXP>)
{
	chomp($exp);
	my @expArr=split(/\t/,$exp);
	next if($expArr[0] eq '');
	if($.==1)
	{
		print WF "ID";
		for(my $i=2;$i<=$#expArr;$i++)
		{
			if($expArr[$i]=~/DATA\.(\d+)/)
			{
				my $expSample=$1;
				print $expSample . "\n";
				if(exists $sampleHash{$expSample})
				{
					print WF "\t$expSample";
					push(@sampleIndex,$i);
					push(@sampleArr,$expSample);
					delete($sampleHash{$expSample});
				}
			}
		}
		print WF "\n";
	}
	else
	{
		print WF $expArr[0];
		foreach my $index(@sampleIndex)
		{
			print WF "\t$expArr[$index]";
		}
		print WF "\n";
	}
}
close(WF);
close(EXP);

open(WFDRUG,">drug.txt") or die $!;
print WFDRUG "ID";
foreach my $sample(@sampleArr)
{
	print WFDRUG "\t$sample";
}
print WFDRUG "\n";
foreach my $key(keys %hash)
{
	print WFDRUG $key;
	foreach my $sample(@sampleArr)
	{
		if(exists ${$hash{$key}}{$sample})
		{
			print WFDRUG "\t" . ${$hash{$key}}{$sample};
		}
		else
		{
			print WFDRUG "\t-100";
		}
	}
	print WFDRUG "\n";
}
close(WFDRUG);
