use strict;
use warnings;

my $blastOutFile=$ARGV[0];
my $allGeneFile=$ARGV[1];
my $diffGeneFile=$ARGV[2];
my $goFastaFile=$ARGV[3];
my $gene2goFile=$ARGV[4];

my %blastHash=();
open(RF,"$blastOutFile") or die $!;
while(my $line=<RF>)
{
    chomp($line);
    my @arr=split(/\t/,$line);
    push(@{$blastHash{$arr[1]}},$arr[0]);
}
close(RF);

my %goHash=();
open(RF,"$goFastaFile") or die $!;
while(my $line=<RF>)
{
    if($line=~/^>(.+?)\s+/)
    {
	chomp($line);
	my $seqId=$1;
	if(exists $blastHash{$seqId})
	{
	    my @arr=split(/\[|\]/,$line);
	    foreach my $go(@arr)
	    {
		if($go=~/^(GO.+?)\s+\"(.+?)\"/)
		{
		    my $goId=$1;
		    foreach (@{$blastHash{$seqId}})
		    {
			${$goHash{$goId}}{$_}=1;
		    }
		}
	    }
	}
    }
}
close(RF);

my %allHash=();
open(RF,"$allGeneFile") or die $!;
while(my $line=<RF>)
{
    if($line=~/^>(.+?)(\s+|\n)/)
    {
	$allHash{$1}=1;
    }
}
close(RF);
my $allLength=keys %allHash;

my %diffHash=();
open(RF,"$diffGeneFile") or die $!;
while(my $line=<RF>)
{
    if($line=~/^>(.+?)(\s+|\n)/)
    {
	$diffHash{$1}=1;
    }
}
close(RF);
my $diffLength=keys %diffHash;

my %gene2goHash=();
open(RF,"$gene2goFile") or die $!;
while(my $line=<RF>)
{
    next if($.==1);
    chomp($line);
    my @arr=split(/\t/,$line);
    if($arr[7] eq 'Process')
    {
	$arr[7]='Biological_Process';
    }
    elsif($arr[7] eq 'Function')
    {
	$arr[7]='Molecular_Function';
    }
    elsif($arr[7] eq 'Component')
    {
	$arr[7]='Cellular_Component';
    }
    $gene2goHash{$arr[2]}=$arr[5] . "\t" . $arr[7];
}
close(RF);

open(WF,">goAnnResult.txt") or die $!;
print WF "GOId\tTerm\tCategory\tGene\tX\tN\tx\tn\n";
foreach my $goKey(keys %goHash)
{
    my %goValue=%{$goHash{$goKey}};
    my $goValueLength = keys %goValue;
    my $allGoCount=0;
    my $diffGoCount=0;
    my @diffGeneArr=();
    foreach my $goValueKey(keys %goValue)
    {
	if(exists $allHash{$goValueKey})
	{
	    $allGoCount++;
	}
	if(exists $diffHash{$goValueKey})
	{
	    $diffGoCount++;
	    push(@diffGeneArr,$goValueKey);
	}
    }
    if(($diffGoCount>0) && (exists $gene2goHash{$goKey}))
    {
	print WF "$goKey\t$gene2goHash{$goKey}\t" . join("|",@diffGeneArr) . "\t$allGoCount\t$allLength\t$diffGoCount\t$diffLength\n";
    }
}
close(WF);
