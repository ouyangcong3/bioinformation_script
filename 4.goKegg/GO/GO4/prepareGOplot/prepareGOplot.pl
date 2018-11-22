use strict;
use warnings;

my $number=$ARGV[0];

my %hash=();
open(DIFF,"diff.txt") or die $!;
while(my $diff=<DIFF>){
	chomp($diff);
	my @arr=split(/\t/,$diff);
	$hash{$arr[0]}=$arr[1];
}
close(DIFF);

open(RF,"go.txt") or die $!;
open(DAVID,">david.txt") or die $!;
open(GENE,">gene.txt") or die $!;
print DAVID "Category\tID\tTerm\tGenes\tadj_pval\n";
print GENE "ID\tlogFC\n";
while(my $line=<RF>){
	next if($.==1);
	next if($.>($number+1));
	chomp($line);
	my @arr=split(/\t/,$line);
	my @category=split(/\_/,$arr[0]);
	my @term=split(/\~/,$arr[1]);
	print DAVID "$category[1]\t$term[0]\t$term[1]\t$arr[5]\t$arr[$#arr]\n";
	my @genes=split(/\,/,$arr[5]);
	foreach my $gene(@genes){
		$gene=~s/^\s+|\s+$//g;
		print GENE "$gene\t$hash{$gene}\n";
	}
}
close(GENE);
close(DAVID);
close(RF);

######ÊÓÆµ×ÊÁÏÍøÖ·£ºwww.xixibio.com
######ÌÔ±¦µêÆÌÍøÖ·£ºhttps://shop119322454.taobao.com
