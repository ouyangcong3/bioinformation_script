use strict;
use warnings;

my %hash=();
open(RF,"id.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	$hash{$arr[2]}="$arr[0]\t$arr[1]";
}
close(RF);

open(KEGG,"KEGG.txt") or die $!;
open(WF,">KEGG.xls") or die $!;
open(NODE,">node.txt") or die $!;
open(NET,">network.txt") or die $!;
print NODE "Id\tProperty\n";
print NET "Pathway\tGene\tInteraction\n";
while(my $line=<KEGG>){
	if($.==1){
		print WF $line;
		next;
	}
	chomp($line);
	my @arr=split(/\t/,$line);
	print NODE "$arr[0]\tkegg\n";
	my @idArr=split(/\//,$arr[$#arr-1]);
	my @symbols=();
	foreach my $id(@idArr){
		my @geneFC=split(/\t/,$hash{$id});
		print NET "$arr[0]\t$geneFC[0]\tKEGG\n";
		if($geneFC[1]>0){
		  print NODE "$geneFC[0]\tup\n";
		}
		else{
			print NODE "$geneFC[0]\tdown\n";
		}
		push(@symbols,$geneFC[0]);
	}
	$arr[$#arr-1]=join("/",@symbols);
	print WF join("\t",@arr) . "\n";
}
close(WF);
close(NODE);
close(NET);
close(KEGG);

###### ”∆µ◊ ¡œÕ¯÷∑£∫www.xixibio.com
######Ã‘±¶µÍ∆ÃÕ¯÷∑£∫https://shop119322454.taobao.com
