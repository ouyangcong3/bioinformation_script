#!/usr/bin/perl -w
use strict;
use warnings;

my $file=$ARGV[0];

#use Data::Dumper;
use JSON;
my $json = new JSON;
my $js;

open JFILE, "$file";
while(<JFILE>) {
	$js .= "$_";
}
my $obj = $json->decode($js);

open(WF,">time.txt") or  die $!;
print WF "id\tfutime\tfustat\n";
for my $i(@{$obj})
{
	my $vitalsStatus=$i->{'diagnoses'}->[0]->{'vital_status'};
	my $submitterId=$i->{'demographic'}->{'submitter_id'};
	my @subId=split(/\_/,$submitterId);
	if($vitalsStatus eq 'alive')
	{
		my $days_to_last_follow_up=$i->{'diagnoses'}->[0]->{'days_to_last_follow_up'};
		if(defined $days_to_last_follow_up)
		{
			print WF "$subId[0]\t$days_to_last_follow_up\t0\n";
		}
	}
	else
	{
		my $days_to_death=$i->{'diagnoses'}->[0]->{'days_to_death'};
		if(defined $days_to_death)
		{
			print WF "$subId[0]\t$days_to_death\t1\n";
		}
	}
}
close(WF);
#print Dumper $obj
