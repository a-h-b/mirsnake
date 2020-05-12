#!/usr/bin/perl
#
use strict;

my $reads=$ARGV[0];
my $minlen=$ARGV[1];
my $maxlen=$ARGV[2];

my %reads=();

my ($id);
my ($seq);
my ($thisID);
my ($count);

open(CON,$reads) or die $!;
while(my $str=<CON>){
    chomp($str);
    next if length($str)==0;
    #print STDERR "str: $str\n";
    if($str=~/>(\S+)/){
	$id=$1;
	my @idcount=split("-",$id);
	$thisID=$idcount[0];
	$count=$idcount[1];
	#print STDERR "$thisID\n";
	$reads{$id}=1;
	die "$thisID already exists\n" if exists($reads{$thisID});
    }else{
	$seq.=$str;
	if (length($str)>=$minlen) {
	    if (length($str)<=$maxlen) {
	    print "$seq\t$count\n";
	    }
	}
	$count=0;	
	$seq="";
    }
} 
close(CON);
