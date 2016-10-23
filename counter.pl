#!/usr/bin/perl
use warnings;
use strict;

print "This program will count the bases per sequence in a multifasta file.\n";
print "MultiFasta: ";
my $INPUT = <STDIN> ;
open(VAR, "<", $INPUT) or die("Incorrect path: $INPUT\n");
my %Sequences;
my $seqid;
while(<VAR>){
    my $line = $_; chomp($line);
    if ($line =~ m/^>(\S+)/){   
        $seqid = $1;             
        $Sequences{$seqid} = "";  
        }
    else {
        $Sequences{$seqid} = $Sequences{$seqid} . $line;
        }
   }
foreach my $sequence_entry (keys %Sequences){   
    my $currentSequence = $Sequences{$sequence_entry};
    my $lengthSequence = length($currentSequence);
    print $sequence_entry . " : " . $lengthSequence . " bases" . "\n";
}
