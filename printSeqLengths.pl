#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $help = 0;
my $inFile;
my $outFile;

GetOptions  ("in=s"      => \$inFile,
             "out=s"      => \$outFile,
             "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$inFile or !$outFile or $help) {
    die "Must supply --in and --out.\n";
}

my $seqIn = Bio::SeqIO->new(-file => $inFile,
                            -format => 'fasta');

open(my $outFH, ">", $outFile) or die "Couldn't open $outFile for writing: $!\n";

while (my $seq = $seqIn->next_seq()) {
    print $outFH $seq->display_id . "\t" . $seq->length . "\n";
}



