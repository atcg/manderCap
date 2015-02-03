#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;
use Cwd;

my $help = 0;
my $configFile;

GetOptions  ("config=s"      => \$configFile,
             "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$configFile or $help) {
    die "Must supply --config <file>.\n";
}
my $startingDir = getcwd();
open(my $configFH, "<", $configFile) or die "Couldn't open $configFile for reading: $!\n";

my %params;
while (my $line = <$configFH>) {
    if ($line =~ /^AdapterFolder\=(.*)/) {
        $params{"adapterFolder"} = $1;
    } elsif ($line =~ /^ReadsFolder\=(.*)/) {
        $params{"readsFolder"} = $1;
    } elsif ($line =~ /^QCoutputFolder\=(.*)/) {
        $params{"qcFolder"} = $1;
    } elsif ($line =~ /^ARCworkingFolder\=(.*)/) {
        $params{"arcWorkingFolder"} = $1;
    } elsif ($line =~ /^ARCoutputFolder\=(.*)/) {
        $params{"arcOutputFolder"} = $1;
    } elsif ($line =~ /^threads\=(.*)/) {
        $params{"threads"} = $1;
    } elsif ($line =~ /^Sample\=(.*)/) {
        push(@{$params{"samples"}}, $1);
    } elsif ($line =~ /^\s*$/) { # This is just a line that contains only whitespace (blank line)
        next;
    } else {
        die "Invalid config file line:\n$line\n";
    }
}
close($configFH);

unless (-d $params{"qcFolder"}) {
    mkdir $params{"qcFolder"} or die qq{Couldn't make output directory $params{"qcFolder"}: $!\n};
}

############################################################################
# Remove all CASAVA-filtered reads:
############################################################################
my $casavaFilteredDir = $params{"qcFolder"} . "/casavaFilt";
unless (-d $casavaFilteredDir) {
    mkdir $casavaFilteredDir;
}
my @fastqIlluminaFilterCommands;
foreach my $sample (@{$params{"samples"}}) {
    my $inFileNameR1 = $params{"readsFolder"} . "/" . $sample . ".R1.fastq";
    my $inFileNameR2 = $params{"readsFolder"} . "/" . $sample . ".R2.fastq";
    my $outFileNameR1 = $casavaFilteredDir . "/" . $sample . ".Ns.R1.fastq";
    my $outFileNameR2 = $casavaFilteredDir . "/" . $sample . ".Ns.R2.fastq";
    push(@fastqIlluminaFilterCommands, "fastq_illumina_filter -vN $inFileNameR1 > $outFileNameR1");
    push(@fastqIlluminaFilterCommands, "fastq_illumina_filter -vN $inFileNameR2 > $outFileNameR2");
}
my $forkManagerCasava = new Parallel::ForkManager($params{"threads"});
foreach my $casavaCommand (@fastqIlluminaFilterCommands) {
    print "Running the following command: \n$casavaCommand\n";
    sleep 5;
    $forkManagerCasava->start and next;
    system("$casavaCommand");
    print "Finished running the following:\n$casavaCommand\n\n";
    $forkManagerCasava->finish;
}
$forkManagerCasava->wait_all_children;
print "--------------------------------------------------\n";
print "Finished running all fastq_illumina_filter (CASAVA filter) commands\n";
print "--------------------------------------------------\n\n";


############################################################################
# Run Trimmomatic on all of the samples.
############################################################################
my $trimmomaticDir = $params{"qcFolder"} . "/trimmomatic";
unless (-d $trimmomaticDir) {
    mkdir $trimmomaticDir or die "Couldn't make directory $trimmomaticDir: $!\n";
}

my @trimmomaticCommands;
foreach my $sample (@{$params{"samples"}}) {
    my $R1 = $casavaFilteredDir . "/$sample" . ".Ns.R1.fastq";
    my $R2 = $casavaFilteredDir . "/$sample" . ".Ns.R2.fastq";
    my $adapterFile = $params{"adapterFolder"} . "/" . $sample . ".adapters";
    
    my $R1OutFilePaired = "$trimmomaticDir/$sample" . ".Ns.R1_paired_trimmed.fastq";
    my $R2OutFilePaired = "$trimmomaticDir/$sample" . ".Ns.R2_paired_trimmed.fastq";
    my $R1OutFileSingles = "$trimmomaticDir/$sample" . ".Ns.R1_singles_trimmed.fastq";
    my $R2OutFileSingles = "$trimmomaticDir/$sample" . ".Ns.R2_singles_trimmed.fastq";
    
    push(@trimmomaticCommands, "java -Xmx8G -jar ~/bin/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads $params{'threads'} -phred33 $R1 $R2 $R1OutFilePaired $R1OutFileSingles $R2OutFilePaired $R2OutFileSingles ILLUMINACLIP:$adapterFile:2:30:10 LEADING:5 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:40");
}
print "Running all Trimmomatic commands\n";
foreach my $trimmomaticCommand(@trimmomaticCommands) {
    print "Running the following command: \n$trimmomaticCommand\n";
    system($trimmomaticCommand);
}
print "Finished running all Trimmomatic commands\n";


############################################################################
# Run fastq-join on all the samples.
############################################################################
my $fastqJoinDir = $params{"qcFolder"} . "/fastq-join";
unless (-d $fastqJoinDir) {
    mkdir $fastqJoinDir or die "Couldn't make directory $fastqJoinDir: $!\n";
}

my @fastqJoinCommands;
chdir($fastqJoinDir);
foreach my $sample (@{$params{"samples"}}) {
    my $fastqJoinOutputName = $sample . ".Ns.%.fastq";
    push(@fastqJoinCommands,"fastq-join -v ' ' ../trimmomatic/$sample.Ns.R1_paired_trimmed.fastq ../trimmomatic/$sample.Ns.R2_paired_trimmed.fastq -o $fastqJoinOutputName");
}
my $forkManagerFastqJoin = new Parallel::ForkManager($params{"threads"});
foreach my $fastqJoinCommand (@fastqJoinCommands) {
    print "Running the following command: \n$fastqJoinCommand\n";
    sleep 5;    
    $forkManagerFastqJoin->start and next;
    system($fastqJoinCommand);
    print "Finished running the following command: \n$fastqJoinCommand\n";
    $forkManagerFastqJoin->finish;
}
$forkManagerFastqJoin->wait_all_children;

chdir($startingDir); # Go back to the starting directory


############################################################################
# Combine the singleton and joined read files
############################################################################
foreach my $sample (@{$params{"samples"}}) { # Don't need threading here because this is fast
    my $combinedFile = $fastqJoinDir . "/$sample.Ns.combinedJoinedAndSingles_trimmed.fastq";
    my $singlesFile1 = "$trimmomaticDir/$sample.Ns.R1_singles_trimmed.fastq"; 
    my $singlesFile2 = "$trimmomaticDir/$sample.Ns.R2_singles_trimmed.fastq";
    my $joinedFile   = "$fastqJoinDir/$sample.Ns.join.fastq";
    print "Running the following command: \ncat $joinedFile $singlesFile1 $singlesFile2 > $combinedFile\n\n";
    system("cat $joinedFile $singlesFile1 $singlesFile2 > $combinedFile");
    unlink($singlesFile1, $singlesFile2, $joinedFile);
}





































