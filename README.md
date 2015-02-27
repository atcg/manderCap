*Summary:*
========
This repository holds the code used to analyze data from an exon capture experiment.

The general strategy is as follows:

  1. Quality control data (remove reads not passing CASAVA filter, remove adapter contamination, window-trim for sequence quality, merge overlapping reads)
  2. Combine data for all of the libraries for a single species into a single set of fastq files
  3. Perform targeted assemblies of our exon data using the [ARC assembly pipeline from iBest](http://ibest.github.io/ARC/)
  4. Map reads from individual libraries back to the assemblies made with ARC to assess how different laboratory conditions affected capture efficiency
  
The pipeline is run from a driver script called `runAnalyses.pl`. You must provide it with a configuration file as follows: `perl runAnalyses.pl -c configurationFile.txt`. The configuration file contains information regarding the names of the sample libraries and where they are located. Its form should be as follows (it includes one header row and should contain no whitespace):

```
AdapterFolder=./adapters
ReadsFolder=./concatenatedReads
QCoutputFolder=./qc
ARCworkingFolder=/home/evan/ramdisk
threads=30
targets=targets.fasta
cycles=6
timeout=180

Sample=SampleID1   
Sample=SampleID2
Sample=SampleID3
...
Sample=SampleIDlast

```
  
  
*Initial housekeeping:*
=====================

The manderCap repository should be in the same directory as the directories containing the read data (in this case those directories are called Project_Shaffer and Undetermined_indices). We'll call this top-level directory that contains the manderCap repo and the Project_Shaffer and Undetermined_indeces folders "topDirectory" from here on out.


Move all the reads from their Sample folders into a single directory:

`cd topDirectory`

`mkdir allReads`

`mv Project_Shaffer/Sample_0*/*.fastq.gz allReads/`

We don't need the empty sample folders, Sample Sheets, or Basecall_Stats. We also won't use the reads that weren't demultiplexed. So nuke those:

`rm -r Project_Shaffer`

`rm -r Undetermined_indices`

Now we want to make a single R1 and single R2 file for each sample. We'll put those into a new directory:

`mkdir concatenatedReads`

`cd allReads`

`bash ../manderCap/concatenateSampleReads.sh > ../logs/concatenateSampleReads.log 2>&1`

We also want the read files to be unzipped:

`cd ../concatenatedReads; gunzip *.gz > ../logs/gunzipSampleReads.log 2>&1`


Finally, we want to produce the sample-specific adapter FASTAs that we'll use to trim adapter contamination:

`cd topDirectory`

`mkdir adapters`

`perl manderCap/generate_adapter_fastas.pl --in manderCap/HSEM020_adaptersKey.txt --out adapters --adapters itru > logs/generate_adapter_fastas.log 2>&1`


Run the QC and assemblies:
-------------------------
The following command runs all the read QC and ARC assemblies (see that script for details):

`perl manderCap/runAnalyses.pl -c manderCap/manderCap.config > logs/runAnalyses.stdout 2>logs/runAnalyses.stderr`


Assembly preparation:
---------------------
After ARC has run, we'll process the ARC assembled contigs to arrive at our final
reference that we'll use to map all of our individual libraries to. We first want
to designate a single assembled contig that is representative of each target. To do
this, we'll find the Reciprocal Best Blast Hits (RBBHs) for each target. This is
performed by manderCap/findRBBHs.pl. Briefly, it blasts the contigs in the assembly
against the sequences in the target fasta file. Then it blasts the targets in the
fasta file against the assembly. Then it goes through these blast reports and finds
all instances where the targets had at least one hit to the assembly. For all of the
best hits, it checks to make sure that that target was also found as the best hit
when we blasted the assembly against the targets (hence the reciprocal part).


`cd ARC/finished_allCTS/`

`cp contigs.fasta contigs.CTSonly6iter.fasta`

`../../manderCap/findRBBHs.pl --assembly contigs.CTSonly6iter.fasta --targets ../../targets.fasta --out RBBHs.CTSonly6iter.fasta > ../../logs/findRBBHs.log 2>&1`


That process finds a total of 8386 reciprocal best blast hits (over 96%). If we were
to map reads to this assembly and visualize coverage, we'd see a fair amount of targets
that have spikes in read coverage at the ends of the targets. This can be the result
of several things. It's possible that these are true repetitive regions that exist
at the edges of our target (at the periphery of the exons, for instance). They may
also be due to chimerism in our assemblies, whereby repetitive, non-contiguous genomic
regions are grafted onto the edges targets due to challenges in <em>de novo</em> assembly.
Here is one example of what what I'm talking about (the black bar represents the portion
of the target that blasted to that assembled contig):
![chimeraSpike](images/chimeraSpike.png)

To minimize this effect, we'll try to mask some of these repetitive regions prior to
mapping. Since a given repetitive region is typically present in more than one
representative target contig, we'll do this by masking regions that blast to multiple
targets via self-blast:

`makeblastdb -in RBBHs.CTSonly6iter.fasta -dbtype nucl > ../../logs/makeblastdb.log 2>&1`

`blastn -db RBBHs.CTSonly6iter.fasta -query RBBHs.CTSonly6iter.fasta -out RBBHs.CTSonly6iter.selfBlast -evalue 1e-20 -outfmt 6`

`perl ../../manderCap/maskChimeras.pl --in RBBHs.CTSonly6iter.selfBlast --sequences RBBHs.CTSonly6iter.fasta --out RBBHs.CTSonly6iter.chimeraMasked.fasta > ../../logs/maskChimeras.log 2>&1`

```
get_fasta_lengths.py --input RBBHs.CTSonly6iter.chimeraMasked.fasta
  Reads:		8,386
  Bp:		11,813,530
  Avg. len:	1,408.72048653
  STDERR len:	3.4372379005
  Min. len:	210
  Max. len:	6,284
  Median len:	1,398.0
  Contigs > 1kb:	7,816
```
We'll use this as our assembly for downstream analyses:

`cd ../..;mkdir assembly;cp ARC/finished_allCTS/RBBHs.CTSonly6iter.chimeraMasked.fasta assembly/;`

`bwa index RBBHs.CTSonly6iter.chimeraMasked.fasta > ../logs/indexAssembly.log 2>&1`

Individual read mapping:
------------------------
Now we'll map reads to the assembly we prepared in the previous step. The general strategy is
to map single reads, then paired end reads. Then we merge the bam files (samtools), clean them
(picard), sort them (picard), and mark duplicates (picard). The MarkDuped files are then
analyzed using flagstat (samtools), and depth-calculated (samtools depth). Those depth files
will be used to evaluate the individual loci. All of these steps are implemented in manderCap/mapReads.pl.

```
cd topDir; mkdir mapping; mkdir mapping/depths
perl manderCap/mapReads.pl --config manderCap/manderCap.config --reference /mnt/Data1/HSEM020/assembly/RBBHs.CTSonly6iter.chimeraMasked.fasta > logs/readMapping.log 2>&1
```

We want to be clear that these bam files and depth files are the result of mapping
to the CTS-only 6 ARC iteration assembly. So we'll rename the mapping and depths folders:
```
mv mapping mappingCTS6iter
cd mappingCTS6iter
mv depths depthsCTS6iter
```

Generate the locus-specific read depths:<br>
```
cd depthsCTS6iter
perl ../../manderCap/createLocusDepthFiles.pl --out locusDepthFilesCTSonly6iter > ../../logs/createLocusDepthFiles.CTSonly6iter.log 2>&1
```


Generate the locus-specific plots of sequencing depth across the loci:<br>
```
cd locusDepthFilesCTSonly6iter;
ls -al | grep "\.depth" | awk '{print $(NF)}' > dpthFilesCTSonly6iter.txt
```



Evaluation of different libraries:
----------------------------------
Here we process the depth files generated by samtools depth on the MarkDuped files.






Target-level metrics:
---------------------

`cd ARC/finished_allCTS/`<br>
`makeblastdb -in RBBHs.CTSonly6iter.chimeraMasked.fasta -dbtype nucl`<br>
`blastn -db RBBHs.CTSonly6iter.chimeraMasked.fasta -query ../../targets.fasta -outfmt 6 -evalue 1e-20 -out targets_bl2_chimeraMaskRBBHs.1e-20.blast`

`perl ../../manderCap/countSuccesfulTargets.pl --blast targets_bl2_chimeraMaskRBBHs.1e-20.blast --assembly RBBHs.CTSonly6iter.chimeraMasked.fasta --depthFileFolder ../../mapping/depths/locusDepthFiles --out` 





Alternative assembly creation and evaluation:
---------------------------------------------
We also wanted to generate assemblies using 6 ARC iterations on all of the read
data from the 3 individuals. Here is the ARC config file for that:
```
## Name=value pairs:
## reference: contains reference sequences in fasta format
## numcycles: maximum number of times to try remapping
## mapper: the mapper to use (blat/bowtie2)
## assembler: the assembler to use (newbler/spades)
## nprocs: number of cores to use
## format: fasta or fasta, all must be the same
## verbose: control mapping/assembly log generation (True/False)
## urt: For Newbler, enable use read tips mode (True/False)
## map_against_reads: On iteration 1, skip assembly, map against mapped reads (True/False)
## assemblytimeout: kill assemblies and discard targets if they take longer than N minutes
##
## Columns:
## Sample_ID:Sample_ID
## FileName: path for fasta/fasta file
## FileType: PE1, PE2, or SE
## FileFormat: fasta or fasta
# reference=/mnt/Data1/HSEM020/targets.fasta
# numcycles=6
# mapper=bowtie2
# assembler=spades
# nprocs=20
# format=fastq
# verbose=True
# urt=True
# map_against_reads=False
# assemblytimeout=180
# bowtie2_k=5
# rip=True
# cdna=False
# subsample=1
# maskrepeats=True
# workingdirectory=/home/evan/manderReads/
Sample_ID	FileName	FileType
all3	All3.Ns.un1.fastq	PE1
all3	All3.Ns.un2.fastq	PE2
all3	All3.Ns.combinedJoinedAndSingles_trimmed.fastq	SE
```










