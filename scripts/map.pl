#!/usr/bin/env perl

use strict;
use Getopt::Std;
use File::Basename;
use warnings;
use threads;

die(qq/

Usage: seqCapture map  [options]

Options:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-f    FILE    Reference fasta sequence file 
              (for example: ABC_123.fasta)
-r    DIR     Path to the directory of cleaned reads
-o    DIR     Path to the results directory
-p    INT     number of threads [1]
-e    INT     if single end reads e = 1
              if paired end reads e = 2 [2]  
-c    INT     only keep concordant mapping for PE reads?
              1 = yes
              0 = no [1]
-h    INT     remove reads that are mapped non-uniquely, even 
              only the best is reported?
              1 = yes
              0 = no [1]
-k    INT     keep priginal sam files? 
              1 = yes
              0 = no [1]
-m    INT     Method for alignment
              1 = novoalign
              2 = bwa [1]

additional options if use novoalign (-m 1): 
-i    INT     if PE, Avg. insert size for the libraries 
-v    INT     if PE, STD insert size for the libraries (usually 0.1*i)
-l    INT     read length [100] 
-t    INT     Maximum alignment score acceptable for the best 
              alignment. Use 90 for population samples (roughly 
              3 mismatches are allowed per pair). If use 150, 
              roughly 5 mismatches are allowed per read in PE. 
              For very divergent genomes, use default value 
              by not defining -t [null]


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

\n\n/) unless (@ARGV);


my %opts = (f => undef,r => undef, v=>undef, o=>undef, c=>1, i=>undef, t=>undef, l=>100, p => 1, e=>2, m=>1, h=>1, k=>1);
getopts('f:r:o:i:t:v:l:P:G:p:e:c:m:h:k:', \%opts);

my $data_dir = redir ($opts{r});
my $res_dir = redir ($opts{o});
mkdir $res_dir unless -e $res_dir;

my $keep = $opts{k};
my $unique = $opts{h};
my $method = $opts{m};

my $nb_process = $opts{p};
my $in_target = $1 if  $opts{f} =~ m /(\S+)\.fa/;
my $ref = $opts{f};
my $pair = $opts{e};
my $con = $opts{c};

## only for novoalign:
my $score;
my $insert;
my $std;
my $length;

if ($method == 1) {
  $score = $opts{t} if $opts{t};
  $score =0 if ! $opts{t};
  $insert = $opts{i};
  $std = $opts{v};
  $length = 2 * $opts{l};
}

my @reads1 = <$data_dir*_1_final.fq>;
my $nb_compute = scalar @reads1;

my $indexed_assemblies_in_target = $in_target  . ".nix";

if ($method == 1) {
  system ("novoindex $indexed_assemblies_in_target $ref")  ;
  
}

if ($method == 2) {
  system ("bwa index $ref");
}

multithread ();  

foreach (@reads1) { 
 alignmentruns ($_);
}

sub alignmentruns {
  my ($read1) = @_;
  my $lib = $1 if basename($read1) =~ m/(\S+)_1_final/;
  my $sorted_in_target_bams = $res_dir . $lib . "_sorted";
  my $sorted_in_target_bams2 = $sorted_in_target_bams. '.bam';
  
  if ($pair == 2){  
    my $read2 = $read1;
    $read2 =~ s/_1_final/_2_final/;
    my $readSolo = $read1; $readSolo =~ s/_1_final/_u_final/;   
    my $outPairedSam1 = $res_dir . $lib . ".outPairedSam1";
    my $outSoloSam1 = $res_dir . $lib . ".outSoloSam1";
    my $paired_in_target =  $res_dir . $lib . ".paired_in_target.sam";
    my $solo_in_target = $res_dir . $lib . ".solo_in_target.sam";
    my $paired_in_target_bam = $res_dir . $lib ."paired_in_target.bam";
    my $solo_in_target_bam = $res_dir . $lib ."solo_in_target.bam";
    my $rawbam = $res_dir . $lib . ".raw.bam";
    
    if ($method == 1) {
      my $soloscore = $score + 30;
      if ($score > 0) {
	my $alignment_paired_to_in_target_assemblies= system("novoalign -R 30 -t $score -n $length  -d $indexed_assemblies_in_target -f $read1  $read2  -i PE $insert, $std  -a AGATCGGAAGAGCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG  AGATCGGAAGAGCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $outPairedSam1");
	
	my $alignment_solo_to_in_target_assemblies = system("novoalign -R 30  -t $soloscore    -d $indexed_assemblies_in_target -f $readSolo  -a AGATCGGAAGAGCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG  AGATCGGAAGAGCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -n $length -o SAM > $outSoloSam1");
	
      }
      
      if ($score == 0) {     
	my $alignment_paired_to_in_target_assemblies= system("novoalign -R 30 -n $length   -d $indexed_assemblies_in_target -f $read1  $read2  -i PE $insert, $std  -a AGATCGGAAGAGCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG  AGATCGGAAGAGCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA  -F STDFQ -o SAM > $outPairedSam1 ");
	
	my $alignment_solo_to_in_target_assemblies = system("novoalign -R 30  -d $indexed_assemblies_in_target -f $readSolo  -a AGATCGGAAGAGCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG  AGATCGGAAGAGCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -n $length -o SAM >  $outSoloSam1 ");
      }
    }
    
    if ($method == 2) {
      system ("bwa mem -t $nb_process $ref $read1  $read2  > $outPairedSam1");
      system ("bwa mem -t $nb_process $ref  $readSolo > $outSoloSam1");
    }

    
    if ($con == 1) {
      my $outPairedSam2 =  $res_dir . $lib . ".outPairedSam2";
      open (IN, "<", $outPairedSam1);
      open (OUT, ">", $outPairedSam2 );
      
      while (<IN>) {
	chomp (my $line = $_);
	if ($line =~ /^@/) {
	  print OUT $line, "\n";
	}
	else {
	  chomp (my @a = split /\s+/,$line);
	  if ( $a[1] == 99 || $a[1] == 147 || $a[1] == 83 || $a[1] == 163 ) {
	    print OUT $line, "\n";
	  }
	  
	}
	
      }
      close IN;
      close OUT;      
      system ("mv $outPairedSam2 $outPairedSam1");      
    }
    
    if ($unique == 1) {
      system("grep -v ZS:Z:NM $outPairedSam1 | grep -v NH:i: | grep -v XA:Z: |  grep -v SA:Z: > $paired_in_target");
      system("grep -v ZS:Z:NM $outSoloSam1 | grep -v NH:i: | grep -v XA:Z: |  grep -v SA:Z: > $solo_in_target");
    }
    if ($unique == 0) {
      system("cp $outPairedSam1 $paired_in_target");
      system("cp $outSoloSam1 $solo_in_target");
    }
    if ($keep == 0) {
      unlink ($outPairedSam1, $outSoloSam1) if ($outPairedSam1 && $outSoloSam1);
    }
    
    system("samtools view -bS $paired_in_target >  $paired_in_target_bam");
    system("samtools view -bS $solo_in_target > $solo_in_target_bam");

    unlink ($paired_in_target, $solo_in_target);
    system("samtools merge  --reference $ref $rawbam $paired_in_target_bam $solo_in_target_bam");
    unlink ($paired_in_target_bam, $solo_in_target_bam);

    system("samtools sort $rawbam -o $sorted_in_target_bams2 ");
    unlink ($rawbam);   
    system("samtools index $sorted_in_target_bams2");      
  }


  if ($pair == 1){  
    my $lib = $1 if basename($read1) =~ m/(\S+)_1_final/;  
    my $outSoloSam1 = $res_dir . $lib . ".outSoloSam1";   
    my $solo_in_target = $res_dir . $lib . ".solo_in_target.sam";
    my $solo_in_target_bam = $res_dir . $lib ."solo_in_target.bam";


    if ($method == 1) {
      if ($score > 0) {     
	my $alignment_solo_to_in_target_assemblies = system("novoalign -R 30  -t $score   -d $indexed_assemblies_in_target -f $read1  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -n $length -o SAM > $outSoloSam1");     
      }
      
      if ($score == 0) {     
	my $alignment_solo_to_in_target_assemblies = system("novoalign -R 30  -d $indexed_assemblies_in_target -f $read1  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -n $length -o SAM > $outSoloSam1 ");
      }
    }
    if ($method == 2) {
      system ("bwa mem -t $nb_process $ref $read1 > $outSoloSam1");
    }
 

    if ($unique == 0) {
      system("cp $outSoloSam1 $solo_in_target");
    }
    if ($unique == 1) {
      system("grep -v ZS:Z:NM $outSoloSam1 | grep -v NH:i: | grep -v XA:Z: |  grep -v SA:Z: > $solo_in_target");
    }
    
    if ($keep == 0) {
      unlink ($outSoloSam1);
    }

    
    system("samtools view -bS $solo_in_target > $solo_in_target_bam");
    
    unlink ($solo_in_target);
    
    system("samtools sort $solo_in_target_bam -o $sorted_in_target_bams2");
     
    unlink ($solo_in_target_bam);

    system("samtools index $sorted_in_target_bams2");  
    
  }
  
  
  my $indref = $in_target . "_" . $lib. ".fasta";
  system ("cp $ref $indref");
  my $index_ref =  $in_target . "_" . $lib.  '.dict';
  
  system ("picard -Xms512m -Xmx4g CreateSequenceDictionary R=$indref O=$index_ref");  
  system ("samtools faidx $indref");
  
  my $addgroupbam = $sorted_in_target_bams2. ".rg.bam"; 
  system("picard -Xms512m -Xmx4g AddOrReplaceReadGroups INPUT=$sorted_in_target_bams2 OUTPUT=$addgroupbam RGID=$lib RGLB=RAD RGPL=illumina RGPU=lane1 RGSM=$lib");
  
  my $addgroupbamsort = $addgroupbam . "_sorted";
  my $addgroupbamsorted = $addgroupbamsort . ".bam";
  
  system("samtools sort $addgroupbam -o $addgroupbamsorted ");
  system ("samtools index $addgroupbamsorted");
  
  my $intervals = $addgroupbamsorted .".intervals";
  system ("gatk -Xms512m -Xmx4g -T RealignerTargetCreator -R $indref -I $addgroupbamsorted -o $intervals --num_threads $nb_process ");
  my $final = $addgroupbamsorted. "1";
  my $final_sort = $addgroupbamsorted . "2";
  
  system ("gatk -Xms512m -Xmx4g -I  $addgroupbamsorted -R $indref -T IndelRealigner -targetIntervals $intervals -o $final ");
  
  my $final_sorted = $final_sort . ".bam";
  system ("samtools sort $final -o $final_sorted");
  
  system ("mv $final_sorted $sorted_in_target_bams2");
  system ("samtools index $sorted_in_target_bams2");
  system ("rm $addgroupbamsorted* $addgroupbam $indref* ");
  
}

sub redir {
  my ($dir) = @_;
  my $out;
  if ($dir =~ m/\/$/ ){
    $out = $dir; 
  }
  else {
    $out = $dir . "/";
  }
  return ($out);
}

sub multithread {  
  my @running = ();
  my @Threads;
  while (scalar @Threads < $nb_compute) {
    
    @running = threads->list(threads::running); 
    if (scalar @running < $nb_process) {
      my $line = shift @reads1;
      my $libs = $1 if basename ($line) =~ m/(\S+)_[1|2|u]_final\.fq/;
      my $thread = threads -> new(\&alignmentruns, $line);
      push (@Threads, $thread);
      my $tid = $thread->tid;
      print "  - starting alignment $libs!\n";
    }
    @running = threads->list(threads::running);
    foreach my $thr (@Threads) {
      if ($thr->is_running()) {
	my $tid = $thr->tid;
      }
      
      elsif ($thr->is_joinable()) {
	my $tid = $thr->tid;
	$thr->join;
      }
      
    }
    @running = threads->list(threads::running);
  }
  while (scalar @running != 0) {
    @running = threads->list(threads::running);
    foreach my $thr (@Threads) {
      $thr->join if ($thr->is_joinable());
    }
  }
  
  
  
}

