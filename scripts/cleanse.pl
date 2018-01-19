#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use List::Util qw(sum);
use Cwd 'abs_path';

die (qq/

Usage: seqCapture cleanse  [options]


Options:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-f     DIR      Path to the folder of raw fastq files
-o     DIR      Path to the results folder
-c     FILE     Contaminant file
-m     FLOAT    percent similarity to contaminant file to
                consider a read to be a contaminant
                [0.90]
-k     INT      Number of threads needed for bowtie2 alignment 
                (remove contamination), skewer(trimming adapters),
                trimmomatic (trimming adapters and low qual), read 
                error correction (for de novo assemblies) [4]
-q     INT      Quality trimming threhold [20]
-d     CHAR     The particular library that you like to process,
                if "all" is specified, all libraries in the folder
                (-f) will be processed sequenctially [all]
-h     INT      Trimmomatic trimming length cutoff [36]
-z              If z is supplied, use fastQC to evaluate
                cleaned sequence reads [null]

####  Error correction is only recommended for de novo assemblies  ####
####  DO NOT run error correction for alignment based analyses ####

-N     INT      error correction
                0 = no error correction -- default
                1 = using bfc (when number of reads is < 30 million)
                2 = using Rcorrector (when number of reads is >= 30 million)
                [0]
-G     INT      if N=1 approx genome size (k\/m\/g allowed)
                for transcriptomic data set this to 50m [50m] 


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


\n/) if !@ARGV;

my %opts = (f=>undef, o=>undef, c=>undef,  k=> 4, q=> 20,  d=>"all", h=>36, N=>0, G =>"50m", m=>0.9);
getopts('f:o:c:d:h:k:q:N:G:zv', \%opts);


my $rawdir = redir ($opts{f});
my $outdir = redir ($opts{o});

mkdir ($outdir) unless -d ($outdir);

my $contam = $opts{c};
die(qq/\nHmm...I could not find where the contaminant file is.. Check "-c"? \n\n/) unless -e $contam;

###checking all dependencies required
my $spath = $1 . 'dependencies/morePackages/'  if  (dirname(abs_path($0)) =~ /(\S+)scripts/) or die "scripts path is incorrect!\n";
my $trimmomatic = 'trimmomatic';
my $bowtie = 'bowtie2';
my $cutadapt = 'cutadapt';
####


#universal adapter P5 if no barcode is used in it
my $uniad = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT';

my $P1 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA';
my $P2 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC';
my $P1rc = 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT';
my $P2rc = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT';

my $cpu = $opts{k};
my $quality = $opts{q};
my $minLength = $opts{h}; 
my $index = $opts{w};
my $maxmis = $opts{g};
my $mismat = $opts{m};

my $rr = $opts{N};
my $gz = $opts{G};
my $rc = '';
my $dir = $rawdir . 'pre-clean/';
if ( -s $dir) {
  system ("rm -r $dir");   
}

if (! -s $dir) {
  rawclean ($rawdir, $cpu);   
}

my @files;
if ($opts{d} eq "all") {
  @files = <$dir*_R1.fq.gz>; 
  print "\n","OK! Now processing all data files!", "\n";
}

else {
  @files = <$dir*$opts{d}_R1.fq.gz>; 
  print "\n","OK! Now processing library $opts{d}!", "\n";
}


foreach my $file1 (@files) {
  open (IN, "gunzip -c $file1 | ");
  my $firstline = <IN>;
  #$InID = $1 if ($firstline =~ m/^\@([^:]+):/ );
  #print "can not find instrument ID!\n" if ($firstline !~ m/^\@[^:]+:/);
  #exit if ($firstline !~ m/^\@[^:]+:/);
  #close IN;
  my $seqline = <IN>;
  my $readLength = length $seqline;
  close IN;
  
  my $lib = $1 if basename($file1) =~ m/(\S+)_R1\./i;   
  
  
  ### use bfc to do error correction 
  if ($rr == 1) {
    my $bfc = 'bfc';
    my $bfc_cor = $dir . $lib . "_bfccor.fq";
    system ("$bfc -s $gz -t $cpu $file1 > $bfc_cor");
    my $r1 =  $dir . $lib . "_R1.fq1";
    
    open (IN, "<", $bfc_cor);
    open (OUT1, ">", $r1);
    
    unlink ($bfc_cor);
    while (<IN>) {
      chomp (my $line = $_);
      if ($line =~ m/^\@/) {
	chomp (my @a = split /\s+/, $line);
	chomp (my $seq1 = <IN>);
	chomp (my $quali1 = <IN>);
	chomp (my $qual1 = <IN>);
	
	print OUT1 $a[0], "\n";
	print OUT1 $seq1 , "\n";
	print OUT1 $quali1, "\n";
	print OUT1 $qual1, "\n";	  
      }
    }
    close IN;
    close OUT1;

    system ("pigz -p $cpu $r1");
    my $r1gz =  $r1. ".gz";
    system ("mv $r1gz $file1");
   
  }
  
  ## use Rcorrector for error correction
  if ($rr == 2) {
    my $run = $spath . "Rcorrector/run_rcorrector.pl"  || die "run_rcorrector.pl can not be located!\n";
    
    system ("perl $run -s $file1  -k 31 -t $cpu -od $dir");
    my $r1 =  $dir . $lib . "_R1.cor.fq.gz";
    open (IN1, "gunzip -c $r1 | ");
    
    my $R1 =  $dir . $lib . "_R1.cor.fq1";
    open (OUT1, ">", $R1);
    
    
    
    while (<IN1>) {
      chomp (my $line = $_);
      if ($line =~ m/^\@/) {
	chomp (my @a = split /\s+/, $line);
	chomp (my $seq1 = <IN1>);
	chomp (my $quali1 = <IN1>);
	chomp (my $qual1 = <IN1>);
	
	print OUT1 $a[0], "\n";
	print OUT1 $seq1 , "\n";
	print OUT1 $quali1, "\n";
	print OUT1 $qual1, "\n";
      }
    }
    close IN1;
    close OUT1;
    unlink ($r1);
    
    system ("pigz -p $cpu $R1");	
    my $r1gz =  $R1. ".gz";	
    system ("mv $r1gz $file1");
    
    
  }

  
  ##############trim adapter and merging###################
  my $start3 = time;
  my %reads = ('1' => $file1);
  my $final1;  
  #my ($outp1) = skewersingle ($dir, $lib, $file1, $P1, $P2, $P1rc,  $P2rc, '0','0','0','0', $minLength, $cpu, "1");
  my $outp1 = cutadaptsingle ($cutadapt, $dir, $lib, $file1);
  my $out1 = trimmomatic($lib, $P1, $P2, $P1rc,  $P2rc, $outp1, $outp1, 'trim1',$trimmomatic,$minLength, $cpu, $quality);
  
  $final1 = $file1 . '_cleaned1';
  my $call = system("mv $out1 $final1");
  
  
  my $time3 = int((time - $start3)/60);
  print "Trimmed in $lib in $time3 minutes! Whew, almost there!\n";
  #########################################################
  
  
  #################remove contaminaitons###################
  my $start4 = time;
  my $contaminants = $outdir . $lib . '.contam.out';
  my $suf = 'se';
  %reads = ('1' => $final1);
  removeContamination(\%reads,$contam,$contaminants,$bowtie,'se', $cpu, $mismat);
  my $time4 = int((time - $start4)/60);
  print "Removed contamination in $lib in $time4 minutes! finishing...\n";	
  #########################################################
  
  
  ######################make final files###################
  makeFinalSE ($outdir,$lib,$final1,$contaminants);
  
  my $call_rm1 = system("rm $dir*$lib*clean*");
}

if ($opts{z}) {
  print "\n","start evaluating ... OH WOW! Your sequence data are looking much better now!" , "\n\n";
  my @clean;
  
  if ($opts{d} eq "all") {
    @clean = < $outdir*_final.fq> ;
  }
  
  else {
    @clean = <$outdir$opts{d}*_final.fq>; 
  }
  
  my $resdir =  $outdir.'evaluation/';
  mkdir $resdir unless -e $resdir;
  foreach (<@clean>) {
    my $lib = $1 if basename($_) =~ m/(\S+)_[1|2|u]_final.fq/; 
    my $call1 = system("fastqc -t 2 $_ -o $resdir");
    die(qq/\nThe program "fastQC" is not in path! \n\n/) if ($call1 == -1 );
    system ("rm $resdir$lib*fastqc.zip");
  }
  
}	




sub cutadaptsingle {
  my ($cutadapt, $dir, $lib, $file1) = @_;
  my $outp1round1 = $dir . $lib . "_read1_cutadapt1.fq";
  my $outp1round2 = $dir . $lib . "_read1_cutadapt2.fq";

  system ("$cutadapt -b AGATCGGAAGAGC -n 5 -f fastq -o  $outp1round1 $file1 ");
  system ("$cutadapt -b GCTCTTCCGATCT -n 5 -f fastq -o  $outp1round2 $outp1round1");

  unlink ($outp1round1);
  return ($outp1round2);
}

sub makeFinalSE {
  my ($outdir, $lib,$reads,$contam) = @_;
  
  my %junk;
 
 # open(IN, "<$low");
 # while(<IN>){
 #   chomp(my $line = $_);
 #   $junk{$line}++;
 # }
 # close(IN);
  
  open(IN, "<$contam");
  while(<IN>) {
    chomp(my $line = $_);
    $junk{$line}++;
  }
  close(IN);
  
  my $new1 = $outdir . $lib . "_1_final.fq";
      
      open(IN, "<", $reads);
      open(OUT, ">", $new1);
      while(<IN>) {
	chomp(my $line = $_);
	my @line =split (/\t+/, $line);
	if (scalar(@line) == 1) {
	  if ($line[0] =~ m/^@(\S+)\/1$/) {
	    my $id = $1;
	    my $seq = <IN>; my $qualid = <IN>; my $qual = <IN>;
	    unless($junk{$id}){
	      print OUT "@", $id , "/1", "\n", $seq,$qualid,$qual;
	    }
	  }
	}	    
      }
      close(IN); close(OUT);
}


sub cutadatpair {
  my ($cutadapt, $dir, $lib, $file1, $file2, $type) = @_;
  my $outp1round1 = $dir . $lib . "_read1_cutadapt1.fq";
  my $outp2round1 = $dir . $lib . "_read2_cutadapt1.fq";
  my $outp1round2 = $dir . $lib . "_read1_cutadapt2.fq";
  my $outp2round2 = $dir . $lib . "_read2_cutadapt2.fq";
  if ($type eq "i") {
    system ("$cutadapt -b AGATCGGAAGAGC -B AGATCGGAAGAGC -n 5 -f fastq -o  $outp1round1 -p $outp2round1 $file1 $file2");
    system ("$cutadapt -b GCTCTTCCGATCT -B GCTCTTCCGATCT -n 5 -f fastq -o  $outp1round2 -p $outp2round2 $outp1round1 $outp2round1");
  }
  if ($type eq "n") {
    system ("$cutadapt -b AGATGTGTATAAGAGACAG -B AGATGTGTATAAGAGACAG -n 5 -f fastq -o  $outp1round1 -p $outp2round1 $file1 $file2");
    system ("$cutadapt -b CTGTCTCTTATACACATCT -B CTGTCTCTTATACACATCT -n 5 -f fastq -o  $outp1round2 -p $outp2round2 $outp1round1 $outp2round1");
  }
  
  unlink ($outp1round1, $outp2round1);
  return ($outp1round2, $outp2round2);
}


sub  extractreads {
  my ($newreads1, $rescaledBamPe, $rescaledBamSe, $ref) = @_;
  ##$reads4 = extractreads(\%newreads1, $rescaledBamPe, $rescaledBamSe,);
  my %newreads1 = %{$newreads1};
  my $read1 = $newreads1{'1'};
  my $read2 = $newreads1{'2'};
  my $readu = $newreads1{'u'};
  system (" bedtools bamtofastq -i $rescaledBamPe -fq $read1 -fq2 $read2");
  system (" bedtools bamtofastq -i $rescaledBamSe -fq $readu");
  my %reads = ('1' => $read1,'2' => $read2, 'u' => $readu); 
 
  return (\%reads);
}

sub skewer {
 my ($outdir, $lib, $file1, $file2,$P1, $P2, $minLength, $cpu) = @_;
 my $outp1 = $outdir . $lib . "_R1-trimmed-pair1.fastq";
 my $outp2 = $outdir . $lib . "_R1-trimmed-pair2.fastq";

 system ("skewer -t $cpu  -x $P1 -y $P2 $file1 $file2 -m pe -l $minLength -f sanger ");
 #system ("rm $outdir*$lib*masked*");
 return ($outp1, $outp2);
}

sub skewersingle {
  my ($outdir, $lib, $file1, $P1,$P2, $P1rc,$P2rc, $n1,$n2,$n3,$n4, $minLength, $cpu, $ps) = @_;
  my $outp1 = $outdir . $lib . "_R". $ps. "-trimmed.fastq";

  my $adfile = $file1 . ".adfile.fa";
  
  open(OUT, ">" ,$adfile);
  print OUT  ">adapter1", "\n";
  print OUT $P1, "\n";
  print OUT ">adapter2", "\n";
  print OUT $P2, "\n";
  print OUT  ">adapter3", "\n";
  print OUT $P1rc, "\n";
  print OUT ">adapter4", "\n";
  print OUT $P2rc, "\n";
  if ($n1 ne '0') {
    print OUT  ">adapter5", "\n";
    print OUT $n1, "\n";
    print OUT ">adapter6", "\n";
    print OUT $n2, "\n";
    print OUT  ">adapter7", "\n";
    print OUT $n3, "\n";
    print OUT ">adapter8", "\n";
    print OUT $n4, "\n";
  }
  close OUT;
  
  my $outp1round1 = $outdir . $lib . "_read1_cutadapt1.fq";

  system ("cutadapt -b AGATCGGAAGAGC -n 5 -f fastq -o  $outp1round1 $file1 ");
  system ("cutadapt -b GCTCTTCCGATCT -n 5 -f fastq -o  $outp1 $outp1round1");

  #system ("skewer -t $cpu  -x $adfile  $file1  -m any -l $minLength -f sanger ");
  system ("rm  $adfile");
  
  unlink ($outp1round1);
  
  return ($outp1);
}

sub makeFinal {
  my ($outdir, $lib,$reads,$contam) = @_;
  
  my %junk;
  
  #open(IN, "<$low");
  #while(<IN>){
  #  chomp(my $line = $_);
  #  $junk{$line}++;
  #}
  #close(IN);
  open(IN, "<$contam");
  while(<IN>) {
    chomp(my $line = $_);
    $junk{$line}++;
  }
  close(IN);
  
  my $new1 = $outdir . $lib . "_1_final.fq";
  my $new2 = $outdir . $lib . "_2_final.fq";
  my $newu = $outdir . $lib . "_u_final.fq";
  my %new = ('1' => $new1, '2' => $new2, 'u' => $newu);
  my %reads = %{$reads};
  foreach my $type (keys %reads) {
    open(OUT, ">$new{$type}");
    open(IN, "<$reads{$type}");
    while(<IN>) {
      chomp(my $line = $_);
      my @line =split (/\t+/, $line);
      if (scalar(@line) == 1) {
	if ($line[0] =~ m/^@(\S+)\/[1|2]$/) {
	  my $id = $1;
	  my $seq = <IN>; my $qualid = <IN>; my $qual = <IN>;
	  unless($junk{$id}){
	    print OUT $line, "\n", $seq,$qualid,$qual;
	  }
	}
      }	    
    }
    close(IN); close(OUT);
  }
}

sub getAdaptors {
  my ($lib,$adaptorFile,$libInfo,$index,$uniad) = @_;
  
  my %lib;
  open(IN, "<$libInfo");
  
  my $header = <IN>;
  
  while(<IN>) {
    chomp(my $line = $_);
    my @d = split(/\s+/,$line);
    next if $d[0] =~ /^$/;
    if ($d[2]) {
      $lib{$d[0]} = {'P7'=>$d[1], 'P5'=>$d[2]};
    }
    if (!$d[2]) {
      $lib{$d[0]} = {'P7'=>$d[1]};
    }
  }
  close(IN);
  
  my %P7;
  my %P5;
  
  
  open(IN, "<$adaptorFile");
  while(<IN>) {
    chomp(my $line = $_);
    next if $line =~ /^$/;
    if ($line =~ m/>P7_index(\d+)/) {
      my $bc = $1;
      chomp(my $seq = <IN>);
      $P7{$bc} = $seq;
    }
    if ($line =~ m/>P5_index(\d+)/) {
      my $bc2 = $1;
      chomp(my $seq2 = <IN>);
      $P5{$bc2} = $seq2;
    }
    
  } 
  close(IN);	  
  my %ad;
  if ($index == 1) {
    %ad = ("uni" => $uniad, "uni_rc" => rc($uniad), "index" => $P7{$lib{$lib}{'P7'}}, "index_rc" => rc($P7{$lib{$lib}{'P7'}}));
    die(qq/\nHmm...I could not find the adapter sequences for $lib ... Check the naming of the adapter sequences in "-a". \n\n/) if (!$ad{"index"});
  }
  if ($index == 2) {
    %ad = ("uni" =>$P5{$lib{$lib}{'P5'}}, "uni_rc" => rc($P5{$lib{$lib}{'P5'}}), "index" => $P7{$lib{$lib}{'P7'}}, "index_rc" => rc($P7{$lib{$lib}{'P7'}}));
    die(qq/\nHmm...I could not find the adapter sequences for $lib ... Check the naming of the adapter sequences in "-a". \n\n/) if (!$ad{"index"} || !$ad{"uni"});
    
  } 
  
  return(\%ad);
}

sub rc {
  my ($seq) = @_;
  my $rc = $seq;
  $rc = reverse($rc);
  $rc =~ tr/ATGCatgc/TACGtacg/;
  return($rc);
}

sub removeContamination {
  my ($reads,$contam,$contaminants, $bowtie, $suf, $cpu,  $mismat) = @_;
  my %reads = %{$reads};
  unless (-f $contam . ".3.bt2") {
    my $bw2build = $bowtie . "-build";
    my $call1 = system("$bw2build $contam $contam");
  }
  if ($suf eq 'pe')  {    
    my $contamout1 = $reads{'1'} . ".contam.sam1";
    my $call2 = system("$bowtie -x $contam -1 $reads{'1'} -2 $reads{'2'} --fast -S $contamout1 -p $cpu --sam-nohead --sam-nosq");
    die(qq/\nThe program "bowtie2" is not in path! \n\n/) if ($call2 == -1 );
    my $contamout2 = $reads{'1'} . ".contam.sam2";
    my $call3 = system("$bowtie -x $contam -U $reads{'u'} --fast -S $contamout2 -p $cpu --sam-nohead --sam-nosq");
    my $contamout_all = $reads{'1'} . ".contam_all.sam";
    system ("cat $contamout1 $contamout2 > $contamout_all");
   
    parseSAM($contaminants,$contamout_all, $mismat);
    system ("rm $contamout1 $contamout2");  
  }
  if ($suf eq 'se') {
    
    my $contamout1 = $reads{'1'} . "contam.sam1";  
    my $call3 = system("$bowtie -x $contam -U $reads{'1'} --fast -S $contamout1 -p $cpu --sam-nohead --sam-nosq");
    parseSAM ($contaminants,$contamout1, $mismat);
    system ("rm $contamout1 ");
  }
  
}
#1H56M1D43M

sub parseSAM {
  my ($contaminants,$contam,  $mismat) = @_;
  open(OUT, ">$contaminants");
  open(IN, "<$contam");
  
  while(<IN>) {
    chomp(my $line = $_);
    my @d = split(/\t/,$line);
    if ($d[2] !~ m/\*/) {
      if ($d[5] =~ m/\d+M/) {
	my @sum = ($d[5] =~ m/(\d+)/);
	my $length = sum (@sum);
	my $md = $1 if $line =~ m/(MD\:Z\S+)/;
	my @a = ($md =~ m/([ATGC])/g);
	if (scalar(@a)/$length < 1 - $mismat) {
	  $d[0] =~ s/\/[1|2]$//;
	  print OUT $d[0], "\n";
	}
      }
    }
  }
  close(IN); close(OUT);  
  #unlink($contam);
}       


sub fixMatePair1 {
  my ($lib,$read,$readarray,$base) = @_;
  my @trim = @{$readarray};
  my %pair;	
  foreach my $reads (@trim) {
    open(IN, "<$reads");
    while(<IN>) {
      chomp(my $line = $_);
      if ($line =~ m/^@(\S+)\/[1|2]$/) {
	$pair{$1}++;
	chomp(my $seq = <IN>); chomp(my $qualid = <IN>); chomp(my $qual = <IN>);
      }
    }
    close(IN);	
  }
  my %reads = %{$read};
  my $out1 = $reads{'1'} . '_' . $base . '_p1';
  my $out2 = $reads{'2'} . '_' . $base . '_p2';
  my $outu = $reads{'1'} . '_' . $base . '_u';
  print $out1, "\n"; 
  print $out2, "\n";
  print $outu, "\n";
  open(OUT1, ">$out1"); 
  open(OUT2, ">$out2"); 
  open(OUTU, ">$outu"); 
  my %newpairs = ('1' => $out1, '2' => $out2, 'u' => $outu);
  foreach my $reads (@trim) {
    open(IN, "<$reads");
    my $file = $1 if $reads =~ m/_R(\d+)\./;
    while(<IN>) {
      chomp(my $line = $_);	
      if ($line =~ m/^@(\S+)\/[1|2]$/) {
	my $id = $1;
	my $seq = <IN>;
	my $qualid = <IN>;
	my $qual = <IN>;
	if ($pair{$id} == 2) {
	  if ($file == 1) {
	    print OUT1 $line . "\n" . $seq . $qualid . $qual;
	  }
	  else {
	    print OUT2 $line . "\n" . $seq . $qualid . $qual;
	  }
	}
	else {
	  print OUTU $line . "\n" . $seq . $qualid . $qual;
	}
      }	
    }
    close(IN);	
  }	
  close(OUT1); close(OUT2); close(OUTU);	
  return(\%newpairs);
}

sub trimmomatic {
  my ($lib,$P1,$P2, $P1rc,$P2rc, $in,$base,$suffix, $trimmomatic,$minLength, $cpu, $quality) = @_;
  my $out = $base . '_' . $suffix;	
  my $adfile = $in . "_adfile.fa";
  open(OUT, ">" ,$adfile);
  print OUT  ">single1", "\n";
  print OUT $P1, "\n";
  print OUT ">single2", "\n";
  print OUT $P2, "\n";
  print OUT ">single3", "\n";
  print OUT $P1rc, "\n";
  print OUT ">single4", "\n";
  print OUT $P2rc, "\n";
  
  close OUT;
  
  my $call1 = system(" $trimmomatic -Xms512m -Xmx4g  SE -threads $cpu  -phred33 $in $out ILLUMINACLIP:$adfile:2:30:10 SLIDINGWINDOW:4:$quality MINLEN:$minLength LEADING:3 TRAILING:3");
  unlink($adfile);
  return($out);
}

sub cutadapt {
  my ($ad,$in,$base,$suffix, $cutadapt,$minLength) = @_;
  my $out  = $base . '_' . $suffix;	
  my $curRead = $in;
  my $tracker = 1;
  my %ad = %{$ad};
  foreach my $key (keys %ad) {
    my $out = $curRead . $tracker;
    my $call = system("$cutadapt -b $ad{$key} -O 4 -n 5 -e 0.15 -f fastq $curRead -o $out -m $minLength");
    die(qq/\nThe program "cutadapt" is not in path! \n\n/) if ($call == -1 );
    unlink($curRead) unless($curRead eq $in);
    $curRead = $out;
    $tracker++;
  }
  my $call2 = system("mv $curRead $out");
  return($out);
}

sub bowtie {
  my ($lib, $ad,$in,$base,$suffix,$bowtie,$minLength, $cpu) = @_;
  my $out  = $base . '_' . $suffix;	
  my $file = $lib."out.sam";
  
  my $adfile = $lib . "_adfile_norev.fa";
  open(OUT, ">$adfile");
  foreach my $name (keys %{$ad}) {
    print OUT ">", $name, "\n", $ad->{$name}, "\n" unless $name =~ m/rc/;
  }
  close(OUT);	
  
  my $bw2build = $bowtie . "-build";
  my $call1 = system("$bw2build $adfile $adfile");
  my $call2 = system("$bowtie --local -D 15 -R 2 -N 1 -L 10 -i S,1,0.75 -p $cpu -k 1 -x $adfile -U $in -S $file");
  die(qq/\nThe program "bowtie2" is not in path! \n\n/) if ($call2 == -1 );
  my $call3 = system("rm $adfile" . "*");
  
  open(IN, "<$file");
  open(OUT, ">$out");
  
  while(<IN>) {
    chomp(my $line1 = $_);
    my @d1 = split(/\t/,$line1);
    
    my $seq1 = $d1[9];
    my $qual1 = $d1[10];
    
    unless($line1 =~ m/^@/) {
      if ($line1 !~ m/\*/) {			
	if ($d1[5] =~ m/^(\d+)S\d+M$/) {
	  my $l = $1;
	  $seq1 = substr $seq1, 0, $l;
	  $qual1 = substr $qual1, 0, $l;
	}
	elsif ($d1[5] =~ m/^(\d+)M\d+S$/) {
	  my $start = $1;
	  $seq1 = substr $seq1, $start;
	  $qual1 = substr $qual1, $start;
	}
	else {
	  my @s;
	  while ($d1[5] =~ m/(\d+)S/g) {
	    push(@s,$1);
	  }
	  @s = sort {$a <=> $b} @s;
	  if ($s[$#s] >= $minLength) {	
	    if ($d1[5] =~ m/^(\S*)$s[$#s]/) {
	      my $match = $1;
	      my $length = $s[$#s];
	      my $start = 0;
	      while ($match =~ m/(\d+)/g) {
		$start += $1;
	      }
	      $seq1 = substr $seq1, $length;
	      $qual1 = substr $qual1, $length;	
	    }													
	  }
	  else {
	    $seq1 = 'N'; $qual1 = 'N';
	  }
	}
      }
      
      if (length($seq1) >= $minLength) {	
	print OUT "@" . $d1[0] . "\n" . $seq1 . "\n" . '+' . "\n" . $qual1 . "\n";	
      }
    }	
  }	
  unlink($file);
  close(IN); close(OUT);	
  return($out);
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



#now convert the raw data files to format compatible to downstream analyses
sub rawclean {
  my ($rawdir, $cpu) = @_;
  my $Result_dir1 = $rawdir . 'combined/';
  mkdir $Result_dir1 unless -e  $Result_dir1;
  my $dir = $rawdir . 'pre-clean/';
  mkdir $dir unless -e  $dir;
  
  my @orig_files = <$rawdir*.gz>; 
  
  die(qq/\nHmm...Is it the right folder? \n\n/) if (scalar (@orig_files) == 0);
  my %lib;
  
  foreach (<@orig_files>) {   
    my $file = basename($_) =~ /(\S+)_L\d+_(R[1|2])_\d+\.fastq\.gz/;
    my $libs = $1. $2;
    $lib{$libs}++; 
  }
  
  @orig_files = <$rawdir*.gz>;

 foreach (<@orig_files>) {   
   my $file = basename($_) =~ /(\S+)_L\d+_(R[1|2])_\d+\.fastq\.gz/;
   my $lib_name = $1;
   my $d = $2;
   my $libs = $lib_name . $d;
   if ($lib{$libs}) {
     
     my $file2 = $rawdir . $lib_name .  '*'.  $d . '*' . '.fastq.gz';
     my $combined = $Result_dir1 . $lib_name . "_". $d . '.fastq.gz'; 

     system ("cat $file2 > $combined ");
     
   }
   delete $lib{$libs} if $lib{$libs};
 }
  
  my @merged_files = < $Result_dir1*.fastq.gz> ;
  
  foreach my $file (<@merged_files>) {
    my $out = $dir .  $1 .".fq.gz" if basename($file)  =~ /(\S+_R[1|2])\.fastq\.gz/;
    my $redundancy = '^--$' ; 
    print "cleaning","\t",$file,"\n";
    if ($file =~ m/_R([1|2])\.fastq\.gz/) {
      if ($1 == 1) {
	system ("gunzip -c $file | grep -A 3 \'^@.* [^:]*:N:[^:]*:\'  | grep -v $redundancy | sed \'s/ 1:N:0:.*/\\/1/g\' |  pigz -p $cpu > $out");
      }
      if ($1 == 2) {
	system ("gunzip -c $file | grep -A 3 \'^@.* [^:]*:N:[^:]*:\' | grep -v $redundancy | sed  \'s/ 2:N:0:.*/\\/2/g\' |  pigz -p $cpu > $out");
      }  
    }
  }
  system ("rm -r $Result_dir1");
  return ($dir);
}





