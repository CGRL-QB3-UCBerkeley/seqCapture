#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];
use List::Util qw(sum);
use File::Temp;

&main;
exit;


sub main {
  &usage if (@ARGV<1);
  my $command = shift(@ARGV);
  my %p = (pop=>\&pop, phylo=>\&phylo);
  die("Unknown command \"$command\"\n") if (!defined($p{$command}));
  &{$p{$command}};
}

sub usage {
  die(qq/
Usage: seqCapture evaluation <command> [<arguments>]\n


Command: 

phylo:   evaluation for phylogenetic datasets

pop:     evaluation for population genetic datasets

\n/);
}


sub pop {

die(qq/

Usage: seqCapture evaluation pop  [options]

Options: 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-genome       FILE    Reference genome; Required 
 
-cleanDir     DIR     Path to a folder with all cleaned reads (XXX_1_final.fq, 
                      XXX_2_final.fq, XXX_u_final.fq; Required

-rawDir       DIR     Path to a folder with all raw reads (XXX_R1.fq.gz, XXX_R2.fq.gz); Required

-bamDir       DIR     Path to a folder with all sorted bam files (XXX_sorted.bam); Required

-resDir       DIR     Path to the results directory; Required

-bedFile      FILE    A bed file that contain regions that probes are tiled
                      skip this option if no bed file 
                    
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                       


\n/) if (!@ARGV);

my ( $genome, $cleanDir, $rawDir,$bamDir, $resDir, $bedFile,$window) = (undef, undef, undef, undef, undef,undef,undef,undef);

#my %opts = (readLen => [100]);
#'readLen=s@{1,1}' => \$opts{readLen}
GetOptions('genome=s@{1,1}' => \$genome,'cleanDir=s@{1,1}' => \$cleanDir,'rawDir=s@{1,1}' => \$rawDir, 'bamDir=s@{1,1}' => \$bamDir, 'bedFile=s@{1,1}' => \$bedFile, 'resDir=s@{1,1}' => \$resDir,  'window=s@{1,1}' => \$window);

my $len;
if (@{$genome}[0]){
  my $file = @{$genome}[0];
  open (IN, "<", $file);    
  while (<IN>) {
    chomp (my $line = $_);
    if ($line !~ m /^>/) {
      $line =~ s/N//g;
      $len += length ($line);
    }
  }
  print "The total number of unmaksed bases in this reference genome is ", $len, ".\n";
}

my $rawD = redir (@{$rawDir}[0]);
my $cleanD = redir (@{$cleanDir}[0]);
my $bamD = redir (@{$bamDir}[0]);
my $resD = redir (@{$resDir}[0]);
mkdir $resD unless -e $resD;
my $out = $resD . "evaluationResults.txt";

open (OUT, ">", $out);
print OUT "Library", "\t", "rawData(Mb)", "\t", "cleanedData(Mb)", "\t", "DataMapped(Mb)", "\t","Specificity(%)", "\t", "Sensitivity(%)", "\t", "AvgCoverage(X)",  "\t","CoverageVar",  "\t",  "AvgCoverageStd", "\t", "Sites2X(%)","\t", "Sites5X(%)","\t", "Sites10X(%)","\t", "Sites20X(%)", " \n";
my @bam = <$bamD*sorted.bam>;
my @raw = <$rawD*_R1.fq.gz>;
my @clean = <$cleanD*_1_final.fq>;

my $mas1;
$mas1 = data (\@raw, \@clean);

my $mas2;
my $mas3;
my $mas4;

if (@{$bedFile}[0]) {
  $mas2 = sensitivitypop (\@bam, @{$bedFile}[0],$len);
  $mas3 = specificitypop (\@bam, @{$bedFile}[0]);
  open (IN, "<", @{$bedFile}[0]);
  my $totalSite;
  while (<IN>) {
    chomp (my @line = split /\s+/, $_);
    if ($line[2] && $line[1]) {
      $totalSite += $line[2] - $line[1];
    }
  }
  #print "The total number of targeted bp is ", $totalSite , ".\n";
}

if (!@{$bedFile}[0]) {
  $mas2 = sensitivitypop (\@bam, 1,$len);  
  $mas3 = specificitypop (\@bam, 1);
}

my %mas1 = %{$mas1};
my %mas2 = %{$mas2};
my %mas3 = %{$mas3};

foreach my $id (sort {$a cmp $b} keys %mas1) {
  print OUT $id, "\t", $mas1{$id}{'raw'}, "\t", $mas1{$id}{'clean'},"\t", $mas3{$id}{'mapped'},"\t", sprintf("%.2f", $mas3{$id}{'mapped'}/$mas1{$id}{'clean'}*100), "\t", $mas2{$id}{'sen'}, "\t", $mas2{$id}{'cov'}, "\t", $mas2{$id}{'covvar'},"\t", $mas2{$id}{'std'}, "\t", $mas2{$id}{'2'}, "\t", $mas2{$id}{'5'},"\t", $mas2{$id}{'10'},"\t", $mas2{$id}{'20'},"\n";
  
}
close OUT;

}

sub phylo {

die(qq/

Usage: seqCapture evaluation phylo  [options]

Options: 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-genomeDir    DIR    Rerence genome Dir (XXX_targetedRegionAndFlanking.fasta)

-cleanDir     DIR     Path to a folder with all cleaned reads (XXX_1_final.fq, 
                      XXX_2_final.fq, XXX_u_final.fq)

-rawDir       DIR     Path to a folder with all raw reads (XXX_R1.fq.gz, XXX_R2.fq.gz)

-bamDir       DIR     Path to a folder with all sorted bam files (XXX_sorted.bam)

-resDir       DIR     Path to the results directory

-bedDir       DIR     Path to a folder with all bed files (XXX_targetedRegionforExonCapEval.bed)
                     
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                       


\n/) if (!@ARGV);
  
  my ($genomeDir, $cleanDir, $rawDir, $bamDir, $resDir, $bedDir,$window) = (undef, undef, undef, undef, undef,undef,undef);
  GetOptions('genomeDir=s@{1,1}' => \$genomeDir,'cleanDir=s@{1,1}' => \$cleanDir,'rawDir=s@{1,1}' => \$rawDir, 'bamDir=s@{1,1}' => \$bamDir, 'bedDir=s@{1,1}' => \$bedDir, 'resDir=s@{1,1}' => \$resDir,  'window=s@{1,1}' => \$window);

my $genomeD = redir (@{$genomeDir}[0]);
my $bedD = redir (@{$bedDir}[0]);
my $rawD = redir (@{$rawDir}[0]);
my $cleanD = redir (@{$cleanDir}[0]);
my $bamD = redir (@{$bamDir}[0]);
my $resD = redir (@{$resDir}[0]);
mkdir $resD unless -e $resD;
my $out = $resD . "evaluationResults.txt";
open (OUT, ">", $out);
print OUT "Library", "\t", "rawData(Mb)", "\t", "cleanedData(Mb)", "\t", "DataMapped(Mb)", "\t","Specificity(%)", "\t", "Sensitivity(%)", "\t", "AvgCoverage(X)",  "\t","CoverageVar",  "\t",  "AvgCoverageStd", "\t", "Sites2X(%)","\t", "Sites5X(%)","\t", "Sites10X(%)","\t", "Sites20X(%)", " \n";
my @bam = <$bamD*sorted.bam>;
my @raw = <$rawD*_R1.fq.gz>;
my @clean = <$cleanD*_1_final.*>;
my @bed = <$bedD*targetedRegionforExonCapEval.bed> if $bedD;
my @genome = <$genomeD*targetedRegionAndFlanking.fasta>;

my $mas1;

$mas1 = data (\@raw, \@clean);
my %mas1 = %{$mas1};
my %mas2;
my %mas3;


foreach (@genome) {
  my $ref = $_;
  my $bed;
  if ($bedD) {
    $bed = $bedD . $1. "_targetedRegionforExonCapEval.bed" if basename($ref) =~ m/(\S+)_targetedRegionAndFlanking\.fasta/;
  }
  
  my $bam = $bamD . $1 . "_sorted.bam"  if basename ($ref) =~ m/(\S+)_targetedRegionAndFlanking\.fasta/;
  my $mas2;
  my $mas3;
  
  my $lib = $1 if  basename ($ref) =~ m/(\S+)_targetedRegionAndFlanking\.fasta/;
  
  
  my $len;
  if ($ref){
    my $file = $ref;
    open (IN, "<", $file);    
    while (<IN>) {
      chomp (my $line = $_);
      if ($line !~ m /^>/) {
	$line =~ s/N//g;
	$len += length ($line);
      }
    }
    print "The total number of unmaksed bases in $lib reference is ", $len, ".\n";
  }
  
  if ($bed) {
    
    open (IN, "<", $bed);
    my $totalSite;
    while (<IN>) {
      chomp (my @line = split /\s+/, $_);
      
      $totalSite += $line[2] - $line[1];
      
      
    }
    print "The total number of targeted bp of $lib is ", $totalSite , ".\n";
    
    
    $mas2 = sensitivity ($bam, $bed,$len, \%mas1);
    %mas2 = %{$mas2};
    $mas3 = specificity ($bam, $bed,\%mas2);
    %mas3 = %{$mas3};
    
    
  } ##if ($bed) {
  
  if (!$bed) {
    $mas2 = sensitivity ($bam, 1,$len,\%mas1 );
    %mas2 = %{$mas2};
    $mas3 = specificity ($bam, 1,\%mas2);
    %mas3 = %{$mas3};
  } ##if (!$bed) {
   
} ## foreach (@genome) {



foreach my $id (sort {$a cmp $b} keys %mas3) {
  print OUT $id, "\t", $mas3{$id}{'raw'}, "\t", $mas3{$id}{'clean'},"\t", $mas3{$id}{'mapped'},"\t", sprintf("%.2f", $mas3{$id}{'mapped'}/$mas3{$id}{'clean'}*100), "\t", $mas3{$id}{'sen'}, "\t", $mas3{$id}{'cov'}, "\t", $mas3{$id}{'covvar'},"\t", $mas3{$id}{'std'}, "\t", $mas3{$id}{'2'}, "\t", $mas3{$id}{'5'},"\t", $mas3{$id}{'10'},"\t", $mas3{$id}{'20'},"\n";
  
}
close OUT;

}


sub data {
  my ($raw, $clean)= @_;    
  my @raw = @{$raw};
  my @clean = @{$clean};
  my %master;   
  foreach (@raw) {
    my $lib = $1 if basename($_) =~ m/(\S+)_R1\.fq\.gz/;
    print "\nNow calculating the size of raw data for ", $lib, "!\n";
    
    open (IN, "gunzip -c $_ | ");
    my $first = <IN>;
    chomp (my $seq1 = <IN>);
    my $length = length ($seq1);
    seek IN, 0, 0;
    
    my $count=0;
    while (<IN>) {
      chomp (my $line = $_);
      if ($line =~ m /^\@/) {
	chomp(my $seq = <IN>);
	chomp(my $qualID = <IN>);
	chomp(my $qual = <IN>);
	$count++;
      }
    }
    close IN;
    my $result = sprintf("%.2f",2*$count*$length/1000000);
    $master{$lib}{'raw'} = $result; 
  }
  
  foreach my $file1 (@clean) {
    my $file2 = $file1;
    $file2 =~ s/_1_final/_2_final/;
    my $fileu = $file1;
    $fileu =~ s/_1_final/_u_final/;
    my $lib = $1 if basename($file1) =~ m/(\S+)_1_final\./i;   
    print "\nNow calculating the size of cleaned data for ", $lib, "!\n";
    my $total_length = tot_length($file1) + tot_length($file2) + tot_length($fileu); 
    my $result = sprintf("%.2f",$total_length/1000000);
    $master{$lib}{'clean'} = $result;
  }
  
  sub tot_length {
    my ($file) = @_;
    open (IN, "<", $file);
    my $totlength = 0;
    while (<IN>) {
      chomp(my $line = $_);
      if ($line =~ m/^\@/) {
	chomp(my $seq = <IN>);
	chomp(my $qualID = <IN>);
	chomp(my $qual = <IN>);
	$totlength += length($seq);
      }
    }
    close IN;
    return ($totlength);
  }
  return (\%master);    
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


sub sensitivitypop {
  my ($bam, $bed, $len) = @_;
  my %master;
  
  my @all = @{$bam};
  foreach (@all) {
    my $sites = 0;
    my $cov = 0;
    my $expxsqr = 0;
    my $file = $_;
    my $one = 0;
    my $five = 0;
    my $ten = 0;
    my $twenty = 0;
    
    my $lib = $1 if basename($file) =~ m/(\S+)_sorted.bam/;
    print "\nNow calculating sensitivity and coverage for ", $lib, "!\n";
    if ($bed =~ m /\b1\b/) {      
      my @call = `samtools depth -q 0 $file`;
      foreach (@call) {
	$sites++;
	chomp (my @line = split /\s+/, $_);
	if ($line[2] >= 2 && $line[2] < 5 ) {
	  $one++;	  
	}
	if ($line[2] >= 5 && $line[2] < 10 ) {
	  $one++;$five++;	  
	}
	if ($line[2] >= 10 && $line[2] < 20 ) {
	  $one++;$five++;$ten++;	  
	}
	if ($line[2] >= 20 ) {
	  $one++;$five++;$ten++;$twenty++;	  
	}
	$cov += $line[2];
	$expxsqr = $expxsqr + ($line[2] ** 2);
      }
      my $result =  sprintf("%.2f",$sites/$len*100);
      my $final_coverage = sprintf("%.2f", $cov/$len);
      my $final1 =  sprintf("%.2f",$one/$len*100);
      my $final5 =  sprintf("%.2f",$five/$len*100);
      my $final10 =  sprintf("%.2f",$ten/$len*100);
      my $final20 =  sprintf("%.2f",$twenty/$len*100);   
      my $varq = ($expxsqr - (($cov ** 2) / $sites)) / ($sites - 1); 
      my $std = sqrt ($varq);
      $master{$lib}{'sen'} = $result; 
      $master{$lib}{'cov'} = $final_coverage;	
      $master{$lib}{'covvar'} =  sprintf("%.2f",$varq);
      $master{$lib}{'std'} =  sprintf("%.2f",$std);
      $master{$lib}{'2'} =  sprintf("%.2f",$final1);
      $master{$lib}{'5'} =  sprintf("%.2f",$final5);
      $master{$lib}{'10'} =  sprintf("%.2f",$final10);
      $master{$lib}{'20'} =  sprintf("%.2f",$final20);
    } ###if ($bed =~ m /\b1\b/) { 
    
    else {
      my @call = `samtools depth  -b $bed -q 0 $file`;
      foreach (@call) {
	$sites++;
	chomp (my @line = split /\s+/, $_);
	if ($line[2] >= 2 && $line[2] < 5 ) {
	  $one++;	  
	}
	if ($line[2] >= 5 && $line[2] < 10 ) {
	  $one++;$five++;	  
	}
	if ($line[2] >= 10 && $line[2] < 20 ) {
	  $one++;$five++;$ten++;	  
	}
	if ($line[2] >= 20 ) {
	  $one++;$five++;$ten++;$twenty++;	  
	}
	$cov += $line[2];
	$expxsqr = $expxsqr + ($line[2] ** 2);
      }
      my $totalS;
      open (IN, "<", $bed);
      while (<IN>) {
	chomp (my @line = split /\s+/, $_);	  
	$totalS += $line[2] - $line[1];
      }
      close IN;
      my $result =  sprintf("%.2f",$sites/$totalS*100);
      my $final_coverage = sprintf("%.2f",$cov/$totalS);
      my $final1 =  sprintf("%.2f",$one/$totalS*100);
      my $final5 =  sprintf("%.2f",$five/$totalS*100);
      my $final10 =  sprintf("%.2f",$ten/$totalS*100);
      my $final20 =  sprintf("%.2f",$twenty/$totalS*100);   
      my $varq = ($expxsqr - (($cov ** 2) / $totalS)) / ($totalS - 1); 
      my $std = sqrt ($varq);
      $master{$lib}{'sen'} = $result; 
      $master{$lib}{'cov'} = $final_coverage;
      $master{$lib}{'covvar'} =  sprintf("%.2f",$varq);
      $master{$lib}{'std'} =  sprintf("%.2f",$std);
      $master{$lib}{'2'} =  sprintf("%.2f",$final1);
      $master{$lib}{'5'} =  sprintf("%.2f",$final5);
      $master{$lib}{'10'} =  sprintf("%.2f",$final10);
      $master{$lib}{'20'} =  sprintf("%.2f",$final20);	
    } ##else
    
    
  } ##foreach (@call)
  return (\%master);
}




sub specificitypop {
  my ($bam, $bed) = @_;
  my @bam = @{$bam}; 
  my %master;
  foreach (@bam) {
    
    my $file = $_;
    my $id = $1 if basename($file) =~ m /(\S+)_sorted.bam/;
    
    print "\nNow calculating specificity for ", $id, "!\n";
    
    if ($bed =~ m /\b1\b/) { 
      my @call = `samtools view $file`;
      my $data;
      foreach (@call) {
	my @line = split /\s+/, $_;
	my $d = $line[5];
	next if $line[5] =~ /\*/;
	my @a = ($d =~ m/(\d+)M/g);
	my $sum = sum(@a);
	$data += $sum;
      }
      my $final = sprintf("%.2f", $data/1000000);
      $master{$id}{'mapped'} = $final;   
    }
    else {
      my @call = `samtools view -L $bed $file`;
      my $data;
      foreach (@call) {
	my @line = split /\s+/, $_;
	my $d = $line[5];
	next if $line[5] =~ /\*/;
	my @a = ($d =~ m/(\d+)M/g);
	my $sum = sum(@a);
	$data += $sum;
      }
      my $final = sprintf("%.2f", $data/1000000);
      $master{$id}{'mapped'} = $final;   	
    }
  }
  return (\%master);
}

sub specificityphylo {
  my ($bam, $bed, $master) = @_;  
  my %master = %{$master};
  my $id = $1 if basename($bam) =~ m /(\S+)_sorted.bam/;
  
  print "\nNow calculating specificity for ", $id, "!\n";
  
  if ($bed =~ m /\b1\b/) { 
    my @call = `samtools view $bam`;
    my $data;
    foreach (@call) {
      my @line = split /\s+/, $_;
      my $d = $line[5];
      next if $line[5] =~ /\*/;
      my @a = ($d =~ m/(\d+)M/g);
      my $sum = sum(@a);
      $data += $sum;
    }
    my $final = sprintf("%.2f", $data/1000000);
    $master{$id}{'mapped'} = $final;   
  }
  else {
    my @call = `samtools view -L $bed $bam`;
    my $data;
    foreach (@call) {
      my @line = split /\s+/, $_;
      my $d = $line[5];
      next if $line[5] =~ /\*/;
      my @a = ($d =~ m/(\d+)M/g);
      my $sum = sum(@a);
      $data += $sum;
    }
    my $final = sprintf("%.2f", $data/1000000);
    $master{$id}{'mapped'} = $final;   	
  }
  
  return (\%master);
}


sub sensitivityphylo {
  my ($bam, $bed, $len, $master) = @_;
  my %master = %{$master};
  my $sites = 0;
  my $cov = 0;
  my $expxsqr = 0;
  my $file = $_;
  my $one = 0;
  my $five = 0;
  my $ten = 0;
  my $twenty = 0;
  
  my $lib = $1 if basename($bam) =~ m/(\S+)_sorted.bam/;
  print "\nNow calculating sensitivity and coverage for ", $lib, "!\n";
  if ($bed =~ m /\b1\b/) {      
    my @call = `samtools depth -q 0 $bam`;
    foreach (@call) {
      $sites++;
      chomp (my @line = split /\s+/, $_);
      if ($line[2] >= 2 && $line[2] < 5 ) {
	$one++;	  
      }
      if ($line[2] >= 5 && $line[2] < 10 ) {
	$one++;$five++;	  
      }
      if ($line[2] >= 10 && $line[2] < 20 ) {
	$one++;$five++;$ten++;	  
      }
      if ($line[2] >= 20 ) {
	$one++;$five++;$ten++;$twenty++;	  
      }
      $cov += $line[2];
      $expxsqr = $expxsqr + ($line[2] ** 2);
    }
    my $result =  sprintf("%.2f",$sites/$len*100);
    my $final_coverage = sprintf("%.2f", $cov/$len);
    my $final1 =  sprintf("%.2f",$one/$len*100);
    my $final5 =  sprintf("%.2f",$five/$len*100);
    my $final10 =  sprintf("%.2f",$ten/$len*100);
    my $final20 =  sprintf("%.2f",$twenty/$len*100);   
    my $varq = ($expxsqr - (($cov ** 2) / $sites)) / ($sites - 1); 
    my $std = sqrt ($varq);
    $master{$lib}{'sen'} = $result; 
    $master{$lib}{'cov'} = $final_coverage;	
    $master{$lib}{'covvar'} =  sprintf("%.2f",$varq);
    $master{$lib}{'std'} =  sprintf("%.2f",$std);
    $master{$lib}{'2'} =  sprintf("%.2f",$final1);
    $master{$lib}{'5'} =  sprintf("%.2f",$final5);
    $master{$lib}{'10'} =  sprintf("%.2f",$final10);
    $master{$lib}{'20'} =  sprintf("%.2f",$final20);
  }
  
  else {
    my @call = `samtools depth -b $bed -q 0 $bam`;
    foreach (@call) {
      $sites++;  
      chomp (my @line = split /\s+/, $_);
      if ($line[2] >= 2 && $line[2] < 5 ) {
	$one++;	  
      }
      if ($line[2] >= 5 && $line[2] < 10 ) {
	$one++;$five++;	  
      }
      if ($line[2] >= 10 && $line[2] < 20 ) {
	$one++;$five++;$ten++;	  
      }
      if ($line[2] >= 20 ) {
	$one++;$five++;$ten++;$twenty++;	  
      }
      $cov += $line[2];
	$expxsqr = $expxsqr + ($line[2] ** 2);
    }	
    my $totalS;
    open (IN, "<", $bed);
    while (<IN>) {
      chomp (my @line = split /\s+/, $_);	
      $totalS += $line[2] - $line[1];	
    }
    close IN;
    my $result =  sprintf("%.2f",$sites/$totalS*100);
    my $final_coverage = sprintf("%.2f",$cov/$totalS);
    my $final1 =  sprintf("%.2f",$one/$totalS*100);
    my $final5 =  sprintf("%.2f",$five/$totalS*100);
    my $final10 =  sprintf("%.2f",$ten/$totalS*100);
    my $final20 =  sprintf("%.2f",$twenty/$totalS*100);   
    my $varq = ($expxsqr - (($cov ** 2) / $totalS)) / ($totalS - 1); 
    my $std = sqrt ($varq);
    $master{$lib}{'sen'} = $result; 
    $master{$lib}{'cov'} = $final_coverage;
    $master{$lib}{'covvar'} =  sprintf("%.2f",$varq);
    $master{$lib}{'std'} =  sprintf("%.2f",$std);
    $master{$lib}{'2'} =  sprintf("%.2f",$final1);
    $master{$lib}{'5'} =  sprintf("%.2f",$final5);
    $master{$lib}{'10'} =  sprintf("%.2f",$final10);
    $master{$lib}{'20'} =  sprintf("%.2f",$final20);	
  }    
  return (\%master);
}
