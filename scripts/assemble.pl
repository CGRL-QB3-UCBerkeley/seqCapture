#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;
use File::Basename;

die(qq/

Usage: seqCapture assemble  [options]


Options: 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-reads    DIR             Directory with all sequence reads
-kmer     INT,INT,INT...  Kmer lengths chosen for SPAde assemblies
                          [21,33,55] (no space)
-lib      CHAR ...        Particular libraries to process? 
                          (e.g. AAA BBB CCC). If -lib is not 
                          used then process all libraries in
                          The folder (-reads)  
-out      CHAR            Directory where results will go
-np       INT             number of processors used for assembly
         
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
\n/) if (!@ARGV);


my ($reads, $out, $np, $kmer, $lib) = (undef, undef, undef, undef, undef);
GetOptions('reads=s@{1,1}' => \$reads,'out=s@{1,1}' => \$out, 'np=s@{1,1}' => \$np, 'lib=s@{,}' => \$lib ,'kmer=s@{,}' => \$kmer);

my $spades = "spades.py";

my $dir = redir (@{$reads}[0]);
my $resDir = redir (@{$out}[0]) ;
mkdir $resDir unless -d $resDir;

my @files = ();
if (!$lib) {
  @files= <$dir*_1_final.fq>;
}

if ($lib) { 
  foreach (@{$lib}) {
    my $file = $dir . $_ .'_1_final.fq';	
      push (@files, $file);      
  }
} 

foreach my $file1 (@files) {
  my $file2 = $file1;
  my $fileu = $file1;
  $file2 =~ s/_1_/_2_/;
  $fileu =~ s/_1_/_u_/;
  my $lib = $1 if  basename($file1) =~ m/(\S+)_1_final.fq/;
  
  my $resultDir = $resDir. $lib. "/";
  mkdir $resultDir unless -d $resultDir;
  
  if ($kmer) {
      system ("$spades -1 $file1 -2 $file2 -s $fileu -k @{$kmer} -t @{$np}[0] -o $resultDir");
  }
  else {
      system ("$spades -1 $file1 -2 $file2 -s $fileu  -t @{$np}[0] -o $resultDir");
  }
  my $raw = $resultDir . "contigs.fasta";
  my $raw2 = $resDir . $lib. ".fasta";
  system ("mv $raw $raw2");
  system ("rm -r $resultDir");
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
