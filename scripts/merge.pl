#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;

die(qq/

Usage: seqCapture mergealign  [options]

-f      DIR         a folder with all subfolders, 
                    and the subfolders are named 
                    as sub1, sub2,... subn;
-o      DIR         result folder 

\n\n/) unless (@ARGV);

my %opts = (f=>undef, o=>undef);
getopts('f:o:', \%opts);

my $fadir = redir ($opts{f});
my $outdir = redir ($opts{o});
mkdir $outdir unless -e $outdir;


my $cov = $outdir . "Coverage_filtered/";
mkdir $cov unless -e $cov;
my $het = $outdir . "HetSites_filtered/";
mkdir $het unless -e $het;
my $aln = $outdir . "Individual_ALNs/";
mkdir $aln unless -e $aln;
my $geno = $outdir . "Individual_GENOs/";
mkdir $geno unless -e $geno;
my $H  = $outdir . "IndividualH_filtered/";
mkdir $H unless -e $H;        
my $Non = $outdir . "Individual_Non_diallelic/";
mkdir $Non unless -e $Non;
my $snpid = $outdir . "Individual_SNPID/";
mkdir $snpid unless -e $snpid;
my $snp = $outdir . "Individual_SNPs/";
mkdir $snp unless -e $snp;
my $miss = $outdir . "missingData_filtered/";
mkdir $miss unless -e $miss;
my $sampleID = $outdir . "sampleID.txt";

system ("cp $fadir" . "sub*/fasta/alignment/Individual_ALNs/*" . " $aln");
system ("cp $fadir" .  "sub*/fasta/alignment/Individual_GENOs/*" . " $geno");
system ("cp $fadir" .  "sub*/fasta/alignment/Individual_Non_diallelic/*" . " $Non");
system ("cp $fadir" .  "sub*/fasta/alignment/Individual_SNPID/*" . " $snpid");
system ("cp $fadir" .  "sub*/fasta/alignment/Individual_SNPs/*" .  " $snp");
system ("cp $fadir" .  "sub*/fasta/alignment/missingData_filtered/*" . " $miss");
system ("cp $fadir" .  "sub*/fasta/alignment/HetSites_filtered/*" . " $het");
system ("cp $fadir" .  "sub1/fasta/alignment/sampleID.txt $outdir");
system ("cp $fadir" .  "sub1/fasta/alignment/sampleID.txt $aln");
system ("cp $fadir" .  "sub1/fasta/alignment/sampleID.txt $geno");
system ("cp $fadir" .  "sub1/fasta/alignment/sampleID.txt $H");
system ("cp $fadir" .  "sub1/fasta/alignment/sampleID.txt $Non");
system ("cp $fadir" .  "sub1/fasta/alignment/sampleID.txt $snp");


my @cov = <$fadir/sub1/fasta/alignment/Coverage_filtered/*>;
foreach (@cov) {
  my $file = basename ($_);
  system ("cat $fadir" . "sub*/fasta/alignment/Coverage_filtered/*$file > $cov$file");
}

my @H = <$fadir*sub1/fasta/alignment/IndividualH_filtered/*>;
foreach (@H) {
  my $file = basename ($_);
  system ("cat $fadir" . "sub*/fasta/alignment/IndividualH_filtered/*$file > $H$file");
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
