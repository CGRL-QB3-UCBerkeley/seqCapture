#!/usr/bin/env perl
use warnings;
use strict;

die(qq/

Usage: seqCapture prealign  [options]

options:

-f     DIR      folder with all filtered individual fasta files from target regions
                (AAA.filtered.fasta, BBB.filtered.fasta, CCC.filtered.fasta...)
-s     FILE     target name file (XXXX_rename or XXXX_rename_compared.txt)
                generate by seqCap intarget;
-o     DIR      result folder
-n     INT      how many subsets do you want to have [10]
    
\n\n/) unless (@ARGV);
my %opts = (f=>undef, o=>undef, n=>10, s=>undef);
getopts('f:o:n:s:', \%opts);

my $fadir = redir ($opts{f});
my $outdir = redir ($opts{o});
mkdir $outdir unless -e $outdir;

my $name = $opts{s};
my $n = $opts{n};

my @names;
open (my $in, "<",  $name);
while (<$in>){
  chomp (my $a = $_);
  if ($a =~ /^>(Contig\d+)/) {
    push @names, $1;
  } 
}
close $in;

my $ii = 0;
my $d = 1;
foreach ( part {$ii++ % $n} @names ) {  
  my @array = @{$_};
  my $resdir = $outdir . "sub" . $d . "/";
  $d++;
  mkdir $resdir unless -e $resdir;
  my $subdir = $resdir . "fasta/";
  mkdir $subdir unless -e $subdir;
  
  my @files = <$fadir*fasta>;
  foreach (@files) {     
    my $file = $_;
    next if $file =~ /\S+_filtered\.fasta\.2/;
    my $name = $1 if basename ($file) =~ /(\S+)_filtered\.fasta$/;
    my $newout = $subdir . $name . "_filtered.fasta";
    open (my $out, ">", $newout);
    
    open (my $in, "<", $file);
    my %subseq;
    while (<$in>) {
      chomp (my $a = $_);
      if ($a =~ /^>\S+_(Contig\d+)/) {
	my $id  = $1;
	chomp (my $seq = <$in>);
	$subseq{$id} = $seq;
      }
    }
    close $in;
    
    foreach my $contig (sort {$a cmp $b} keys %subseq) {
      if (grep {$_ eq $contig} @array) {
	print $out ">", $name, "_", $contig, "\n";
	print $out $subseq{$contig}, "\n";
      }
    }
    close $out;
  }    
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
