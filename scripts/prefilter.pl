#!/usr/bin/env perl
use warnings;
use strict;
use Tie::Array::Packed;
use Getopt::Std;

die (qq/

Usage: seqCapture prefilter  [options]

-a    INT     a folder with all bams (make sure to only 
              use bams to be ananlyzed 
              in the downstream analysis)
-r    DIR     result dir 
-m    INT     Minimal percentile threthrold 
              to keep a contig [1]            
-n    INT     Maximal percentile threthrold 
              to keep a contig [99] 
-s    INT     S standard deviations of the mean 
              to keep a contig [3]
    
\n/) if !@ARGV;

my %opts = (a=>undef, r=>undef, m=>1, n=>99,s => 3);  
getopts('a:r:m:n:s:', \%opts);

my $resdir = redir ($opts{r});
mkdir $resdir unless -e $resdir;
my $bamdir = redir ($opts{a});
my $min = $opts{m};
my $max = $opts{n}; 
my $s = $opts{s};


my $file;
if (! -f $resdir . 'depth.txt') {
  $file = $resdir . 'depth.txt';
  my @depthall = `samtools depth  $bamdir*bam`;
  open (OUT, ">",  $file);
  foreach (@depthall) {
    chomp (my @a = split /\s+/, $_);
    my $chr = $a[0];
    my $pos = $a[1];
    my @s = @a[2..$#a];
    my $add = 0;
    my $count = 0;
    foreach (@s) {
      my $n = $_;
      $add = $add + $n;
      $count++ if $n > 0; 
    }
    my $avg = 0;
    $avg = sprintf("%.3f", $add/$count) if $count > 0;
    
    print OUT $chr, "\t", $pos, "\t", $avg, "\n"; 
  }
  close OUT;
}
else {
  $file = $resdir . 'depth.txt';
}


my $per_gene_depth = $resdir . "summary_gene_depth.txt";
my $gene_depth_percentile = $resdir . "summary_gene_depth_percentile.txt";
my $site_depth_percentile = $resdir . "summary_site_depth_percentile.txt";
my $gene_outside_percentile = $resdir . "summary_gene_outside_percentile.txt";
my $gene_outside_sd = $resdir . "summary_gene_outside_sd_filter.txt";
my $site_depth = $resdir . "summary_site_depth.txt";
my $site_gene_filtered = $resdir . "ifFiltered.txt";

tie my @site, 'Tie::Array::Packed::Number';
tie my @gene, 'Tie::Array::Packed::Number';

my $stotalq = 0;
my $scount = 0;
my $sexpxsqr = 0;

open (IN, "<", $file);
open (S, ">", $site_depth);

my %gene;
while (<IN>) {
  chomp (my @l =split /\s+/, $_);
  $gene{$l[0]}{'dep'} += $l[2];
  print S $l[2], "\n";
  $gene{$l[0]}{'count'} ++;
  push @site, $l[2]; 
  $stotalq = $stotalq + $l[2];
  $sexpxsqr = $sexpxsqr + ($l[2] ** 2);
  $scount++;
}
close S;


my $sexpq = $stotalq/$scount; 
my $svarq = ($sexpxsqr - (($stotalq ** 2) / $scount)) / ($scount - 1); 
my $sstd = sqrt ($svarq);


my $totalq = 0;
my $count = 0;
my $expxsqr = 0;

open (OUT1, ">", $per_gene_depth);

foreach my $g (sort { $a cmp $b} keys %gene) {
  print OUT1 $g, "\t",  $gene{$g}{'dep'}/$gene{$g}{'count'},"\n";
  $totalq = $totalq + $gene{$g}{'dep'}/$gene{$g}{'count'};
  $expxsqr = $expxsqr + (($gene{$g}{'dep'}/$gene{$g}{'count'}) ** 2);
  $count++;
  push @gene, $gene{$g}{'dep'}/$gene{$g}{'count'} ;
}
close OUT1;

my $expq = $totalq/$count; 
my $varq = ($expxsqr - (($totalq ** 2) / $count)) / ($count - 1); 
my $std = sqrt ($varq);

my $high;
my $low;

open (OUT2, ">", $gene_depth_percentile);
tied(@gene)->sort;
print OUT2 "gene depth mean: ", sprintf("%.2f", $expq), "\n";
print OUT2 "gene depth standard deviation: ", sprintf("%.2f", $std), "\n";
print "gene depth mean: ", sprintf("%.2f", $expq), "\n";
print "gene depth standard deviation: ", sprintf("%.2f", $std), "\n";

print OUT2 "gene percentile:", "\n";
print "gene percentile:", "\n";


open (IF, ">", $site_gene_filtered);
print IF "coverage\tpercentile\tnumberOfSitesRemoved\tpercentSitesRemoved\n"; 

foreach my $id (0, 0.5, 1, 1.5, 2, 2.5, 5, 10,20,30,40,50,60,70,80,90,95, 97.5,98,98.5,99,99.5,100) {
  seek IN , 0, 0;
  print OUT2 $id . " percentile: ",  sprintf("%.5f", $gene[$#gene*$id/100]), "\n";
  $low = sprintf("%.5f", $gene[$#gene*$id/100]) if $id == $min; 
  $high = sprintf("%.5f", $gene[$#gene*$id/100]) if $id == $max; 
  #$scount
  my $totallost = 0;
  if ($id >= 80) {
   
    my %contighigh;
    while (<IN>) {
      chomp (my @a = split /\s+/,$_);   
      my $cov = $a[2];
      $contighigh{$a[0]}{$a[1]}++ if $cov >  $gene[$#gene*$id/100]; 
      $totallost ++ if $cov >  $gene[$#gene*$id/100];    
    }
    my $contig_high_to_remove = $resdir . "sites_to_remove_higherThan_" . $id . "th_percentile". ".txt";
    open (HIGH, ">", $contig_high_to_remove);
    foreach my $contig (sort {$a cmp $b} keys %contighigh) {
      foreach my $pos (sort {$a <=> $b} keys  %{$contighigh{$contig}} ) {	
	print HIGH $contig,"\t", $pos, "\n";
      }
    }
    close HIGH;
    undef %contighigh;
    
    print IF sprintf("%.5f", $gene[$#gene*$id/100]), ">=", $id, "\t", $totallost, "\t", sprintf("%.5f",$totallost/$scount ), "\n";
  }
  if ($id <= 20) {
    my %contiglow;
    while (<IN>) {
      chomp (my @a = split /\s+/,$_);
      
      my $cov = $a[2];
      $contiglow{$a[0]}{$a[1]}++  if $cov <   $gene[$#gene*$id/100];
      $totallost ++ if $cov <   $gene[$#gene*$id/100];    
    }

    my $contig_low_to_remove = $resdir . "sites_to_remove_lowerThan_" . $id . "th_percentile". ".txt";
    open (LOW, ">", $contig_low_to_remove);
    foreach my $contig (sort {$a cmp $b} keys %contiglow) {
      foreach my $pos (sort {$a <=> $b} keys  %{$contiglow{$contig}} ) {
	print LOW $contig,"\t", $pos, "\n";
      }
    }
    close LOW;
    undef %contiglow;
    
    print IF sprintf("%.5f", $gene[$#gene*$id/100]), "<=", $id, "\t", $totallost, "\t", sprintf("%.5f",$totallost/$scount ), "\n";
  }
  print $id . " percentile: ",  sprintf("%.5f", $gene[$#gene*$id/100]), "\n";
    
}
close IN;
close IF;
close OUT2;
print "~~~~~~~~~~~~~~~~~~~", "\n";
open (OUT3, ">", $site_depth_percentile);

print OUT3 "site depth mean: ", sprintf("%.5f", $sexpq), "\n";
print OUT3 "site depth standard deviation: ", sprintf("%.5f", $sstd), "\n";
print "site depth mean: ", sprintf("%.5f", $sexpq), "\n";
print "site depth standard deviation: ", sprintf("%.5f", $sstd), "\n";


print OUT3 "site_percentile:", "\n"; 
print "site_percentile:", "\n";


tied(@site)->sort;
for (0, 0.5, 1, 1.5, 2, 2.5, 5, 10,20,30,40,50,60,70,80,90,95, 97.5,98,98.5,99,99.5,100) {
  
  print OUT3 $_ . " percentile: ",  sprintf("%.5f", $site[$#site*$_/100]), "\n";
  print $_ . " percentile: ",  sprintf("%.5f", $site[$#site*$_/100]), "\n";
  
}
close OUT3;
open (GENE, "<", $per_gene_depth);
open (OUT4, ">", $gene_outside_percentile);
open (OUT5, ">", $gene_outside_sd);


while (<GENE>) {
  chomp (my @a = split /\s+/,$_);
  print OUT4 $a[0], "\n" if ($a[1] < $low || $a[1] > $high);
  print OUT5 $a[0], "\n" if ($a[1] < $expq -  $std * $s || $a[1] >  $expq + $std * $s);
}

close OUT4;
close OUT5;
close GENE;


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
