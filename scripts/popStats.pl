#!/usr/bin/env perl

###########################################################################
## Version 0.1.0 Jul 31 2016
## Author:Ke Bi (kebi@berkeley.edu)
## Bioinformatics Scientist
## Computational Genomics Resource Laboratory
## California Institute for Quantitative Biosciences
## University of California, Berkeley
## 238 Koshland Hall
## Berkeley, CA 94720-3102
## http://qb3.berkeley.edu/qb3/cgrl/


### The codes for MK test (function "mk") are from Alisha Holloway at UC Davis
### I made some modifications but please notice that this is Alisha Holloway's work!
###########################################################################

use warnings;
use strict;
use Bio::PopGen::Statistics;
use Bio::PopGen::Utilities;
use Bio::AlignIO;
use Bio::Tools::CodonTable;
use File::Basename;
use Getopt::Std;
use Getopt::Long;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );
use List::Util 'shuffle';
use List::Util qw< sum >;
use Statistics::Basic qw(:all);
use Math::Cephes qw(:all);
use Statistics::Distributions;
use Text::NSP::Measures::2D::Fisher2::twotailed;



&main;
exit;

sub main {
        &usage if (@ARGV<1);
        my $command = shift(@ARGV);
        my %TRANS = (stats=>\&stats, mk=>\&mk, Dstatsi=>\&Dstatsi, Dstatsp=>\&Dstatsp, Dxyp =>\&Dxyp,Dxyi =>\&Dxyi,fst => \&fst, dnds=> \&dnds);
        die("Unknown command \"$command\"\n") if (!defined($TRANS{$command}));
        &{$TRANS{$command}};
      }


sub usage {
  die(qq/

Usage: seqCapture popstats <command> [<arguments>]\n

Command: 

stats:         calculate Pi, Watterson theta, Tajima's D, 
               Fu and Li's D (if outgroup is provided)
mk:            perform McDonald-Kreitman tests for intraspecific
               samples with one or two outgroup samples
Dstatsi:       compute D statistics (ABBA-BABA test) for four
               samples at a time for each locus. It 
               also bootstrap the data to calculate 
               standard deviation, a z score and a p-value
Dstatsp:       compute D statistics (ABBA-BABA test) for population
               samples. 
Dxyp:          For each pair of populations, calculate global Dxy
               and Dxy for each marker 
Dxyi:          For each pair of individuals, calculate global Dxy
               and Dxy for each marker 
fst:           For each pair of populations, calculate pairwise Fst 
               for each site and globally (using Reynold's estimator)
dnds:          Infer positive selection based on branch-site model using
               codeml (paml package)

\n/);
}


sub fst {
   die (qq/
seqCapture popstats fst [options]

Basic Options:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-f    DIR     A folder with SNP file for each 
              population (example below)
              files must have extension ".SNP";
-u    DIR     outdir
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



example for SNP fileS in the folder specified by -f. 
        each line is a site. The first column is chromosome\/contigID, 
        the second column is position and each of the following column is 
        genotype call of each individual: "0" means homo for ref alleles, "1" means 
        hetero and "2" means homo for alt allele. "-1" means missing data. 
        Note: all SNP files must have exactly the same number of sites even though 
        they can contain different number of individuals. 
         

Pop1.SNP:
  
contig1  1       1       2       1       0       0       0       0       0       0             
contig1  2       0       0       0       0      -1       1       1       1       2
contig1  3       1       2       1       0       0       0       0       0       0             
contig1  4       0       0       0       0      -1       1       1       1       2 
contig2  1       1       2       1       0       0       0       0       0       0             
contig2  2      -1       0       0       0      -1       1       1       1       2 
contig2  3       1       2       1       0       0       0       0       0       0             
contig2  4       2       2       1       0      -1       1       1       1       2              
....

Pop2.SNP:
contig1  1       0       0       0       0       0       0       0                  
contig1  2       0       0       0       2      -1       0       0      
contig1  3       1       2       1       0       0       0       0                  
contig1  4       0       0       0       2      -1       1       1 
contig2  1       1       0       1       0       0       0       0                 
contig2  2      -1       0       0       2      -1       1       1 
contig2  3       1       0       1       0       0       0       0          
contig2  4       2       2       1       0      -1       1       1      

Pop3.SNP:
contig1  1       0       1       0       1       0                         
contig1  2       0       0       0       0      -1             
contig1  3       1       2       1       2       0                       
contig1  4       0       0       0       0      -1        
contig2  1       1       2       1       0       0                       
contig2  2      -1       0       0       0      -1        
contig2  3       1       2       1       0       0               
contig2  4       2       2       1       0      -1            
         
\n\n/) unless (@ARGV);
    
   my %opts = (f => undef, o=> undef);
   getopts('f:u:', \%opts);
   my $dir = redir ($opts{f});
   my $outdir = redir ($opts{u});
   mkdir $outdir unless -e $outdir;
   
   my @pop = <$dir*SNP>;
   my $out2 = $outdir . "pairwise_global.Fst";
   open (my $global, ">", $out2);
   
   for (my $i = 0; $i < scalar @pop ; $i++) {
     my $pop1 = $pop[$i];
     for (my $j= $i+1; $j < scalar @pop ; $j++) {
       my $pop2 = $pop[$j];
       my %freq;
      
       my $lib1 = $1 if basename ($pop1) =~ /(\S+)\.SNP/;
       my $lib2 = $1 if basename ($pop2) =~ /(\S+)\.SNP/;
       open (my $freq1, "<", $pop1);
       open (my $freq2, "<", $pop2);

       while (<$freq1>) {
	 chomp (my @l = split /\s+/, $_);
	 my $chr = $l[0];
	 my $pos = $l[1];
	 my @geno = @l[2 .. $#l];
	 my $alt = 0;
	 my $total = 0;
	 foreach my $g (@geno) {
	   $alt = $alt + 2 if $g == 2; 
	   $alt = $alt + 1 if $g == 1;
	   $total = $total + 2 unless $g == -1;
	 }
	 if ($total > 0 ) {
	   $freq{$chr}{$pos}{$lib1}{'alt'} = $alt;  
	   $freq{$chr}{$pos}{$lib1}{'total'} = $total;
	   $freq{$chr}{$pos}{'all'}{'alt'} += $alt;  
	   $freq{$chr}{$pos}{'all'}{'total'} += $total;
	 }
       }
       close $freq1;
       
       while (<$freq2>) {
	 chomp (my @l = split /\s+/, $_);
	 my $chr = $l[0];
	 my $pos = $l[1];
	 my @geno = @l[2 .. $#l];
	 my $alt = 0;
	 my $total = 0;
	 foreach my $g (@geno) {
	   $alt = $alt + 2 if $g == 2; 
	   $alt = $alt + 1 if $g == 1;
	   $total = $total + 2 unless $g == -1;
	 }
	 if ($total > 0 && $freq{$chr}{$pos}{$lib1}) {
	   $freq{$chr}{$pos}{$lib2}{'alt'} = $alt;  
	   $freq{$chr}{$pos}{$lib2}{'total'} = $total;
	   $freq{$chr}{$pos}{'all'}{'alt'} += $alt;  
	   $freq{$chr}{$pos}{'all'}{'total'} += $total;
	 }
	 if ($freq{$chr}{$pos}{$lib1} && $total == 0) {
	   delete $freq{$chr}{$pos};
	 }
       }
       close $freq2;

       
       my $out1 = $outdir . $lib1 . "_" . $lib2 . "_persite.Fst";
       open (my $persitefst, ">", $out1);       
       
       my $out3 = $outdir . $lib1 . "_" . $lib2 . "_contig.Fst";
       open (my $contigfst, ">", $out3);
       
       my $globala = 0;
       my $globalb = 0;
       
       foreach my $chr (sort {$a cmp $b} keys %freq) {
	 my $contiga = 0;
	 my $contigb = 0;    
	 foreach my $pos (sort {$a <=> $b} keys %{$freq{$chr}}) {
	   my $n1 = $freq{$chr}{$pos}{$lib1}{'total'}/2;
	   my $n2 = $freq{$chr}{$pos}{$lib2}{'total'}/2;
	   my $ps1 = 0;
	   my $ps2 = 0;
	   my $pst = 0;
	   $ps1 = $freq{$chr}{$pos}{$lib1}{'alt'}/$freq{$chr}{$pos}{$lib1}{'total'} if $freq{$chr}{$pos}{$lib1}{'total'} != 0;
	   $ps2 = $freq{$chr}{$pos}{$lib2}{'alt'}/$freq{$chr}{$pos}{$lib2}{'total'} if $freq{$chr}{$pos}{$lib2}{'total'} != 0;
	   $pst = $freq{$chr}{$pos}{'all'}{'alt'}/$freq{$chr}{$pos}{'all'}{'total'} if $freq{$chr}{$pos}{'all'}{'total'} != 0;	   
	   #print $chr, "\t", $pos, "\t", $n1, "\t", $n2, "\t", $ps1, "\t", $ps2, "\t", $pst, "\n";
	   my $a = 0;
	   my $b = 0;
	   $b = ($n1 * 2 * $ps1 * (1 - $ps1) + $n2 * 2 * $ps2 * (1 - $ps2)) / ($n1 + $n2 -1) if  ($n1 + $n2 -1) != 0;	    
	   $a = (4 * $n1 * (($ps1 - $pst) ** 2) + 4 * $n2 * (($ps2 -$pst) ** 2)  - $b)/(2 * ( (2 * $n1 * $n2)/($n1 + $n2) )) if $n1 != 0 && $n2 != 0;
	   
	   #$b =  4*0.071*(1-0.071) / 8 = 0.0355
	   #$a = 28x(0.071-0.0555)** + 8x(0.0555)*** / 6.22 = 0.006727 + 0.024688 -  0.1161/ 6.22 = -0.0008379

	   $a = 0 if $a < 0;
	   $b = 0 if $b < 0;
	   $globala += $a ;	   
	   $globalb += $a + $b ;	 
	   $contiga += $a;
	   $contigb += $a + $b;	   
	   my $sitefst = 0;
	   $sitefst = $a / ($a + $b) if ($a + $b) != 0;
	   print $persitefst $chr, "\t", $pos, "\t", sprintf("%.4f",$a), "\t", sprintf("%.4f",$b), "\t", sprintf("%.4f",$sitefst), "\n";
	 }
	 
	 print $contigfst $chr, "\t", sprintf("%.4f", $contiga/$contigb),  "\n" if $contigb !=0;
	 print $contigfst $chr, "\t", "0",  "\n" if $contigb == 0;
       }
       close $persitefst;
       close $contigfst;
       
       print $global $lib1, "_vs_", $lib2, "\t", sprintf("%.4f",$globala/$globalb), "\n" if $globalb !=0;
       print $global $lib1, "_vs_", $lib2, "\t", "0", "\n" if $globalb == 0;
     }
   }
   close $global;
 }


sub dnds {
  die (qq/
seqCapture popstats dnds [options]

Basic Options:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-f    DIR     A folder with alignment files in phy format
              AND the corresponding tree files in newick format, 
              with the branch of interest specified by "#1".
              The sequences in phy file must be aligned by codons
              If the name for phy file is: 
              ABC.phy
              Then the name for the corresponding newick file 
              muse be:
              ABC.tree 
-u    FILE    Outfile dir
-p    FLOAT   significant level for LRT? [0.05]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         
\n\n/) unless (@ARGV);
    
    my %opts = (f => undef,  u => undef, p => 0.05);
    getopts('f:u:p:', \%opts);
    my $dir = redir ($opts{f});
    my $outdir = redir ($opts{u});
    mkdir $outdir unless -e $outdir;
    my $p = $opts{p};

    
    my @phy = <$dir*.phy>;
    
    foreach (@phy) {
      my $phy = $_;
      my $lib = $1 if basename ($phy) =~/(\S+)\.phy/;
      my $newick = $dir . $lib . ".tree";
      my $in1 =  $outdir . $lib . ".fixed.ctl";
      my $in2 =  $outdir . $lib . ".ctl";
      my $out1 = $outdir . $lib . ".fixed.mlc";
      my $out2 = $outdir . $lib . ".mlc";

      open (OUT1, ">", $in1);
      open (OUT2, ">", $in2);

      print OUT1 'seqfile  = ', $phy,"\n";
      print OUT1 'treefile = ', $newick, "\n";
      print OUT1 'outfile  = ', $out1,"\n";
      print OUT1 'noisy = 9      * 0,1,2,3,9: how much rubbish on the screen',"\n";     
      print OUT1 'verbose = 1    * 1: detailed output, 0: concise output', "\n";    
      print OUT1 'runmode = 0    * 0: user tree;  1: semi-automatic;  2: automatic; 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise', "\n";    
      print OUT1 'seqtype = 1    * 1:codons; 2:AAs; 3:codons-->AAs', "\n";  
      print OUT1 'CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table', "\n";  
      print OUT1 'clock = 0       * 0: no clock, unrooted tree, 1: clock, rooted tree', "\n";  
      print OUT1 'aaDist = 0      * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}', "\n";  
      print OUT1 'model = 2       * models for codons: 0:one, 1:b, 2:2 or more dN/dS ratios for branches', "\n";                  
      print OUT1 'NSsites = 2     * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete; 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal', "\n";
      
      print OUT1 'icode = 0       * 0:standard genetic code; 1:mammalian mt; 2-10:see below', "\n";
      print OUT1 'Mgene = 0       * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all',"\n";

      print OUT1 'fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated', "\n";
      print OUT1   'kappa = 2     * initial or fixed kappa', "\n";
      print OUT1 'fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate', "\n";
      print OUT1   'omega = 1     * initial or fixed omega, for codons or codon-based AAs', "\n";

      print OUT1  'getSE = 0      * 0: dont want them, 1: want S.E.s of estimates', "\n";
      print OUT1  'RateAncestor = 0       * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)', "\n";
      print OUT1  'Small_Diff = .45e-6    * Default value', "\n";
      print OUT1  'cleandata = 1          * remove sites with ambiguity data (1:yes, 0:no)?', "\n";
      print OUT1 'fix_blength = 0         * 0: ignore, -1: random, 1: initial, 2: fixed',"\n"; 
      
      close OUT1;

      
      print OUT2 'seqfile  = ', $phy,"\n";
      print OUT2 'treefile = ', $newick, "\n";
      print OUT2 'outfile  = ', $out2,"\n";
      print OUT2 'noisy = 9      * 0,1,2,3,9: how much rubbish on the screen',"\n";     
      print OUT2 'verbose = 1    * 1: detailed output, 0: concise output', "\n";    
      print OUT2 'runmode = 0    * 0: user tree;  1: semi-automatic;  2: automatic; 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise', "\n";    
      print OUT2 'seqtype = 1    * 1:codons; 2:AAs; 3:codons-->AAs', "\n";  
      print OUT2 'CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table', "\n";  
      print OUT2 'clock = 0       * 0: no clock, unrooted tree, 1: clock, rooted tree', "\n";  
      print OUT2 'aaDist = 0      * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}', "\n";  
      print OUT2 'model = 2       * models for codons: 0:one, 1:b, 2:2 or more dN/dS ratios for branches', "\n";                  
      print OUT2 'NSsites = 2     * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete; 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal', "\n";
      
      print OUT2 'icode = 0       * 0:standard genetic code; 1:mammalian mt; 2-10:see below', "\n";
      print OUT2 'Mgene = 0       * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all',"\n";

      print OUT2 'fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated', "\n";
      print OUT2 'kappa = 2       * initial or fixed kappa', "\n";
      print OUT2 'fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate', "\n";
      print OUT2 'omega = 1       * initial or fixed omega, for codons or codon-based AAs', "\n";

      print OUT2  'getSE = 0      * 0: dont want them, 1: want S.E.s of estimates', "\n";
      print OUT2  'RateAncestor = 0       * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)', "\n";
      print OUT2  'Small_Diff = .45e-6    * Default value', "\n";
      print OUT2  'cleandata = 1          * remove sites with ambiguity data (1:yes, 0:no)?', "\n";
      print OUT2 'fix_blength = 0         * 0: ignore, -1: random, 1: initial, 2: fixed',"\n"; 

      close OUT2;

      system ("codeml $in1");
      system ("codeml $in2");

      my $ln1 = 0;
      my $ln2 = 0;
      my $sam1 = 0;
      my $sam2 = 0;
      
      open (my $null, "<", $out1);
      while (<$null>) {
	chomp (my $line = $_);
	if ($line =~ /lnL\(ntime/) {
	  #lnL(ntime: 41  np: 45):  -4710.222252      +0.000000
	  my @a = split /\s+/, $line;
	  my @b = split /\)/, $a[3];
	  $ln1 = $a[4];
	  $sam1 = $b[0];
	}
      }
      close $null;

      open (my $model, "<", $out2);
      while (<$model>) {
	chomp (my $line = $_);
	if ($line =~ /lnL\(ntime/) {
	  #lnL(ntime: 41  np: 45):  -4710.222252      +0.000000
	  my @a = split /\s+/, $line;
	  my @b = split /\)/, $a[3];
	  $ln2 = $a[4];
	  $sam2 = $b[0];
	}
      }
      close $model;
      
      my $outfile = $outdir . $lib . ".dnds.out";
      open (my $of, ">", $outfile);
      my $LRT = 2 * ($ln2 - $ln1); 
      my $df = $sam2 - $sam1;
      my $chisprob = 1;
      $chisprob = Statistics::Distributions::chisqrprob ($df, $LRT);
      print $of "p-value: ", $chisprob, "\n";
      
      if ($chisprob <= $p) {
	open (my $in, "<", $out2);
	while (<$in>) {
	  chomp (my $line = $_);
	  if ($line =~ /Bayes Empirical Bayes/) {	   
	    do { 	      
	      print $of $line,"\n";
	      chomp ($line = <$in>);
	      
	    } until ($line =~ /The grid/);
	    
	  }
	  
	}	
	close $in;    
      }
      else {
	print $of "not significant\n";	
      }      
      close $of;     
    }
  } 


sub Dxyi {
    die (qq/
seqCapture popstats Dxyi [options]

Basic Options:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-f    FILE     A SNP file (example format below);
-u    FILE     Outfile prefix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



example for a SNP file specified by -f. 
        each line is a site. The first column is chromosome\/contigID, 
        the second column is position and each of the following column is 
        genotype call of each individual: "0" means homo for ref alleles, "1" means 
        hetero and "2" means homo for alt allele. "-1" means missing data. 
         

Pop1.SNP:
  
contig1  1       1       2       1       0       0       0       0       0       0             
contig1  2       0       0       0       0      -1       1       1       1       2
contig1  3       1       2       1       0       0       0       0       0       0             
contig1  4       0       0       0       0      -1       1       1       1       2 
contig2  1       1       2       1       0       0       0       0       0       0             
contig2  2      -1       0       0       0      -1       1       1       1       2 
contig2  3       1       2       1       0       0       0       0       0       0             
contig2  4       2       2       1       0      -1       1       1       1       2              

         
\n\n/) unless (@ARGV);
    
    my %opts = (f=>undef,  u=>undef);
    getopts('f:u:', \%opts);
    my $infile = $opts{f};
    my $dxyp_global = $opts{u}. "_individual_dxy.txt";

    open (my $in, "<", $infile);
    open (my $out, ">", $dxyp_global);
    
    chomp (my @first = split /\s+/, <$in>);
    seek $in, 0, 0;
    
    print $out "markerID\t";
    my $numsamples = scalar @first - 2;
    for (my $i = 1; $i <=$numsamples; $i ++ ) {
      for (my $j = $i + 1; $j <=$numsamples; $j ++) {
	print $out "ind", $i, "_ind", $j,"\t"; 
      }
    }
    print $out "\n";
    print "Now computing Dxy between each pair of individuals...\n";
    print "The global estimate of Dxy for each pair of individuals is located at the last line of the output file.\n";

    my %p1;
    my %dxy;
    my %overall;
    
    while (<$in>) {
      chomp (my @line = split /\s+/, $_);
      my $contig = $line[0];
      my $pos = $line[1];
      push @{$p1{$contig}{$pos}},  @line[2..$#line];
    }
    close $in;

    foreach my $id (sort {$a cmp $b} keys %p1) {
      foreach my $pos (sort {$a cmp $b} keys %{$p1{$id}}) {	
	for (my $i = 0; $i < scalar @{$p1{$id}{$pos}}; $i++) {
	  my $geno1 = $p1{$id}{$pos}[$i];
	  next if $geno1 == -1;
	  for (my $j = $i+1; $j < scalar @{$p1{$id}{$pos}}; $j++) {
	    my $geno2 = $p1{$id}{$pos}[$j];
	    next if $geno2 == -1;
	    	    
	    if ($geno1 == $geno2) {
	      if ($geno1 == 1) {
		$dxy{$id}{$i}{$j}{'data'} += 1;
		$overall{$i}{$j}{'data'} += 1
	      }
	    } 
	    if ($geno1 != $geno2) {
	      if (abs($geno1-$geno2) == 2) {
		$dxy{$id}{$i}{$j}{'data'} += 2;
		$overall{$i}{$j}{'data'} += 2;
	      }
	      if (abs($geno1-$geno2) == 1) {
		$dxy{$id}{$i}{$j}{'data'} += 1;
		$overall{$i}{$j}{'data'} += 1;
	      }
	    }
	    $dxy{$id}{$i}{$j}{'site'} ++; ;
	    $overall{$i}{$j}{'site'} ++;
 
	  }
	}
	
      }
    }

    
    foreach my $id (sort {$a cmp $b} keys %dxy) {
       print $out $id, "\t";
      foreach my $i (sort {$a <=> $b} keys %{$dxy{$id}}) {
	foreach my $j (sort {$a <=> $b} keys %{$dxy{$id}{$i}}) {
	  $dxy{$id}{$i}{$j}{'data'} = 0.000 if !$dxy{$id}{$i}{$j}{'data'};
	  print $out sprintf("%.4f",$dxy{$id}{$i}{$j}{'data'}/$dxy{$id}{$i}{$j}{'site'}), "\t" if  ($dxy{$id}{$i}{$j}{'site'});
	  print $out "na", "\t" if (!$dxy{$id}{$i}{$j}{'site'});
	}
      }
      print $out "\n";
     }
    
    print $out "overall\t";
    foreach my $i (sort {$a <=> $b} keys %overall) {
	foreach my $j (sort {$a <=> $b} keys %{$overall{$i}}) {	  
	  print $out sprintf("%.4f",$overall{$i}{$j}{'data'}/$overall{$i}{$j}{'site'}), "\t" if $overall{$i}{$j}{'site'} > 0;
	  print $out "na", "\t" if $overall{$i}{$j}{'site'} <= 0;
	}
      }
    print $out "\n";

    
    close $out;
  }




sub Dxyp {
    die (qq/
seqCapture popstats Dxyp [options]

Basic Options:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-f    DIR     A folder with SNP file for each 
              population (example below)
              files must have extension ".SNP";
-u    FILE    Outfile prefix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



example for SNP fileS in the folder specified by -f. 
        each line is a site. The first column is chromosome\/contigID, 
        the second column is position and each of the following column is 
        genotype call of each individual: "0" means homo for ref alleles, "1" means 
        hetero and "2" means homo for alt allele. "-1" means missing data. 
        Note: all SNP files must have exactly the same number of sites even though 
        they can contain different number of individuals. 
         

Pop1.SNP:
  
contig1  1       1       2       1       0       0       0       0       0       0             
contig1  2       0       0       0       0      -1       1       1       1       2
contig1  3       1       2       1       0       0       0       0       0       0             
contig1  4       0       0       0       0      -1       1       1       1       2 
contig2  1       1       2       1       0       0       0       0       0       0             
contig2  2      -1       0       0       0      -1       1       1       1       2 
contig2  3       1       2       1       0       0       0       0       0       0             
contig2  4       2       2       1       0      -1       1       1       1       2              
....

Pop2.SNP:
contig1  1       0       0       0       0       0       0       0                  
contig1  2       0       0       0       0      -1       0       0      
contig1  3       1       2       1       0       0       0       0                  
contig1  4       0       0       0       0      -1       1       1 
contig2  1       1       2       1       0       0       0       0                 
contig2  2      -1       0       0       0      -1       1       1 
contig2  3       1       2       1       0       0       0       0          
contig2  4       2       2       1       0      -1       1       1             
....

         
\n\n/) unless (@ARGV);
    
    my %opts = (f=>undef,  u=>undef);
    getopts('f:u:', \%opts);
    my $dir = redir ($opts{f});
    my @snpsfiles = <$dir*SNP>;

    my $dxyp_global = $opts{u}. ".global.txt";
    open (my $out, ">", $dxyp_global);
    
    print $out "pop1\tpop2\tDxy\n";

    my %pops;
    for (my $i = 0; $i < scalar (@snpsfiles); $i ++) {
      my $file1 = $snpsfiles[$i];
      my $pop1 = $1 if basename ($file1) =~ /(\S+)\.SNP/;
      for (my $j = $i+1; $j < scalar (@snpsfiles); $j++) {
	my $file2 = $snpsfiles[$j];
	my $pop2 = $1 if basename ($file2) =~ /(\S+)\.SNP/;
	print "Now computing Dxy between $pop1 and $pop2...\n";
	my $dxyp_pair = $opts{u}. "_" . $pop1 . "_". $pop2 . "_marker_by_marker.txt";
	my $global = calDxy ($file1, $file2, $dxyp_pair);

	print $out $pop1, "\t", $pop2, "\t", $global, "\n";
	
      }
    }
    close $out;
  }

sub calDxy {
  my ($pop1, $pop2, $out) = @_;
  open (G1, "<", $pop1);
  open (G2, "<", $pop2);
  my %p1; my %p2;
  while (<G1>) {
    chomp (my @line = split /\s+/, $_);
    my $contig = $line[0];
    my $pos = $line[1];
    push @{$p1{$contig}{$pos}},  @line[2..$#line];
    chomp (my $l = <G2>);
    my @a = split /\s+/, $l;
    push @{$p2{$a[0]}{$a[1]}},  @a[2..$#a];

    if ( max (@line[2..$#line]) == -1 ||  max (@a[2..$#a]) == -1 ) {
      delete $p1{$contig}{$pos};
      delete $p2{$contig}{$pos};
    }
    
  }
  close G1; close G2;

  my %Dxy;
  my $all = 0;
  my $sites = 0;
  my $global = "nan"; 
  foreach my $id (sort {$a cmp $b} keys %p1) {
    foreach my $pos (sort {$a cmp $b} keys %{$p1{$id}}) {
      print "processing ",$id, "\t", $pos, "! \n"; 
      my $n = 0;
      my $count = 0;
      foreach  (@{$p1{$id}{$pos}}) {
	my $geno1 = $_; 	
	next if $geno1 == -1;
	foreach (@{$p2{$id}{$pos}}) {
	  my $geno2 = $_;
	  next if $geno2 == -1;
	  if ($geno1 == $geno2) {
	    $count++;
	    if ($geno1 != 1) {
	      $n += 0;
	    }
	    if ($geno1 == 1) {
	      $n += 1;
	    }
	  } 
	  if ($geno1 != $geno2) {
	    $count++;
	    if (abs($geno1-$geno2) == 2) {
	      $n += 2;
	    }
	    if (abs($geno1-$geno2) == 1) {
	      $n += 1;
	    }
	  }
	} 
      }  
      $all += $n/$count if ($count > 0 ) ;
      $sites ++  if ($count > 0 ) ;
      $Dxy{$id}{'data'} += $n/$count if ($count > 0 ) ;
      $Dxy{$id}{'site'} ++ if ($count > 0 ) ;
    }    
  }
  $global = sprintf("%.4f", $all/$sites);
  open (OUT, ">", $out );
  print OUT "markerID", "\t","Dxy", "\n";
  foreach my $id (sort {$a cmp $b} keys %Dxy) {
    print OUT $id, "\t", sprintf("%.4f",$Dxy{$id}{'data'}/$Dxy{$id}{'site'}), "\n";
  }
  close OUT;

  return ($global);
  
  
  
}


sub Dstatsp {
    die (qq/
seqCapture popstats Dstatsp [options]

Basic Options:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-f    FIR     A folder with fasta files with alignments of population 
              samples from four species.(example below)
-b    INT     Number of replicates for resampling
-u    FILE    Outfile name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note: Sites that contain missing data or ambiguity will be dropped.
      Sites that are segregating in outgroups will be dropped.
    

example for a fasta file in the folder specified by -f. 
        Sequence must not have line breaks: one line for id, 
        and one line for sequence for each sample. 
        must be named as pop1, pop2, pop3, outgroup!

   >pop1
   ATGGCTAACGTAATTAAACCCGTTTTGACTTACCAGTTAGATGGCTCCAATCGTGATTTT
   >pop1
   TTGGCTAACGTAATTAAAACCGTTTTGACTTACCAGTTAGATGGCTCCAATCGTGATTTT
   >pop2
   ATGGCTAAGGTAATTAAAACCGTTTTTACTTACCAGTTAGATGGCTCCAATCGTGATTTT
   >pop2
   TTGGCTATGGTAATTAAAACCGTTTTTACTTACCAGTTAGATGGCTCCAATCGTGATTTT
   >pop2
   ATGGCTAACGTAATTAAAACCGTTTTTACTTACCAGTTAGATGGCTCCAATCGTGATTTT
   >pop3
   ATGGCTATGGTAATTAAAACCGTTTTTACTTACCAGTTAGATGGCTCCAATCGTGATTTT
   >pop3
   ATGGCTAAGGTAATTAAAACCGTTTTTACTTACCAGTTAGATGGCTCCAATCGTGATTTT
   >outgroup
   TTGGCTAACGTAATTAACACCGTTTTGACTTACCAGTTAGATGGCTCCAATCGTGATTTT
   >outgroup
   TTGGCTAACGTAATTAACACCGTTTTGACTTACCAGTTAGATGGCTCCAATCGTGATTTT
 
   The tree structure for the four samples should look like below
          
          #### pop1 
          #      
       #### 
       #  #
    ####  #### pop2  
    #  #
    #  #### pop3   
    # 
    #### outgroup
 
         
\n\n/) unless (@ARGV);
    
    my %opts = (f=>undef, b=>undef, u=>undef);
    getopts('f:u:b:', \%opts);
    
    my $dir = redir ($opts{f});
    my @aln = <$dir*fa>;
    my $boot;
    $boot = $opts{b} if $opts{b};
    
    my $out = $opts{u};
    open (OUT, ">", $out);
    print OUT "Locus\tTotalNumSites\tSitesVariable\tNumABBA-BABASegregting\traw_D\tD_sd\tZ-score\tp-value\n";
    
    
    foreach (@aln) {
      my $file = $_;
      my $lib = $1 if basename ($file) =~ /(\S+)\.fa/;
      open (my $filein , "<", $file);
      my %pop;
      while (<$filein>) {
	chomp (my $line = $_);
	if ($line =~ /^>(\S+)/) {
	  chomp (my $seq = <$filein>);
	  push @{$pop{$1}}, $seq;	  
	}
      }
      close $filein;
      my ($d, $total, $seg) = calDp (\%pop);
      
      if (!$boot) {
      print OUT $lib, "\t", length ($pop{'pop1'}[0]), "\t", $total, "\t", $seg, "\t", sprintf("%.4f", $d),  "\t", "na", "\t", "na", "\t", "na", "\n"; 
      }
      if ($boot) {
	print "Now doing bootstrapping!\n";
	my @bootp;

	my @a = (1 .. length ($pop{'pop1'}[0]));
	my $dd = 0;
	for (1 .. $boot){
	  $dd ++;
	  print "round $dd...\n" if ($dd/50 == round $dd/50);
	  my %popboot;
	  for (my $m = 0 ; $m < length ($pop{'pop1'}[0]); $m++ ) {	    	     
	    my $index = $a[rand @a] - 1;
	    foreach my $p (sort {$a cmp $b} keys %pop ) {	      
	      my $d = 0;
	      foreach (@{$pop{$p}}) {
		my @seq = split //, $_;
		$popboot{$p}[$d] .= $seq[$index];
		$d++;
	      }
	    }
	  }  
	  my ($rawd, $totalS, $segS) = calDp (\%popboot); 
	  push @bootp, $rawd;
	}
	my $sdev = sprintf("%.4f", stddev (@bootp));
	my $mean =  sprintf("%.4f", mean (@bootp));
	my $zscore = 'nan';
	$zscore = sprintf("%.4f", $d /$sdev) unless $sdev == 0;
	my $p_value = "nan";
	$p_value = (1-ndtr(abs($zscore)))*2 unless $zscore eq 'nan'; ##two tailed
	
	print OUT $lib, "\t",  length ($pop{'pop1'}[0]), "\t", $total, "\t", $seg, "\t", sprintf("%.4f", $d),  "\t", $sdev, "\t", $zscore, "\t", sprintf("%.5f",$p_value), "\n"; 
      } #if $boot;
    }#foreach (@aln) {
    close OUT;
  }

sub calDp {
  my ($pops) = @_;
  my %pop = %{$pops};
  my %m = ("N" => 0,"?" => 0, "-" => 0, "Y"=>0, "S"=>0, "K"=>0, "M"=>0, "R"=>0, "W"=>0, "B"=>0, "D"=>0, "V"=>0, "H"=>0,);
  my %am =  ("Y"=>0, "S"=>0, "K"=>0, "M"=>0, "R"=>0, "W"=>0);
  
  my @a = (0 .. (length $pop{'pop1'}[0]) -1);
  
  my $up = 0;
  my $down = 0;
  my $seg = 0;
  my $total = 0;
  foreach my $in (@a) {
    my @all;
    my @pop1;
    my @pop2;
    my @pop3;
    my @out;
    my $drop = 0;
    foreach my $p (sort {$a cmp $b} keys %pop ) {
      foreach (@{$pop{$p}}) {
	my @seq = split //, $_;
	push @all, $seq[$in];
	$drop++ if $m{$seq[$in]};
	push @pop1, $seq[$in] if $p eq 'pop1';
	push @pop2, $seq[$in] if $p eq 'pop2';
	push @pop3, $seq[$in] if $p eq 'pop3';
	push @out, $seq[$in] if $p eq 'outgroup';
      }
    }
	
    next unless (uniq ( @all) == 2);
    next if $drop > 0;
    next if (uniq (@out) > 1);
    
    $total++;
    my $ref = $out[0];
    
    my $refsize = scalar (@out);
    my $alt;
    for (@all) {
      $alt = $_ if $_ ne $ref;
    }
    
    my ($pop1alt, $pop1size) = countallele (\@pop1,$alt);
    my ($pop2alt, $pop2size) = countallele (\@pop2,$alt);
    my ($pop3alt, $pop3size) = countallele (\@pop3,$alt);
    
    
    unless (($pop1alt == $pop1size && $pop2alt == $pop2size) || ($pop1alt == 0 && $pop2alt == 0)) {
      unless ($pop3alt == 0) {
	# convert to referece allele
	$seg++;
	$up += (1-($pop1alt/$pop1size)) * ($pop2alt/$pop2size) * ($pop3alt/$pop3size)  - ($pop1alt/$pop1size) * (1-($pop2alt/$pop2size)) * ($pop3alt/$pop3size);
	$down += (1-($pop1alt/$pop1size)) * ($pop2alt/$pop2size) * ($pop3alt/$pop3size)  + ($pop1alt/$pop1size) * (1-($pop2alt/$pop2size)) * ($pop3alt/$pop3size);
	
      }
    }   
  }#foreach my $in (@a)
  
  my $d = 0;
  $d = $up/$down if $down;
  
  return ($d,$total, $seg);
}


sub countallele {
  my ($a, $alt) = @_;
  my @pop = @{$a};
  
  my $count = 0;
  my $size = 0;
  
  foreach (@pop) {
    $count++ if $_ eq $alt;
    $size++;
  }
  return ($count, $size);
}



sub Dstatsi {
    die (qq/
seqCapture popstats Dstatsi [options]

This script is inspired by the R package evobir

Basic Options:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-f    FIR     A folder with fasta files with alignments of 4 samples
              The order of the four samples matters! (example below)             
-u    FILE    Outfile name

Optional - resampling (recommended)

-b    INT     Bootstrapping: Number of replicates for bootstrapping or jackknife
-j    INT     Jackknife resampling: sizes to drop in jackknife,
              must be smaller than 1\/2 of the alignment length 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note: If -j is not specified then a bootstrapping method will be used. 
      Otherwise jackKnife will be used. 
      If both -j and -b are not specified, then no bootstrapping or jackknife 
      will be performed. 

Note: Sites that contain missing data or ambiguity will be dropped.

An example for a fasta file in the folder specified by -f. 
        Sequence must not have line breaks: one line for id, 
        and one line for sequence for each sample. 

   >sample1
   ATGGCTAACGTAATTAAACCCGTTTTGACTTACCAGTTAGATGGCTCCAATCGTGATTTT
   >sample2
   ATGGCTAACGTAATTAAAACCGTTTTGACTTACCAGTTAGATGGCTCCAATCGTGATTTT
   >sample3
   ATGGCTAACGTAATTAAAACCGTTTTTACTTACCAGTTAGATGGCTCCAATCGTGATTTT
   >sample4
   ATGGCTAACGTAATTAACACCGTTTTGACTTACCAGTTAGATGGCTCCAATCGTGATTTT
 
   The tree structure for the four samples should look like below
          
          #### sample1 
          #      
       #### 
       #  #
    ####  #### sample2  
    #  #
    #  #### sample3   
    # 
    #### sample4 
  
  (isn't this tree nice?) 
         
\n\n/) unless (@ARGV);
    
    my %opts = (f=>undef, b=>undef, u=>undef, j=>undef);
    getopts('f:b:u:j:', \%opts);
    
    my $dir = redir ($opts{f});
    my @aln = <$dir*fa>;

    my $boot;
    $boot = $opts{b} if $opts{b};
    my $jack;
    $jack = $opts{j} if $opts{j};
    
    my $out = $opts{u};
    open (OUT, ">", $out);
    print OUT "Locus\tTotalNumSites\tNumBBAA\tNumABBA\tNumBABA\traw_D\tD_sd\tZ-score\tp-value\n";
    
    foreach (@aln) {
     my $file = $_;
     my $lib = $1 if basename ($file) =~ /(\S+)\.fa/;
     open (my $filein , "<", $file);
     
     my $header = <$filein>;
     chomp (my @sp1 = split //, <$filein>);
     $header = <$filein>;
     chomp (my @sp2 = split //, <$filein>);
     $header = <$filein>;
     chomp (my @sp3 = split //, <$filein>);
     $header = <$filein>;
     chomp (my @sp4 = split //, <$filein>);
     close $filein;
       
     my ($rawd,$abba, $baba, $bbaa) = calD (\@sp1, \@sp2, \@sp3, \@sp4); 
     if (!$boot && !$jack){
       print "no bootstraping!\n";
       print OUT $lib, "\t", scalar @sp1 -1, "\t",$bbaa, "\t",  $abba, "\t", $baba, "\t", sprintf("%.4f", $rawd), "\t","na", "\t", "na", "\t", "na", "\n";
     }
     
     elsif ($boot && !$jack){
       print "Now doing boostrap!\n";
       my @bootd;       
       for (1 .. $boot){       	 
	 my @sps1;
	 my @sps2;
	 my @sps3;
	 my @sps4; 	 
	 my @a = (1 .. scalar @sp1 - 1);	 
	 for (my $m = 0 ; $m < scalar @sp1 - 1; $m++ ) {	    	     
	     my $index = $a[rand @a] - 1;	     
	     push @sps1, $sp1[$index];
	     push @sps2, $sp2[$index];
	     push @sps3, $sp3[$index];
	     push @sps4, $sp4[$index];	      	   
	   }
	 my ($d,$abbasim, $babasim, $bbaasim) = calD (\@sps1, \@sps2, \@sps3, \@sps4);
	 push @bootd, $d;
       } #for (1 .. $boot)

       my $sdev = sprintf("%.4f", stddev (@bootd));
       my $mean =  sprintf("%.4f", mean (@bootd));
       my $zscore = 'nan';
       $zscore = sprintf("%.4f", $rawd /$sdev) unless $sdev == 0;
       my $p_value = "nan";
       $p_value = (1-ndtr($zscore))*2 unless $zscore eq 'nan'; ##two tailed
       
       #print OUT "file\tTotalNumberOfSites\tNumBBAA\tNumABBA\tNumBABA\traw_D\tZ-score\tD_sd\tp-value(two-tailed)\n";
       print OUT $lib, "\t", scalar @sp1 -1, "\t",$bbaa, "\t",  $abba, "\t", $baba, "\t", sprintf("%.4f", $rawd), "\t", $sdev, "\t", $zscore, "\t", sprintf("%.5f",$p_value), "\n";
     }#if ($boot && !$jack){
     
     elsif ($boot && $jack){
       print "Now doing jackknife resampling!\n";
       my @jackd;
       my $left = scalar @sp1 - 1 - $jack;
       my $half = (scalar @sp1 - 1)/2;
       die "the block size musy be shorter than half of the alignment length!\n" if $jack >= $half;
       die "the numer of replicates must be set smaller. Suggest less than $left \n" if $boot > $left;
     
       my @a;
       for (my $in = 1; $in <= $left; $in = $in + $left/$boot) {	
	 push @a, round ($in) if $in <= $left;
       }
	 
       my @sps1;
       my @sps2;
       my @sps3;
       my @sps4;
    
       for my $index (@a) {	    	     
	 push @sps1, @sp1[-($index+$jack) .. -$index];
	 push @sps2, @sp2[-($index+$jack) .. -$index];
	 push @sps3, @sp3[-($index+$jack) .. -$index];
	 push @sps4, @sp4[-($index+$jack) .. -$index];
	 my ($d,$abbasim, $babasim, $bbaasim) = calD (\@sps1, \@sps2, \@sps3, \@sps4);
	 push @jackd, $d;
       }

       my $sdev = sprintf("%.4f", stddev (@jackd));
       my $mean =  sprintf("%.4f", mean (@jackd));
       my $zscore = 'nan';
       $zscore = sprintf("%.4f", $rawd /$sdev) unless $sdev == 0;
       my $p_value = "nan";
       $p_value = (1-ndtr($zscore))*2 unless $zscore eq 'nan'; ##two tailed
       
       #print OUT "file\tTotalNumberOfSites\tNumBBAA\tNumABBA\tNumBABA\traw_D\tZ-score\tD_sd\tp-value(two-tailed)\n";
       print OUT $lib, "\t", scalar @sp1 -1, "\t",$bbaa, "\t",  $abba, "\t", $baba, "\t", sprintf("%.4f", $rawd), "\t", $sdev, "\t", $zscore, "\t", sprintf("%.5f",$p_value), "\n";
     }#if ($boot && !$jack){
     
     
   }#foreach (@aln)
    close OUT;
  } 

sub calD {
  my ($p1, $p2, $p3,$p4) = @_;
  my @sp1 = @{$p1};
  my @sp2 = @{$p2};
  my @sp3 = @{$p3};
  my @sp4 = @{$p4};

  my $abba = 0;
  my $baba = 0;
  my $bbaa = 0;
  my $rawd = 0;
  
  my %m = ("N" => 0,"?" => 0, "-" => 0, "Y"=>0, "S"=>0, "K"=>0, "M"=>0, "R"=>0, "W"=>0, "B"=>0, "D"=>0, "V"=>0, "H"=>0,);
  my %am =  ( "Y"=>0, "S"=>0, "K"=>0, "M"=>0, "R"=>0, "W"=>0);

  for (my $i = 0; $i < scalar @sp1 -1; $i++) {
    if ($m{$sp1[$i]} || $m{$sp2[$i]} || $m{$sp3[$i]} || $m{$sp4[$i]}) {
      next;
    }
    else {
      my @all = ($sp1[$i], $sp2[$i], $sp3[$i], $sp4[$i]); 
      if (scalar (uniq @all) == 2) {
	
	if ($sp1[$i] ne $sp2[$i]) {
	  if ($sp3[$i] ne $sp4[$i]) {
	    if ($sp1[$i] eq $sp3[$i]) {
	      $baba++;		 
	    }
	    if ($sp2[$i] eq $sp3[$i]) {
	      $abba++;		 
	    }
	  }
	}
	if ($sp1[$i] eq $sp2[$i]) {
	  if ($sp3[$i] ne $sp1[$i]) {
	    if ($sp3[$i] eq $sp4[$i]) {
	      $bbaa++;
	    }
	  }
	}
      }#if (scalar (uniq @all) == 2) {
      
    }#else      
  } #for (my $i = 0; $i < $len; $i++) {
  
  $rawd =  ($abba - $baba) / ($abba + $baba) if (($abba + $baba) > 0); 
  return ($rawd,$abba, $baba, $bbaa);
}

sub mk {
   die (qq/
seqCapture popstats mk [options]

Basic Options:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-f    DIR     A folder with fasta files with codon alignments of all samples
              including ingroup and outgroups
-i    CHAR    Ingroup sample ID 
-1    CHAR    First outgroup sample ID (required)
-2    CHAR    Second outgroup sample ID (needed for polarized analysis) [null]
-u    FILE    Outfile name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


example for a fasta file in the folder specified by -f. Sequence must not have line breaks:
        one line for id, one line for sequence for each sample.

	>ingroup; sample1
	ATGGCTAACGTAATTAAACCCGTTTTGACTTACCAGTTAGATGGCTCCAATCGTGATTTT
	>ingroup; sample2
	ATGGCTAACGTAATTAAAACCGTTTTGACTTACCAGTTAGATGGCTCCAATCGTGATTTT
	>ingroup; sample3
	ATGGCTAACGTAATTAAAACCGTTTTTACTTACCAGTTAGATGGCTCCAATCGTGATTTT
	>outgroup1; sample1
	ATGGCTAACGTAATTAACACCGTTTTGACTTACCAGTTAGATGGCTCCAATCGTGATTTT
        >outgroup2; sample2
        ATGGCTAACGTAATTAAAACCGTTTTGACTTACCAGTTAGATGGCTCCAATCGTCATTTT

e.g. for unpolarized analysis

	perl popStats.pl mk -u outfile  -1 outgroup1 -i ingroup -f A_folder_with_fasta_files

e.g. for polarized anlaysis

	perl popStats.pl mk -u outfile  -1 outgroup1 -2 outgroup2 -i ingroup -f A_folder_with_fasta_files
         
\n\n/) unless (@ARGV);

   my %opts = (f=>undef, i=>undef, 1=>undef, 2=>undef, u=>undef);
   getopts('f:i:1:2:u:', \%opts);

  
   my $dir = redir ($opts{f});
   my @aln = <$dir*fa>;

   my $o1 = $opts{1};
   my $o2;
   $o2 = $opts{2} if $opts{2};
   my $inid = $opts{i};

   my $out = $opts{u};
   open (OUT, ">", $out);
   print OUT "file\tnSyspoly\tnSysfix\tSyspoly\tSysfix\tMKcodons\tavgSample\tFETpval\n";

   
   my $mut = initialize();
   my %mut_code = %{$mut};
   
   foreach (@aln) {
     my $file = $_;
     my $lib = $1 if basename ($file) =~ /(\S+)\.fa/;
     open (my $filein , "<", $file);

     my @in;
     my @out;
     
     while (<$filein>) {
       chomp (my $line = $_);
       if ($line =~ /^>(\S+);/) {
	 my $id = $1;
	 chomp (my $seq = <$filein>);
	 if ($inid eq $id) {
	   push @in, $seq;	   
	 }
	 elsif ($o1 eq $id) {
	   push @out, $seq;
	 }
	 elsif ($o2 && $o2 eq $id) {
	   push @out, $seq;   
	 }
	 else {
	   next;
	 }
       }
     } 
     close $filein;    
     calcMK($lib, \@in, \@out, \%mut_code);         
   } ##foreach @aln

   close OUT;
}

sub calcMK{
  my ($lib, $SEQ, $outSEQ, $mut_code) = @_;
  my $outnum = @{$outSEQ};
  my %mut_code = %{$mut_code};
  
  my $POS=0;
  my $numNT = length $SEQ->[0];
  my $numSEQS = scalar @{$SEQ};
  my $numcodons=0;
  my $Nfix=0;
  my $Sfix=0;
  my $numIN=0;
  my $Npoly=0;
  my $Spoly=0;
  
  while($POS < $numNT){	## each window
    my @codons;
    my $ins=0;
    my $sameasOUT=0;
    my $difffromOUT=0;
    my $INsame=0; ## number of different sims will be scalar @codons
    my $OUTcodon1 = substr($outSEQ->[0], $POS, 3);
    my $OUT1AA = codon2aa($OUTcodon1);
    my $OUTcodon2 = $OUTcodon1;
    if($outnum == 2){
      $OUTcodon2 = substr($outSEQ->[1], $POS, 3);
    }
    if($OUTcodon1 eq $OUTcodon2 && $OUT1AA ne 'X'){
      for(my $x=0;$x<$numSEQS;$x++){
	my $tempcodon = substr($SEQ->[$x], $POS, 3);
	my $AAX = codon2aa($tempcodon);
	if($AAX ne 'X'){
	  $ins++;
	  if($tempcodon eq $OUTcodon1){
	    $sameasOUT++;
	  }
	  elsif($codons[0]){ ## there is at least one other sim codon & this sim is diff from mel
	    my $flag=0;
	    $difffromOUT++;
	    my $a=0;
	    while(($a < scalar @codons) && ($flag != 2)){ ##only push each codon type on once
	      if($codons[$a] ne $tempcodon){
		$flag=1;
	      }
	      elsif($codons[$a] eq $tempcodon){
		$flag=2;
		$INsame++;
	      }
	      $a++;
	    }
	    if($flag == 1){
	      push(@codons, $tempcodon);
	    }
	  }
	  else{ ## if 1st codon is different from mel, push it on
	    $difffromOUT++;
	    push(@codons, $tempcodon);
	  }
	}
      }
      if($ins > 1 && scalar @codons){
	if($difffromOUT){
	  if($INsame == $ins-1){	## fixed diff from mel
	    my $pair = $codons[0].$OUTcodon1;
	    my ($tNdiff, $tSdiff) = codon2codon($pair, \%mut_code);
	    $Nfix += $tNdiff;
	    $Sfix += $tSdiff;
	    #print "$tNdiff\t$tSdiff\n";
	  }
	  elsif(scalar @codons <= 2){	## polymorphic site, can have 2 codon types for sim
	    if($sameasOUT){	## 2 sim codons, 1 sameasOUT, 1 diff, calc path for diff add to poly
	      my $tempNdiff=3;
	      my $tempSdiff=0;
	      for(my $y=0;$y < scalar @codons;$y++){	## takes simple route of evol from mel to a sim that minimizes NS subs
		my $pair = $codons[$y].$OUTcodon1;
		my ($tNdiff, $tSdiff) = codon2codon($pair,\%mut_code);
		if($tNdiff < $tempNdiff || ($tNdiff == $tempNdiff && $tSdiff <= $tempSdiff)){
		  $tempNdiff = $tNdiff;
		  $tempSdiff = $tSdiff;
		}
	      }
	      ## now compare the two sim codons if there are 3 total sim codons
	      if(scalar @codons == 2){
		my $pair = $codons[0].$codons[1];
		my ($tNdiff, $tSdiff) = codon2codon($pair,\%mut_code);
		$tempNdiff += $tNdiff;
		$tempSdiff += $tSdiff;
	      }							
	      $Npoly += $tempNdiff;
	      $Spoly += $tempSdiff;
	    }
	    else{	## 2 sim codons, both diff from mel
	      my $tempNfix=3;
	      my $tempSfix=0;
	      my $tempNpoly=3;
	      my $tempSpoly=0;
	      for(my $z=0;$z < scalar @codons;$z++){
		## calc path btn mel - sim1 - sim2, then mel - sim2 - sim1
		my $pair = $codons[$z].$OUTcodon1;
		my ($tNdiff, $tSdiff) = codon2codon($pair,\%mut_code);
		## use min NS change as path for fixed counts
		if($tNdiff < $tempNfix || ($tNdiff == $tempNfix && $tSdiff <= $tempSfix)){
		  $tempNfix = $tNdiff;
		  $tempSfix = $tSdiff;
		}
	      }
	      $Nfix += $tempNfix;
	      $Sfix += $tempSfix;
	      
	      ## calc path btn 2 sim alleles
	      ## use min NS change as path for poly counts							
	      my $polypair = $codons[0].$codons[1];
	      my ($tNp, $tSp) = codon2codon($polypair,\%mut_code);
	      $Npoly += $tNp;
	      $Spoly += $tSp;
	    }
	  }
	}
	
	#print $ins, "\t", $numIN, "\t", $numcodons, "\n";
      }
    }
    if($ins>1){
      $numcodons++; ## if 2 or more alleles/codon
      $numIN += $ins; ## total number of sim alleles over the sequence
    }
    $POS += 3;
  }
  
  if($Npoly || $Spoly || $Nfix || $Sfix){		
    
    my $avgsamplesize = $numIN/$numcodons;
    my $FET_pvalue = calculateStatistic (n11=>$Npoly,
					 n1p=>$Npoly+$Spoly,
					 np1=>$Npoly+$Nfix,
					 npp=>$Npoly+$Nfix+$Spoly+$Sfix);
    
    print OUT $lib, "\t", $Npoly, "\t", $Nfix, "\t", $Spoly, "\t", $Sfix, "\t", $numcodons, "\t";
    printf OUT "%.2f\t", $avgsamplesize; 
    print OUT $FET_pvalue, "\n";
  }
  
}



sub stats {
  
  die (qq/
seqCapture popstats stats [options]

Basic Options:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-a    DIR     Folder with all alignments of ingroups
-o    DIR     Results folder
-b    CHAR    Outgroup sample ID if you have one [null]
              the sequence of outgroup sample must be 
              present in the alignments of ingroups 
-s    FILE    A text file contains samples to be summarized
              (one sample ID per line). If no sample file is 
              provided, the script will calculate these stats
              for all samples

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      
\n\n/) unless (@ARGV);
  
  
  my %opts = (a=>undef, s=>undef, o=>undef, b=>undef);
  getopts('a:s:o:b:', \%opts);

  my $alndir = redir ($opts{a}) || die "can not open alignment folder!\n";
  my $resdir = redir ($opts{o}) || die "must provide results folder!\n";
  mkdir $resdir unless -e $resdir;
  
  my $name = $opts{s} || die "must provide a file with sample IDs!\n";;

  my $outgroup;
  $outgroup = $opts{b} if $opts{b};

  my %names;
  open (NAME, "<", $name);
  while (<NAME>) {
    chomp (my $line = $_);
    $names{$line}++; 
  }
  close NAME;
  
  
  my @data = <$alndir*.aln>;
  
  my %all_aln;
  my %statsind;
  my $outseq;
  
  foreach (@data) {
    my $file = $_;
    my $lib = $1 if basename ($file) =~ /(\S+)\.aln/;
    my $alningroup = $resdir . $lib .".aln";

    my $alnoutgroup;
    if ($outgroup) { 
      $alnoutgroup = $resdir . $lib .".outgroup.aln";
      open (OUTG, ">", $alnoutgroup);
      print OUTG "CLUSTAL W", "\n\n";
    }
    
    open (IN, "<", $file);
    open (OUT, ">", $alningroup);
    print OUT "CLUSTAL W", "\n\n";
    
    my $len;
    my $outpre = 0;
    my $samsize = 0;
    while (<IN>) {
      chomp (my $line = $_);
      if ($line =~ /^>(\S+)/) {
	my $id = $1;
	chomp (my $seq = <IN>);
	
	if ($names{$id}) {
	  $all_aln{$id} .= $seq;
	  print OUT $id, "\t", $seq,"\n" if $seq =~ /A|G|C|T/;
	  $len = length ($seq) if $seq =~ /A|G|C|T/;
	  $samsize++ if $seq =~ /A|G|C|T/;
	}
	elsif ($outgroup && $outgroup eq $id) {
	  $outseq .= $seq;
	  print OUTG $id, "\t", $seq,"\n";
	  $outpre = 1 if $seq =~ /A|G|C|T/;
	}
	else {
	  next;
	}	
      }
    }
    close IN;
    close OUT;
    close OUTG if ($outgroup);

    my $D;
    if ($outgroup && $outpre == 1) {
      $D = fuliD ($alningroup, $alnoutgroup);
    }
    else {
      $D = "NA";
    }
    
    my ($pi, $theta, $tajima_D, $seg)  = thetas ($alningroup, $len, $samsize); 
    $statsind{$lib}{'pi'} = $pi;
    $statsind{$lib}{'theta'} = $theta;
    $statsind{$lib}{'D'} = $D;
    $statsind{$lib}{'td'} = $tajima_D;
    $statsind{$lib}{'seg'} = $seg;
    
    unlink ($alningroup);
    unlink ($alnoutgroup) if $alnoutgroup;
    
  } #foreach @data

  
  my $allalnout;
  if ($outgroup) { 
    $allalnout = $resdir . "alloutgroup.aln";
    open (OUTGA, ">", $allalnout);
    print OUTGA "CLUSTAL W", "\n\n";
    print OUTGA $outgroup, "\t", $outseq, "\n";
    close OUTGA;
  }
  
  my $allalnin =  $resdir . "allingroup.aln";
  open (OUTINALL, ">", $allalnin);
  print OUTINALL "CLUSTAL W", "\n\n";
  
  my $lenall;
  my $samsizeall = 0;
  foreach my $id (sort {$a cmp $b} keys %all_aln) {
    $samsizeall++;
    $lenall = length ($all_aln{$id});
    print OUTINALL $id, "\t", $all_aln{$id},"\n";
  }
  close OUTINALL;
 
  my $Dall;
  if ($outgroup) {
    $Dall = fuliD ($allalnin, $allalnout);
  }
  else {
    $Dall = "NA";
  }
  
  my ($pi2, $theta2, $tajima_D2, $seg2)  = thetas ($allalnin, $lenall,$samsizeall); 
  unlink ($allalnin);
  unlink ($allalnout) if  ($outgroup);
  
  my $finalout = $resdir . "finalstats.txt";
  open (FINAL, ">", $finalout);
  print FINAL "marker\tpi\thetaW\tTajima\tFuLiD\tNumSeg\n";
  print FINAL "global", "\t", $pi2, "\t", $theta2, "\t", $tajima_D2, "\t", $Dall, "\t", $seg2, "\n";
  
  foreach my $lib (sort {$a cmp $b} keys %statsind) {
    print FINAL $lib, "\t",   $statsind{$lib}{'pi'}, "\t",  $statsind{$lib}{'theta'}, "\t", $statsind{$lib}{'td'}, "\t", $statsind{$lib}{'D'}, "\t", $statsind{$lib}{'seg'}, "\n";    
  }
  close FINAL;
  
}

sub fuliD {
  my ($in, $out) = @_;
  my $parser = Bio::AlignIO->new(-format => 'clustalw', -file   => $in);
  my $aln = $parser ->next_aln;
  my $pop = Bio::PopGen::Utilities->aln_to_population(-alignment => $aln);
  
  my $parser1 = Bio::AlignIO->new(-format => 'clustalw', -file   => $out);
  my $aln1 = $parser1->next_aln;
  my $pop1 = Bio::PopGen::Utilities->aln_to_population(-alignment => $aln1);
  my $D = Bio::PopGen::Statistics->fu_and_li_D($pop,$pop1);
  return $D;
}

sub thetas {
  my ($in, $len, $samplesize) = @_;
  
  my $parser = Bio::AlignIO->new(-format => 'clustalw', -file => $in);
  my $aln = $parser ->next_aln;
  my $pop = Bio::PopGen::Utilities->aln_to_population(-alignment => $aln);

  my $seg = Bio::PopGen::Statistics->segregating_sites_count($pop);
  my $pi     = Bio::PopGen::Statistics->pi($pop,$len);
  my $theta  = Bio::PopGen::Statistics->theta($samplesize,$seg,$len);
  my $tajima_D  = Bio::PopGen::Statistics->tajima_D($pop);
  #my $het  = Bio::PopGen::Statistics->heterozygosity(9, 0.1);
  
  return ($pi, $theta, $tajima_D,$seg);
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


sub het {
  my ($samplesize, $freq1) = @_;
  my $freq2 = 1-$freq1;
  my $pq = ($freq1**2) + (($freq2)**2); ##p2+q2
  my $het = ((1-$pq)*$samplesize)/($samplesize-1) ;
  return ($het);
}


sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => 'X',    # Stop
    'TAG' => 'X',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => 'X',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    
    );

    if(exists $genetic_code{$codon}) {
    	return $genetic_code{$codon};
    }else{
		return 'X';
    }
}



############### codon2codon ####################

sub codon2codon {
  my($codonpair, $mut_code) = @_;
  my %mut_code = %{$mut_code};
  
  $codonpair = uc $codonpair;
  
  if(exists $mut_code{$codonpair}) {
    my $pair = $mut_code{$codonpair};
    return split(/,/,$pair);
  }
  else{
    print OUT "Bad codon pair \"$codonpair\"!!\n";
    exit;
  }
}


############### initialization of codon2codon hash ###################

sub initialize{
  my %m = (
	'AAAAAA'	=>	'0,0',
	'AAAAAC'	=>	'1,0',
	'AAAAAG'	=>	'0,1',
	'AAAAAT'	=>	'1,0',
	'AAAACA'	=>	'1,0',
	'AAAACC'	=>	'1,1',
	'AAAACG'	=>	'1,1',
	'AAAACT'	=>	'1,1',
	'AAAAGA'	=>	'1,0',
	'AAAAGC'	=>	'2,0',
	'AAAAGG'	=>	'1,1',
	'AAAAGT'	=>	'2,0',
	'AAAATA'	=>	'1,0',
	'AAAATC'	=>	'1,1',
	'AAAATG'	=>	'1,1',
	'AAAATT'	=>	'1,1',
	'AAACAA'	=>	'1,0',
	'AAACAC'	=>	'2,0',
	'AAACAG'	=>	'1,1',
	'AAACAT'	=>	'2,0',
	'AAACCA'	=>	'2,0',
	'AAACCC'	=>	'2,1',
	'AAACCG'	=>	'2,1',
	'AAACCT'	=>	'2,1',
	'AAACGA'	=>	'1,1',
	'AAACGC'	=>	'1,2',
	'AAACGG'	=>	'1,2',
	'AAACGT'	=>	'1,2',
	'AAACTA'	=>	'2,0',
	'AAACTC'	=>	'2,1',
	'AAACTG'	=>	'2,1',
	'AAACTT'	=>	'2,1',
	'AAAGAA'	=>	'1,0',
	'AAAGAC'	=>	'2,0',
	'AAAGAG'	=>	'1,1',
	'AAAGAT'	=>	'2,0',
	'AAAGCA'	=>	'2,0',
	'AAAGCC'	=>	'2,1',
	'AAAGCG'	=>	'2,1',
	'AAAGCT'	=>	'2,1',
	'AAAGGA'	=>	'2,0',
	'AAAGGC'	=>	'2,1',
	'AAAGGG'	=>	'2,1',
	'AAAGGT'	=>	'2,1',
	'AAAGTA'	=>	'2,0',
	'AAAGTC'	=>	'2,1',
	'AAAGTG'	=>	'2,1',
	'AAAGTT'	=>	'2,1',
	'AAATAC'	=>	'2,0',
	'AAATAT'	=>	'2,0',
	'AAATCA'	=>	'2,0',
	'AAATCC'	=>	'2,1',
	'AAATCG'	=>	'2,1',
	'AAATCT'	=>	'2,1',
	'AAATGC'	=>	'3,0',
	'AAATGG'	=>	'2,1',
	'AAATGT'	=>	'3,0',
	'AAATTA'	=>	'2,0',
	'AAATTC'	=>	'2,1',
	'AAATTG'	=>	'2,1',
	'AAATTT'	=>	'2,1',
	'AACAAA'	=>	'1,0',
	'AACAAC'	=>	'0,0',
	'AACAAG'	=>	'1,0',
	'AACAAT'	=>	'0,1',
	'AACACA'	=>	'1,1',
	'AACACC'	=>	'1,0',
	'AACACG'	=>	'1,1',
	'AACACT'	=>	'1,1',
	'AACAGA'	=>	'2,0',
	'AACAGC'	=>	'1,0',
	'AACAGG'	=>	'2,0',
	'AACAGT'	=>	'1,1',
	'AACATA'	=>	'1,1',
	'AACATC'	=>	'1,0',
	'AACATG'	=>	'2,0',
	'AACATT'	=>	'1,1',
	'AACCAA'	=>	'2,0',
	'AACCAC'	=>	'1,0',
	'AACCAG'	=>	'2,0',
	'AACCAT'	=>	'1,1',
	'AACCCA'	=>	'2,1',
	'AACCCC'	=>	'2,0',
	'AACCCG'	=>	'2,1',
	'AACCCT'	=>	'2,1',
	'AACCGA'	=>	'2,1',
	'AACCGC'	=>	'2,0',
	'AACCGG'	=>	'2,1',
	'AACCGT'	=>	'2,1',
	'AACCTA'	=>	'2,1',
	'AACCTC'	=>	'2,0',
	'AACCTG'	=>	'2,1',
	'AACCTT'	=>	'2,1',
	'AACGAA'	=>	'2,0',
	'AACGAC'	=>	'1,0',
	'AACGAG'	=>	'2,0',
	'AACGAT'	=>	'1,1',
	'AACGCA'	=>	'2,1',
	'AACGCC'	=>	'2,0',
	'AACGCG'	=>	'2,1',
	'AACGCT'	=>	'2,1',
	'AACGGA'	=>	'2,1',
	'AACGGC'	=>	'2,0',
	'AACGGG'	=>	'2,1',
	'AACGGT'	=>	'2,1',
	'AACGTA'	=>	'2,1',
	'AACGTC'	=>	'2,0',
	'AACGTG'	=>	'2,1',
	'AACGTT'	=>	'2,1',
	'AACTAC'	=>	'1,0',
	'AACTAT'	=>	'1,1',
	'AACTCA'	=>	'2,1',
	'AACTCC'	=>	'2,0',
	'AACTCG'	=>	'2,1',
	'AACTCT'	=>	'2,1',
	'AACTGC'	=>	'2,0',
	'AACTGG'	=>	'3,0',
	'AACTGT'	=>	'2,1',
	'AACTTA'	=>	'2,1',
	'AACTTC'	=>	'2,0',
	'AACTTG'	=>	'3,0',
	'AACTTT'	=>	'2,1',
	'AAGAAA'	=>	'0,1',
	'AAGAAC'	=>	'1,0',
	'AAGAAG'	=>	'0,0',
	'AAGAAT'	=>	'1,0',
	'AAGACA'	=>	'1,1',
	'AAGACC'	=>	'1,1',
	'AAGACG'	=>	'1,0',
	'AAGACT'	=>	'1,1',
	'AAGAGA'	=>	'1,1',
	'AAGAGC'	=>	'2,0',
	'AAGAGG'	=>	'1,0',
	'AAGAGT'	=>	'2,0',
	'AAGATA'	=>	'1,1',
	'AAGATC'	=>	'2,0',
	'AAGATG'	=>	'1,0',
	'AAGATT'	=>	'2,0',
	'AAGCAA'	=>	'1,1',
	'AAGCAC'	=>	'2,0',
	'AAGCAG'	=>	'1,0',
	'AAGCAT'	=>	'2,0',
	'AAGCCA'	=>	'2,1',
	'AAGCCC'	=>	'2,1',
	'AAGCCG'	=>	'2,0',
	'AAGCCT'	=>	'2,1',
	'AAGCGA'	=>	'1,2',
	'AAGCGC'	=>	'1,2',
	'AAGCGG'	=>	'1,1',
	'AAGCGT'	=>	'1,2',
	'AAGCTA'	=>	'2,1',
	'AAGCTC'	=>	'2,1',
	'AAGCTG'	=>	'2,0',
	'AAGCTT'	=>	'2,1',
	'AAGGAA'	=>	'1,1',
	'AAGGAC'	=>	'2,0',
	'AAGGAG'	=>	'1,0',
	'AAGGAT'	=>	'2,0',
	'AAGGCA'	=>	'2,1',
	'AAGGCC'	=>	'2,1',
	'AAGGCG'	=>	'2,0',
	'AAGGCT'	=>	'2,1',
	'AAGGGA'	=>	'2,1',
	'AAGGGC'	=>	'2,1',
	'AAGGGG'	=>	'2,0',
	'AAGGGT'	=>	'2,1',
	'AAGGTA'	=>	'2,1',
	'AAGGTC'	=>	'2,1',
	'AAGGTG'	=>	'2,0',
	'AAGGTT'	=>	'2,1',
	'AAGTAC'	=>	'2,0',
	'AAGTAT'	=>	'2,0',
	'AAGTCA'	=>	'2,1',
	'AAGTCC'	=>	'2,1',
	'AAGTCG'	=>	'2,0',
	'AAGTCT'	=>	'2,1',
	'AAGTGC'	=>	'3,0',
	'AAGTGG'	=>	'2,0',
	'AAGTGT'	=>	'3,0',
	'AAGTTA'	=>	'2,1',
	'AAGTTC'	=>	'3,0',
	'AAGTTG'	=>	'2,0',
	'AAGTTT'	=>	'3,0',
	'AATAAA'	=>	'1,0',
	'AATAAC'	=>	'0,1',
	'AATAAG'	=>	'1,0',
	'AATAAT'	=>	'0,0',
	'AATACA'	=>	'1,1',
	'AATACC'	=>	'1,1',
	'AATACG'	=>	'1,1',
	'AATACT'	=>	'1,0',
	'AATAGA'	=>	'2,0',
	'AATAGC'	=>	'1,1',
	'AATAGG'	=>	'2,0',
	'AATAGT'	=>	'1,0',
	'AATATA'	=>	'1,1',
	'AATATC'	=>	'1,1',
	'AATATG'	=>	'2,0',
	'AATATT'	=>	'1,0',
	'AATCAA'	=>	'2,0',
	'AATCAC'	=>	'1,1',
	'AATCAG'	=>	'2,0',
	'AATCAT'	=>	'1,0',
	'AATCCA'	=>	'2,1',
	'AATCCC'	=>	'2,1',
	'AATCCG'	=>	'2,1',
	'AATCCT'	=>	'2,0',
	'AATCGA'	=>	'2,1',
	'AATCGC'	=>	'2,1',
	'AATCGG'	=>	'2,1',
	'AATCGT'	=>	'2,0',
	'AATCTA'	=>	'2,1',
	'AATCTC'	=>	'2,1',
	'AATCTG'	=>	'2,1',
	'AATCTT'	=>	'2,0',
	'AATGAA'	=>	'2,0',
	'AATGAC'	=>	'1,1',
	'AATGAG'	=>	'2,0',
	'AATGAT'	=>	'1,0',
	'AATGCA'	=>	'2,1',
	'AATGCC'	=>	'2,1',
	'AATGCG'	=>	'2,1',
	'AATGCT'	=>	'2,0',
	'AATGGA'	=>	'2,1',
	'AATGGC'	=>	'2,1',
	'AATGGG'	=>	'2,1',
	'AATGGT'	=>	'2,0',
	'AATGTA'	=>	'2,1',
	'AATGTC'	=>	'2,1',
	'AATGTG'	=>	'2,1',
	'AATGTT'	=>	'2,0',
	'AATTAC'	=>	'1,1',
	'AATTAT'	=>	'1,0',
	'AATTCA'	=>	'2,1',
	'AATTCC'	=>	'2,1',
	'AATTCG'	=>	'2,1',
	'AATTCT'	=>	'2,0',
	'AATTGC'	=>	'2,1',
	'AATTGG'	=>	'3,0',
	'AATTGT'	=>	'2,0',
	'AATTTA'	=>	'2,1',
	'AATTTC'	=>	'2,1',
	'AATTTG'	=>	'3,0',
	'AATTTT'	=>	'2,0',
	'ACAAAA'	=>	'1,0',
	'ACAAAC'	=>	'1,1',
	'ACAAAG'	=>	'1,1',
	'ACAAAT'	=>	'1,1',
	'ACAACA'	=>	'0,0',
	'ACAACC'	=>	'0,1',
	'ACAACG'	=>	'0,1',
	'ACAACT'	=>	'0,1',
	'ACAAGA'	=>	'1,0',
	'ACAAGC'	=>	'1,1',
	'ACAAGG'	=>	'1,1',
	'ACAAGT'	=>	'1,1',
	'ACAATA'	=>	'1,0',
	'ACAATC'	=>	'1,1',
	'ACAATG'	=>	'1,1',
	'ACAATT'	=>	'1,1',
	'ACACAA'	=>	'2,0',
	'ACACAC'	=>	'2,1',
	'ACACAG'	=>	'2,1',
	'ACACAT'	=>	'2,1',
	'ACACCA'	=>	'1,0',
	'ACACCC'	=>	'1,1',
	'ACACCG'	=>	'1,1',
	'ACACCT'	=>	'1,1',
	'ACACGA'	=>	'1,1',
	'ACACGC'	=>	'1,2',
	'ACACGG'	=>	'1,2',
	'ACACGT'	=>	'1,2',
	'ACACTA'	=>	'2,0',
	'ACACTC'	=>	'2,1',
	'ACACTG'	=>	'2,1',
	'ACACTT'	=>	'2,1',
	'ACAGAA'	=>	'2,0',
	'ACAGAC'	=>	'2,1',
	'ACAGAG'	=>	'2,1',
	'ACAGAT'	=>	'2,1',
	'ACAGCA'	=>	'1,0',
	'ACAGCC'	=>	'1,1',
	'ACAGCG'	=>	'1,1',
	'ACAGCT'	=>	'1,1',
	'ACAGGA'	=>	'2,0',
	'ACAGGC'	=>	'2,1',
	'ACAGGG'	=>	'2,1',
	'ACAGGT'	=>	'2,1',
	'ACAGTA'	=>	'2,0',
	'ACAGTC'	=>	'2,1',
	'ACAGTG'	=>	'2,1',
	'ACAGTT'	=>	'2,1',
	'ACATAC'	=>	'2,1',
	'ACATAT'	=>	'2,1',
	'ACATCA'	=>	'1,0',
	'ACATCC'	=>	'1,1',
	'ACATCG'	=>	'1,1',
	'ACATCT'	=>	'1,1',
	'ACATGC'	=>	'2,1',
	'ACATGG'	=>	'2,1',
	'ACATGT'	=>	'2,1',
	'ACATTA'	=>	'2,0',
	'ACATTC'	=>	'2,1',
	'ACATTG'	=>	'2,1',
	'ACATTT'	=>	'2,1',
	'ACCAAA'	=>	'1,1',
	'ACCAAC'	=>	'1,0',
	'ACCAAG'	=>	'1,1',
	'ACCAAT'	=>	'1,1',
	'ACCACA'	=>	'0,1',
	'ACCACC'	=>	'0,0',
	'ACCACG'	=>	'0,1',
	'ACCACT'	=>	'0,1',
	'ACCAGA'	=>	'1,1',
	'ACCAGC'	=>	'1,0',
	'ACCAGG'	=>	'1,1',
	'ACCAGT'	=>	'1,1',
	'ACCATA'	=>	'1,1',
	'ACCATC'	=>	'1,0',
	'ACCATG'	=>	'1,1',
	'ACCATT'	=>	'1,1',
	'ACCCAA'	=>	'2,1',
	'ACCCAC'	=>	'2,0',
	'ACCCAG'	=>	'2,1',
	'ACCCAT'	=>	'2,1',
	'ACCCCA'	=>	'1,1',
	'ACCCCC'	=>	'1,0',
	'ACCCCG'	=>	'1,1',
	'ACCCCT'	=>	'1,1',
	'ACCCGA'	=>	'1,2',
	'ACCCGC'	=>	'2,0',
	'ACCCGG'	=>	'1,2',
	'ACCCGT'	=>	'2,1',
	'ACCCTA'	=>	'2,1',
	'ACCCTC'	=>	'2,0',
	'ACCCTG'	=>	'2,1',
	'ACCCTT'	=>	'2,1',
	'ACCGAA'	=>	'2,1',
	'ACCGAC'	=>	'2,0',
	'ACCGAG'	=>	'2,1',
	'ACCGAT'	=>	'2,1',
	'ACCGCA'	=>	'1,1',
	'ACCGCC'	=>	'1,0',
	'ACCGCG'	=>	'1,1',
	'ACCGCT'	=>	'1,1',
	'ACCGGA'	=>	'2,1',
	'ACCGGC'	=>	'2,0',
	'ACCGGG'	=>	'2,1',
	'ACCGGT'	=>	'2,1',
	'ACCGTA'	=>	'2,1',
	'ACCGTC'	=>	'2,0',
	'ACCGTG'	=>	'2,1',
	'ACCGTT'	=>	'2,1',
	'ACCTAC'	=>	'2,0',
	'ACCTAT'	=>	'2,1',
	'ACCTCA'	=>	'1,1',
	'ACCTCC'	=>	'1,0',
	'ACCTCG'	=>	'1,1',
	'ACCTCT'	=>	'1,1',
	'ACCTGC'	=>	'2,0',
	'ACCTGG'	=>	'2,1',
	'ACCTGT'	=>	'2,1',
	'ACCTTA'	=>	'2,1',
	'ACCTTC'	=>	'2,0',
	'ACCTTG'	=>	'2,1',
	'ACCTTT'	=>	'2,1',
	'ACGAAA'	=>	'1,1',
	'ACGAAC'	=>	'1,1',
	'ACGAAG'	=>	'1,0',
	'ACGAAT'	=>	'1,1',
	'ACGACA'	=>	'0,1',
	'ACGACC'	=>	'0,1',
	'ACGACG'	=>	'0,0',
	'ACGACT'	=>	'0,1',
	'ACGAGA'	=>	'1,1',
	'ACGAGC'	=>	'1,1',
	'ACGAGG'	=>	'1,0',
	'ACGAGT'	=>	'1,1',
	'ACGATA'	=>	'1,1',
	'ACGATC'	=>	'1,1',
	'ACGATG'	=>	'1,0',
	'ACGATT'	=>	'1,1',
	'ACGCAA'	=>	'2,1',
	'ACGCAC'	=>	'2,1',
	'ACGCAG'	=>	'2,0',
	'ACGCAT'	=>	'2,1',
	'ACGCCA'	=>	'1,1',
	'ACGCCC'	=>	'1,1',
	'ACGCCG'	=>	'1,0',
	'ACGCCT'	=>	'1,1',
	'ACGCGA'	=>	'1,2',
	'ACGCGC'	=>	'1,2',
	'ACGCGG'	=>	'1,1',
	'ACGCGT'	=>	'1,2',
	'ACGCTA'	=>	'2,1',
	'ACGCTC'	=>	'2,1',
	'ACGCTG'	=>	'2,0',
	'ACGCTT'	=>	'2,1',
	'ACGGAA'	=>	'2,1',
	'ACGGAC'	=>	'2,1',
	'ACGGAG'	=>	'2,0',
	'ACGGAT'	=>	'2,1',
	'ACGGCA'	=>	'1,1',
	'ACGGCC'	=>	'1,1',
	'ACGGCG'	=>	'1,0',
	'ACGGCT'	=>	'1,1',
	'ACGGGA'	=>	'2,1',
	'ACGGGC'	=>	'2,1',
	'ACGGGG'	=>	'2,0',
	'ACGGGT'	=>	'2,1',
	'ACGGTA'	=>	'2,1',
	'ACGGTC'	=>	'2,1',
	'ACGGTG'	=>	'2,0',
	'ACGGTT'	=>	'2,1',
	'ACGTAC'	=>	'2,1',
	'ACGTAT'	=>	'2,1',
	'ACGTCA'	=>	'1,1',
	'ACGTCC'	=>	'1,1',
	'ACGTCG'	=>	'1,0',
	'ACGTCT'	=>	'1,1',
	'ACGTGC'	=>	'2,1',
	'ACGTGG'	=>	'2,0',
	'ACGTGT'	=>	'2,1',
	'ACGTTA'	=>	'2,1',
	'ACGTTC'	=>	'2,1',
	'ACGTTG'	=>	'2,0',
	'ACGTTT'	=>	'2,1',
	'ACTAAA'	=>	'1,1',
	'ACTAAC'	=>	'1,1',
	'ACTAAG'	=>	'1,1',
	'ACTAAT'	=>	'1,0',
	'ACTACA'	=>	'0,1',
	'ACTACC'	=>	'0,1',
	'ACTACG'	=>	'0,1',
	'ACTACT'	=>	'0,0',
	'ACTAGA'	=>	'1,1',
	'ACTAGC'	=>	'1,1',
	'ACTAGG'	=>	'1,1',
	'ACTAGT'	=>	'1,0',
	'ACTATA'	=>	'1,1',
	'ACTATC'	=>	'1,1',
	'ACTATG'	=>	'1,1',
	'ACTATT'	=>	'1,0',
	'ACTCAA'	=>	'2,1',
	'ACTCAC'	=>	'2,1',
	'ACTCAG'	=>	'2,1',
	'ACTCAT'	=>	'2,0',
	'ACTCCA'	=>	'1,1',
	'ACTCCC'	=>	'1,1',
	'ACTCCG'	=>	'1,1',
	'ACTCCT'	=>	'1,0',
	'ACTCGA'	=>	'1,2',
	'ACTCGC'	=>	'2,1',
	'ACTCGG'	=>	'1,2',
	'ACTCGT'	=>	'2,0',
	'ACTCTA'	=>	'2,1',
	'ACTCTC'	=>	'2,1',
	'ACTCTG'	=>	'2,1',
	'ACTCTT'	=>	'2,0',
	'ACTGAA'	=>	'2,1',
	'ACTGAC'	=>	'2,1',
	'ACTGAG'	=>	'2,1',
	'ACTGAT'	=>	'2,0',
	'ACTGCA'	=>	'1,1',
	'ACTGCC'	=>	'1,1',
	'ACTGCG'	=>	'1,1',
	'ACTGCT'	=>	'1,0',
	'ACTGGA'	=>	'2,1',
	'ACTGGC'	=>	'2,1',
	'ACTGGG'	=>	'2,1',
	'ACTGGT'	=>	'2,0',
	'ACTGTA'	=>	'2,1',
	'ACTGTC'	=>	'2,1',
	'ACTGTG'	=>	'2,1',
	'ACTGTT'	=>	'2,0',
	'ACTTAC'	=>	'2,1',
	'ACTTAT'	=>	'2,0',
	'ACTTCA'	=>	'1,1',
	'ACTTCC'	=>	'1,1',
	'ACTTCG'	=>	'1,1',
	'ACTTCT'	=>	'1,0',
	'ACTTGC'	=>	'2,1',
	'ACTTGG'	=>	'2,1',
	'ACTTGT'	=>	'2,0',
	'ACTTTA'	=>	'2,1',
	'ACTTTC'	=>	'2,1',
	'ACTTTG'	=>	'2,1',
	'ACTTTT'	=>	'2,0',
	'AGAAAA'	=>	'1,0',
	'AGAAAC'	=>	'2,0',
	'AGAAAG'	=>	'1,1',
	'AGAAAT'	=>	'2,0',
	'AGAACA'	=>	'1,0',
	'AGAACC'	=>	'1,1',
	'AGAACG'	=>	'1,1',
	'AGAACT'	=>	'1,1',
	'AGAAGA'	=>	'0,0',
	'AGAAGC'	=>	'1,0',
	'AGAAGG'	=>	'0,1',
	'AGAAGT'	=>	'1,0',
	'AGAATA'	=>	'1,0',
	'AGAATC'	=>	'1,1',
	'AGAATG'	=>	'1,1',
	'AGAATT'	=>	'1,1',
	'AGACAA'	=>	'1,1',
	'AGACAC'	=>	'1,2',
	'AGACAG'	=>	'1,2',
	'AGACAT'	=>	'1,2',
	'AGACCA'	=>	'1,1',
	'AGACCC'	=>	'1,2',
	'AGACCG'	=>	'1,2',
	'AGACCT'	=>	'1,2',
	'AGACGA'	=>	'0,1',
	'AGACGC'	=>	'0,2',
	'AGACGG'	=>	'0,2',
	'AGACGT'	=>	'0,2',
	'AGACTA'	=>	'1,1',
	'AGACTC'	=>	'1,2',
	'AGACTG'	=>	'1,2',
	'AGACTT'	=>	'1,2',
	'AGAGAA'	=>	'2,0',
	'AGAGAC'	=>	'2,1',
	'AGAGAG'	=>	'2,1',
	'AGAGAT'	=>	'2,1',
	'AGAGCA'	=>	'2,0',
	'AGAGCC'	=>	'2,1',
	'AGAGCG'	=>	'2,1',
	'AGAGCT'	=>	'2,1',
	'AGAGGA'	=>	'1,0',
	'AGAGGC'	=>	'1,1',
	'AGAGGG'	=>	'1,1',
	'AGAGGT'	=>	'1,1',
	'AGAGTA'	=>	'2,0',
	'AGAGTC'	=>	'2,1',
	'AGAGTG'	=>	'2,1',
	'AGAGTT'	=>	'2,1',
	'AGATAC'	=>	'3,0',
	'AGATAT'	=>	'3,0',
	'AGATCA'	=>	'2,0',
	'AGATCC'	=>	'2,1',
	'AGATCG'	=>	'2,1',
	'AGATCT'	=>	'2,1',
	'AGATGC'	=>	'2,0',
	'AGATGG'	=>	'1,1',
	'AGATGT'	=>	'2,0',
	'AGATTA'	=>	'2,0',
	'AGATTC'	=>	'2,1',
	'AGATTG'	=>	'2,1',
	'AGATTT'	=>	'2,1',
	'AGCAAA'	=>	'2,0',
	'AGCAAC'	=>	'1,0',
	'AGCAAG'	=>	'2,0',
	'AGCAAT'	=>	'1,1',
	'AGCACA'	=>	'1,1',
	'AGCACC'	=>	'1,0',
	'AGCACG'	=>	'1,1',
	'AGCACT'	=>	'1,1',
	'AGCAGA'	=>	'1,0',
	'AGCAGC'	=>	'0,0',
	'AGCAGG'	=>	'1,0',
	'AGCAGT'	=>	'0,1',
	'AGCATA'	=>	'1,1',
	'AGCATC'	=>	'1,0',
	'AGCATG'	=>	'2,0',
	'AGCATT'	=>	'1,1',
	'AGCCAA'	=>	'2,1',
	'AGCCAC'	=>	'2,0',
	'AGCCAG'	=>	'2,1',
	'AGCCAT'	=>	'2,1',
	'AGCCCA'	=>	'2,1',
	'AGCCCC'	=>	'2,0',
	'AGCCCG'	=>	'2,1',
	'AGCCCT'	=>	'2,1',
	'AGCCGA'	=>	'1,1',
	'AGCCGC'	=>	'1,0',
	'AGCCGG'	=>	'1,1',
	'AGCCGT'	=>	'1,1',
	'AGCCTA'	=>	'2,1',
	'AGCCTC'	=>	'2,0',
	'AGCCTG'	=>	'2,1',
	'AGCCTT'	=>	'2,1',
	'AGCGAA'	=>	'2,1',
	'AGCGAC'	=>	'2,0',
	'AGCGAG'	=>	'2,1',
	'AGCGAT'	=>	'2,1',
	'AGCGCA'	=>	'2,1',
	'AGCGCC'	=>	'2,0',
	'AGCGCG'	=>	'2,1',
	'AGCGCT'	=>	'2,1',
	'AGCGGA'	=>	'1,1',
	'AGCGGC'	=>	'1,0',
	'AGCGGG'	=>	'1,1',
	'AGCGGT'	=>	'1,1',
	'AGCGTA'	=>	'2,1',
	'AGCGTC'	=>	'2,0',
	'AGCGTG'	=>	'2,1',
	'AGCGTT'	=>	'2,1',
	'AGCTAC'	=>	'2,0',
	'AGCTAT'	=>	'2,1',
	'AGCTCA'	=>	'2,1',
	'AGCTCC'	=>	'2,0',
	'AGCTCG'	=>	'2,1',
	'AGCTCT'	=>	'2,1',
	'AGCTGC'	=>	'1,0',
	'AGCTGG'	=>	'2,0',
	'AGCTGT'	=>	'1,1',
	'AGCTTA'	=>	'2,1',
	'AGCTTC'	=>	'2,0',
	'AGCTTG'	=>	'3,0',
	'AGCTTT'	=>	'2,1',
	'AGGAAA'	=>	'1,1',
	'AGGAAC'	=>	'2,0',
	'AGGAAG'	=>	'1,0',
	'AGGAAT'	=>	'2,0',
	'AGGACA'	=>	'1,1',
	'AGGACC'	=>	'1,1',
	'AGGACG'	=>	'1,0',
	'AGGACT'	=>	'1,1',
	'AGGAGA'	=>	'0,1',
	'AGGAGC'	=>	'1,0',
	'AGGAGG'	=>	'0,0',
	'AGGAGT'	=>	'1,0',
	'AGGATA'	=>	'1,1',
	'AGGATC'	=>	'2,0',
	'AGGATG'	=>	'1,0',
	'AGGATT'	=>	'2,0',
	'AGGCAA'	=>	'1,2',
	'AGGCAC'	=>	'1,2',
	'AGGCAG'	=>	'1,1',
	'AGGCAT'	=>	'1,2',
	'AGGCCA'	=>	'1,2',
	'AGGCCC'	=>	'1,2',
	'AGGCCG'	=>	'1,1',
	'AGGCCT'	=>	'1,2',
	'AGGCGA'	=>	'0,2',
	'AGGCGC'	=>	'0,2',
	'AGGCGG'	=>	'0,1',
	'AGGCGT'	=>	'0,2',
	'AGGCTA'	=>	'1,2',
	'AGGCTC'	=>	'1,2',
	'AGGCTG'	=>	'1,1',
	'AGGCTT'	=>	'1,2',
	'AGGGAA'	=>	'2,1',
	'AGGGAC'	=>	'2,1',
	'AGGGAG'	=>	'2,0',
	'AGGGAT'	=>	'2,1',
	'AGGGCA'	=>	'2,1',
	'AGGGCC'	=>	'2,1',
	'AGGGCG'	=>	'2,0',
	'AGGGCT'	=>	'2,1',
	'AGGGGA'	=>	'1,1',
	'AGGGGC'	=>	'1,1',
	'AGGGGG'	=>	'1,0',
	'AGGGGT'	=>	'1,1',
	'AGGGTA'	=>	'2,1',
	'AGGGTC'	=>	'2,1',
	'AGGGTG'	=>	'2,0',
	'AGGGTT'	=>	'2,1',
	'AGGTAC'	=>	'3,0',
	'AGGTAT'	=>	'3,0',
	'AGGTCA'	=>	'2,1',
	'AGGTCC'	=>	'2,1',
	'AGGTCG'	=>	'2,0',
	'AGGTCT'	=>	'2,1',
	'AGGTGC'	=>	'2,0',
	'AGGTGG'	=>	'1,0',
	'AGGTGT'	=>	'2,0',
	'AGGTTA'	=>	'2,1',
	'AGGTTC'	=>	'3,0',
	'AGGTTG'	=>	'2,0',
	'AGGTTT'	=>	'3,0',
	'AGTAAA'	=>	'2,0',
	'AGTAAC'	=>	'1,1',
	'AGTAAG'	=>	'2,0',
	'AGTAAT'	=>	'1,0',
	'AGTACA'	=>	'1,1',
	'AGTACC'	=>	'1,1',
	'AGTACG'	=>	'1,1',
	'AGTACT'	=>	'1,0',
	'AGTAGA'	=>	'1,0',
	'AGTAGC'	=>	'0,1',
	'AGTAGG'	=>	'1,0',
	'AGTAGT'	=>	'0,0',
	'AGTATA'	=>	'1,1',
	'AGTATC'	=>	'1,1',
	'AGTATG'	=>	'2,0',
	'AGTATT'	=>	'1,0',
	'AGTCAA'	=>	'2,1',
	'AGTCAC'	=>	'2,1',
	'AGTCAG'	=>	'2,1',
	'AGTCAT'	=>	'2,0',
	'AGTCCA'	=>	'2,1',
	'AGTCCC'	=>	'2,1',
	'AGTCCG'	=>	'2,1',
	'AGTCCT'	=>	'2,0',
	'AGTCGA'	=>	'1,1',
	'AGTCGC'	=>	'1,1',
	'AGTCGG'	=>	'1,1',
	'AGTCGT'	=>	'1,0',
	'AGTCTA'	=>	'2,1',
	'AGTCTC'	=>	'2,1',
	'AGTCTG'	=>	'2,1',
	'AGTCTT'	=>	'2,0',
	'AGTGAA'	=>	'2,1',
	'AGTGAC'	=>	'2,1',
	'AGTGAG'	=>	'2,1',
	'AGTGAT'	=>	'2,0',
	'AGTGCA'	=>	'2,1',
	'AGTGCC'	=>	'2,1',
	'AGTGCG'	=>	'2,1',
	'AGTGCT'	=>	'2,0',
	'AGTGGA'	=>	'1,1',
	'AGTGGC'	=>	'1,1',
	'AGTGGG'	=>	'1,1',
	'AGTGGT'	=>	'1,0',
	'AGTGTA'	=>	'2,1',
	'AGTGTC'	=>	'2,1',
	'AGTGTG'	=>	'2,1',
	'AGTGTT'	=>	'2,0',
	'AGTTAC'	=>	'2,1',
	'AGTTAT'	=>	'2,0',
	'AGTTCA'	=>	'2,1',
	'AGTTCC'	=>	'2,1',
	'AGTTCG'	=>	'2,1',
	'AGTTCT'	=>	'2,0',
	'AGTTGC'	=>	'1,1',
	'AGTTGG'	=>	'2,0',
	'AGTTGT'	=>	'1,0',
	'AGTTTA'	=>	'2,1',
	'AGTTTC'	=>	'2,1',
	'AGTTTG'	=>	'3,0',
	'AGTTTT'	=>	'2,0',
	'ATAAAA'	=>	'1,0',
	'ATAAAC'	=>	'1,1',
	'ATAAAG'	=>	'1,1',
	'ATAAAT'	=>	'1,1',
	'ATAACA'	=>	'1,0',
	'ATAACC'	=>	'1,1',
	'ATAACG'	=>	'1,1',
	'ATAACT'	=>	'1,1',
	'ATAAGA'	=>	'1,0',
	'ATAAGC'	=>	'1,1',
	'ATAAGG'	=>	'1,1',
	'ATAAGT'	=>	'1,1',
	'ATAATA'	=>	'0,0',
	'ATAATC'	=>	'0,1',
	'ATAATG'	=>	'1,0',
	'ATAATT'	=>	'0,1',
	'ATACAA'	=>	'2,0',
	'ATACAC'	=>	'2,1',
	'ATACAG'	=>	'2,1',
	'ATACAT'	=>	'2,1',
	'ATACCA'	=>	'2,0',
	'ATACCC'	=>	'2,1',
	'ATACCG'	=>	'2,1',
	'ATACCT'	=>	'2,1',
	'ATACGA'	=>	'1,1',
	'ATACGC'	=>	'1,2',
	'ATACGG'	=>	'1,2',
	'ATACGT'	=>	'1,2',
	'ATACTA'	=>	'1,0',
	'ATACTC'	=>	'1,1',
	'ATACTG'	=>	'1,1',
	'ATACTT'	=>	'1,1',
	'ATAGAA'	=>	'2,0',
	'ATAGAC'	=>	'2,1',
	'ATAGAG'	=>	'2,1',
	'ATAGAT'	=>	'2,1',
	'ATAGCA'	=>	'2,0',
	'ATAGCC'	=>	'2,1',
	'ATAGCG'	=>	'2,1',
	'ATAGCT'	=>	'2,1',
	'ATAGGA'	=>	'2,0',
	'ATAGGC'	=>	'2,1',
	'ATAGGG'	=>	'2,1',
	'ATAGGT'	=>	'2,1',
	'ATAGTA'	=>	'1,0',
	'ATAGTC'	=>	'1,1',
	'ATAGTG'	=>	'1,1',
	'ATAGTT'	=>	'1,1',
	'ATATAC'	=>	'2,1',
	'ATATAT'	=>	'2,1',
	'ATATCA'	=>	'2,0',
	'ATATCC'	=>	'2,1',
	'ATATCG'	=>	'2,1',
	'ATATCT'	=>	'2,1',
	'ATATGC'	=>	'2,1',
	'ATATGG'	=>	'2,1',
	'ATATGT'	=>	'2,1',
	'ATATTA'	=>	'1,0',
	'ATATTC'	=>	'1,1',
	'ATATTG'	=>	'1,1',
	'ATATTT'	=>	'1,1',
	'ATCAAA'	=>	'1,1',
	'ATCAAC'	=>	'1,0',
	'ATCAAG'	=>	'2,0',
	'ATCAAT'	=>	'1,1',
	'ATCACA'	=>	'1,1',
	'ATCACC'	=>	'1,0',
	'ATCACG'	=>	'1,1',
	'ATCACT'	=>	'1,1',
	'ATCAGA'	=>	'1,1',
	'ATCAGC'	=>	'1,0',
	'ATCAGG'	=>	'2,0',
	'ATCAGT'	=>	'1,1',
	'ATCATA'	=>	'0,1',
	'ATCATC'	=>	'0,0',
	'ATCATG'	=>	'1,0',
	'ATCATT'	=>	'0,1',
	'ATCCAA'	=>	'2,1',
	'ATCCAC'	=>	'2,0',
	'ATCCAG'	=>	'2,1',
	'ATCCAT'	=>	'2,1',
	'ATCCCA'	=>	'2,1',
	'ATCCCC'	=>	'2,0',
	'ATCCCG'	=>	'2,1',
	'ATCCCT'	=>	'2,1',
	'ATCCGA'	=>	'1,2',
	'ATCCGC'	=>	'2,0',
	'ATCCGG'	=>	'2,1',
	'ATCCGT'	=>	'2,1',
	'ATCCTA'	=>	'1,1',
	'ATCCTC'	=>	'1,0',
	'ATCCTG'	=>	'1,1',
	'ATCCTT'	=>	'1,1',
	'ATCGAA'	=>	'2,1',
	'ATCGAC'	=>	'2,0',
	'ATCGAG'	=>	'2,1',
	'ATCGAT'	=>	'2,1',
	'ATCGCA'	=>	'2,1',
	'ATCGCC'	=>	'2,0',
	'ATCGCG'	=>	'2,1',
	'ATCGCT'	=>	'2,1',
	'ATCGGA'	=>	'2,1',
	'ATCGGC'	=>	'2,0',
	'ATCGGG'	=>	'2,1',
	'ATCGGT'	=>	'2,1',
	'ATCGTA'	=>	'1,1',
	'ATCGTC'	=>	'1,0',
	'ATCGTG'	=>	'1,1',
	'ATCGTT'	=>	'1,1',
	'ATCTAC'	=>	'2,0',
	'ATCTAT'	=>	'2,1',
	'ATCTCA'	=>	'2,1',
	'ATCTCC'	=>	'2,0',
	'ATCTCG'	=>	'2,1',
	'ATCTCT'	=>	'2,1',
	'ATCTGC'	=>	'2,0',
	'ATCTGG'	=>	'3,0',
	'ATCTGT'	=>	'2,1',
	'ATCTTA'	=>	'1,1',
	'ATCTTC'	=>	'1,0',
	'ATCTTG'	=>	'2,0',
	'ATCTTT'	=>	'1,1',
	'ATGAAA'	=>	'1,1',
	'ATGAAC'	=>	'2,0',
	'ATGAAG'	=>	'1,0',
	'ATGAAT'	=>	'2,0',
	'ATGACA'	=>	'1,1',
	'ATGACC'	=>	'1,1',
	'ATGACG'	=>	'1,0',
	'ATGACT'	=>	'1,1',
	'ATGAGA'	=>	'1,1',
	'ATGAGC'	=>	'2,0',
	'ATGAGG'	=>	'1,0',
	'ATGAGT'	=>	'2,0',
	'ATGATA'	=>	'1,0',
	'ATGATC'	=>	'1,0',
	'ATGATG'	=>	'0,0',
	'ATGATT'	=>	'1,0',
	'ATGCAA'	=>	'2,1',
	'ATGCAC'	=>	'2,1',
	'ATGCAG'	=>	'2,0',
	'ATGCAT'	=>	'2,1',
	'ATGCCA'	=>	'2,1',
	'ATGCCC'	=>	'2,1',
	'ATGCCG'	=>	'2,0',
	'ATGCCT'	=>	'2,1',
	'ATGCGA'	=>	'1,2',
	'ATGCGC'	=>	'1,2',
	'ATGCGG'	=>	'1,1',
	'ATGCGT'	=>	'1,2',
	'ATGCTA'	=>	'1,1',
	'ATGCTC'	=>	'1,1',
	'ATGCTG'	=>	'1,0',
	'ATGCTT'	=>	'1,1',
	'ATGGAA'	=>	'2,1',
	'ATGGAC'	=>	'2,1',
	'ATGGAG'	=>	'2,0',
	'ATGGAT'	=>	'2,1',
	'ATGGCA'	=>	'2,1',
	'ATGGCC'	=>	'2,1',
	'ATGGCG'	=>	'2,0',
	'ATGGCT'	=>	'2,1',
	'ATGGGA'	=>	'2,1',
	'ATGGGC'	=>	'2,1',
	'ATGGGG'	=>	'2,0',
	'ATGGGT'	=>	'2,1',
	'ATGGTA'	=>	'1,1',
	'ATGGTC'	=>	'1,1',
	'ATGGTG'	=>	'1,0',
	'ATGGTT'	=>	'1,1',
	'ATGTAC'	=>	'3,0',
	'ATGTAT'	=>	'3,0',
	'ATGTCA'	=>	'2,1',
	'ATGTCC'	=>	'2,1',
	'ATGTCG'	=>	'2,0',
	'ATGTCT'	=>	'2,1',
	'ATGTGC'	=>	'3,0',
	'ATGTGG'	=>	'2,0',
	'ATGTGT'	=>	'3,0',
	'ATGTTA'	=>	'1,1',
	'ATGTTC'	=>	'2,0',
	'ATGTTG'	=>	'1,0',
	'ATGTTT'	=>	'2,0',
	'ATTAAA'	=>	'1,1',
	'ATTAAC'	=>	'1,1',
	'ATTAAG'	=>	'2,0',
	'ATTAAT'	=>	'1,0',
	'ATTACA'	=>	'1,1',
	'ATTACC'	=>	'1,1',
	'ATTACG'	=>	'1,1',
	'ATTACT'	=>	'1,0',
	'ATTAGA'	=>	'1,1',
	'ATTAGC'	=>	'1,1',
	'ATTAGG'	=>	'2,0',
	'ATTAGT'	=>	'1,0',
	'ATTATA'	=>	'0,1',
	'ATTATC'	=>	'0,1',
	'ATTATG'	=>	'1,0',
	'ATTATT'	=>	'0,0',
	'ATTCAA'	=>	'2,1',
	'ATTCAC'	=>	'2,1',
	'ATTCAG'	=>	'2,1',
	'ATTCAT'	=>	'2,0',
	'ATTCCA'	=>	'2,1',
	'ATTCCC'	=>	'2,1',
	'ATTCCG'	=>	'2,1',
	'ATTCCT'	=>	'2,0',
	'ATTCGA'	=>	'1,2',
	'ATTCGC'	=>	'2,1',
	'ATTCGG'	=>	'2,1',
	'ATTCGT'	=>	'2,0',
	'ATTCTA'	=>	'1,1',
	'ATTCTC'	=>	'1,1',
	'ATTCTG'	=>	'1,1',
	'ATTCTT'	=>	'1,0',
	'ATTGAA'	=>	'2,1',
	'ATTGAC'	=>	'2,1',
	'ATTGAG'	=>	'2,1',
	'ATTGAT'	=>	'2,0',
	'ATTGCA'	=>	'2,1',
	'ATTGCC'	=>	'2,1',
	'ATTGCG'	=>	'2,1',
	'ATTGCT'	=>	'2,0',
	'ATTGGA'	=>	'2,1',
	'ATTGGC'	=>	'2,1',
	'ATTGGG'	=>	'2,1',
	'ATTGGT'	=>	'2,0',
	'ATTGTA'	=>	'1,1',
	'ATTGTC'	=>	'1,1',
	'ATTGTG'	=>	'1,1',
	'ATTGTT'	=>	'1,0',
	'ATTTAC'	=>	'2,1',
	'ATTTAT'	=>	'2,0',
	'ATTTCA'	=>	'2,1',
	'ATTTCC'	=>	'2,1',
	'ATTTCG'	=>	'2,1',
	'ATTTCT'	=>	'2,0',
	'ATTTGC'	=>	'2,1',
	'ATTTGG'	=>	'3,0',
	'ATTTGT'	=>	'2,0',
	'ATTTTA'	=>	'1,1',
	'ATTTTC'	=>	'1,1',
	'ATTTTG'	=>	'2,0',
	'ATTTTT'	=>	'1,0',
	'CAAAAA'	=>	'1,0',
	'CAAAAC'	=>	'2,0',
	'CAAAAG'	=>	'1,1',
	'CAAAAT'	=>	'2,0',
	'CAAACA'	=>	'2,0',
	'CAAACC'	=>	'2,1',
	'CAAACG'	=>	'2,1',
	'CAAACT'	=>	'2,1',
	'CAAAGA'	=>	'1,1',
	'CAAAGC'	=>	'2,1',
	'CAAAGG'	=>	'1,2',
	'CAAAGT'	=>	'2,1',
	'CAAATA'	=>	'2,0',
	'CAAATC'	=>	'2,1',
	'CAAATG'	=>	'2,1',
	'CAAATT'	=>	'2,1',
	'CAACAA'	=>	'0,0',
	'CAACAC'	=>	'1,0',
	'CAACAG'	=>	'0,1',
	'CAACAT'	=>	'1,0',
	'CAACCA'	=>	'1,0',
	'CAACCC'	=>	'1,1',
	'CAACCG'	=>	'1,1',
	'CAACCT'	=>	'1,1',
	'CAACGA'	=>	'1,0',
	'CAACGC'	=>	'1,1',
	'CAACGG'	=>	'1,1',
	'CAACGT'	=>	'1,1',
	'CAACTA'	=>	'1,0',
	'CAACTC'	=>	'1,1',
	'CAACTG'	=>	'1,1',
	'CAACTT'	=>	'1,1',
	'CAAGAA'	=>	'1,0',
	'CAAGAC'	=>	'2,0',
	'CAAGAG'	=>	'1,1',
	'CAAGAT'	=>	'2,0',
	'CAAGCA'	=>	'2,0',
	'CAAGCC'	=>	'2,1',
	'CAAGCG'	=>	'2,1',
	'CAAGCT'	=>	'2,1',
	'CAAGGA'	=>	'2,0',
	'CAAGGC'	=>	'2,1',
	'CAAGGG'	=>	'2,1',
	'CAAGGT'	=>	'2,1',
	'CAAGTA'	=>	'2,0',
	'CAAGTC'	=>	'2,1',
	'CAAGTG'	=>	'2,1',
	'CAAGTT'	=>	'2,1',
	'CAATAC'	=>	'2,0',
	'CAATAT'	=>	'2,0',
	'CAATCA'	=>	'2,0',
	'CAATCC'	=>	'2,1',
	'CAATCG'	=>	'2,1',
	'CAATCT'	=>	'2,1',
	'CAATGC'	=>	'2,1',
	'CAATGG'	=>	'2,1',
	'CAATGT'	=>	'2,1',
	'CAATTA'	=>	'1,1',
	'CAATTC'	=>	'2,1',
	'CAATTG'	=>	'1,2',
	'CAATTT'	=>	'2,1',
	'CACAAA'	=>	'2,0',
	'CACAAC'	=>	'1,0',
	'CACAAG'	=>	'2,0',
	'CACAAT'	=>	'1,1',
	'CACACA'	=>	'2,1',
	'CACACC'	=>	'2,0',
	'CACACG'	=>	'2,1',
	'CACACT'	=>	'2,1',
	'CACAGA'	=>	'1,2',
	'CACAGC'	=>	'2,0',
	'CACAGG'	=>	'1,2',
	'CACAGT'	=>	'2,1',
	'CACATA'	=>	'2,1',
	'CACATC'	=>	'2,0',
	'CACATG'	=>	'2,1',
	'CACATT'	=>	'2,1',
	'CACCAA'	=>	'1,0',
	'CACCAC'	=>	'0,0',
	'CACCAG'	=>	'1,0',
	'CACCAT'	=>	'0,1',
	'CACCCA'	=>	'1,1',
	'CACCCC'	=>	'1,0',
	'CACCCG'	=>	'1,1',
	'CACCCT'	=>	'1,1',
	'CACCGA'	=>	'1,1',
	'CACCGC'	=>	'1,0',
	'CACCGG'	=>	'1,1',
	'CACCGT'	=>	'1,1',
	'CACCTA'	=>	'1,1',
	'CACCTC'	=>	'1,0',
	'CACCTG'	=>	'1,1',
	'CACCTT'	=>	'1,1',
	'CACGAA'	=>	'2,0',
	'CACGAC'	=>	'1,0',
	'CACGAG'	=>	'2,0',
	'CACGAT'	=>	'1,1',
	'CACGCA'	=>	'2,1',
	'CACGCC'	=>	'2,0',
	'CACGCG'	=>	'2,1',
	'CACGCT'	=>	'2,1',
	'CACGGA'	=>	'2,1',
	'CACGGC'	=>	'2,0',
	'CACGGG'	=>	'2,1',
	'CACGGT'	=>	'2,1',
	'CACGTA'	=>	'2,1',
	'CACGTC'	=>	'2,0',
	'CACGTG'	=>	'2,1',
	'CACGTT'	=>	'2,1',
	'CACTAC'	=>	'1,0',
	'CACTAT'	=>	'1,1',
	'CACTCA'	=>	'2,1',
	'CACTCC'	=>	'2,0',
	'CACTCG'	=>	'2,1',
	'CACTCT'	=>	'2,1',
	'CACTGC'	=>	'2,0',
	'CACTGG'	=>	'2,1',
	'CACTGT'	=>	'2,1',
	'CACTTA'	=>	'1,2',
	'CACTTC'	=>	'2,0',
	'CACTTG'	=>	'1,2',
	'CACTTT'	=>	'2,1',
	'CAGAAA'	=>	'1,1',
	'CAGAAC'	=>	'2,0',
	'CAGAAG'	=>	'1,0',
	'CAGAAT'	=>	'2,0',
	'CAGACA'	=>	'2,1',
	'CAGACC'	=>	'2,1',
	'CAGACG'	=>	'2,0',
	'CAGACT'	=>	'2,1',
	'CAGAGA'	=>	'1,2',
	'CAGAGC'	=>	'2,1',
	'CAGAGG'	=>	'1,1',
	'CAGAGT'	=>	'2,1',
	'CAGATA'	=>	'2,1',
	'CAGATC'	=>	'2,1',
	'CAGATG'	=>	'2,0',
	'CAGATT'	=>	'2,1',
	'CAGCAA'	=>	'0,1',
	'CAGCAC'	=>	'1,0',
	'CAGCAG'	=>	'0,0',
	'CAGCAT'	=>	'1,0',
	'CAGCCA'	=>	'1,1',
	'CAGCCC'	=>	'1,1',
	'CAGCCG'	=>	'1,0',
	'CAGCCT'	=>	'1,1',
	'CAGCGA'	=>	'1,1',
	'CAGCGC'	=>	'1,1',
	'CAGCGG'	=>	'1,0',
	'CAGCGT'	=>	'1,1',
	'CAGCTA'	=>	'1,1',
	'CAGCTC'	=>	'1,1',
	'CAGCTG'	=>	'1,0',
	'CAGCTT'	=>	'1,1',
	'CAGGAA'	=>	'1,1',
	'CAGGAC'	=>	'2,0',
	'CAGGAG'	=>	'1,0',
	'CAGGAT'	=>	'2,0',
	'CAGGCA'	=>	'2,1',
	'CAGGCC'	=>	'2,1',
	'CAGGCG'	=>	'2,0',
	'CAGGCT'	=>	'2,1',
	'CAGGGA'	=>	'2,1',
	'CAGGGC'	=>	'2,1',
	'CAGGGG'	=>	'2,0',
	'CAGGGT'	=>	'2,1',
	'CAGGTA'	=>	'2,1',
	'CAGGTC'	=>	'2,1',
	'CAGGTG'	=>	'2,0',
	'CAGGTT'	=>	'2,1',
	'CAGTAC'	=>	'2,0',
	'CAGTAT'	=>	'2,0',
	'CAGTCA'	=>	'2,1',
	'CAGTCC'	=>	'2,1',
	'CAGTCG'	=>	'2,0',
	'CAGTCT'	=>	'2,1',
	'CAGTGC'	=>	'2,1',
	'CAGTGG'	=>	'2,0',
	'CAGTGT'	=>	'2,1',
	'CAGTTA'	=>	'1,2',
	'CAGTTC'	=>	'2,1',
	'CAGTTG'	=>	'1,1',
	'CAGTTT'	=>	'2,1',
	'CATAAA'	=>	'2,0',
	'CATAAC'	=>	'1,1',
	'CATAAG'	=>	'2,0',
	'CATAAT'	=>	'1,0',
	'CATACA'	=>	'2,1',
	'CATACC'	=>	'2,1',
	'CATACG'	=>	'2,1',
	'CATACT'	=>	'2,0',
	'CATAGA'	=>	'1,2',
	'CATAGC'	=>	'2,1',
	'CATAGG'	=>	'1,2',
	'CATAGT'	=>	'2,0',
	'CATATA'	=>	'2,1',
	'CATATC'	=>	'2,1',
	'CATATG'	=>	'2,1',
	'CATATT'	=>	'2,0',
	'CATCAA'	=>	'1,0',
	'CATCAC'	=>	'0,1',
	'CATCAG'	=>	'1,0',
	'CATCAT'	=>	'0,0',
	'CATCCA'	=>	'1,1',
	'CATCCC'	=>	'1,1',
	'CATCCG'	=>	'1,1',
	'CATCCT'	=>	'1,0',
	'CATCGA'	=>	'1,1',
	'CATCGC'	=>	'1,1',
	'CATCGG'	=>	'1,1',
	'CATCGT'	=>	'1,0',
	'CATCTA'	=>	'1,1',
	'CATCTC'	=>	'1,1',
	'CATCTG'	=>	'1,1',
	'CATCTT'	=>	'1,0',
	'CATGAA'	=>	'2,0',
	'CATGAC'	=>	'1,1',
	'CATGAG'	=>	'2,0',
	'CATGAT'	=>	'1,0',
	'CATGCA'	=>	'2,1',
	'CATGCC'	=>	'2,1',
	'CATGCG'	=>	'2,1',
	'CATGCT'	=>	'2,0',
	'CATGGA'	=>	'2,1',
	'CATGGC'	=>	'2,1',
	'CATGGG'	=>	'2,1',
	'CATGGT'	=>	'2,0',
	'CATGTA'	=>	'2,1',
	'CATGTC'	=>	'2,1',
	'CATGTG'	=>	'2,1',
	'CATGTT'	=>	'2,0',
	'CATTAC'	=>	'1,1',
	'CATTAT'	=>	'1,0',
	'CATTCA'	=>	'2,1',
	'CATTCC'	=>	'2,1',
	'CATTCG'	=>	'2,1',
	'CATTCT'	=>	'2,0',
	'CATTGC'	=>	'2,1',
	'CATTGG'	=>	'2,1',
	'CATTGT'	=>	'2,0',
	'CATTTA'	=>	'1,2',
	'CATTTC'	=>	'2,1',
	'CATTTG'	=>	'1,2',
	'CATTTT'	=>	'2,0',
	'CCAAAA'	=>	'2,0',
	'CCAAAC'	=>	'2,1',
	'CCAAAG'	=>	'2,1',
	'CCAAAT'	=>	'2,1',
	'CCAACA'	=>	'1,0',
	'CCAACC'	=>	'1,1',
	'CCAACG'	=>	'1,1',
	'CCAACT'	=>	'1,1',
	'CCAAGA'	=>	'1,1',
	'CCAAGC'	=>	'2,1',
	'CCAAGG'	=>	'1,2',
	'CCAAGT'	=>	'2,1',
	'CCAATA'	=>	'2,0',
	'CCAATC'	=>	'2,1',
	'CCAATG'	=>	'2,1',
	'CCAATT'	=>	'2,1',
	'CCACAA'	=>	'1,0',
	'CCACAC'	=>	'1,1',
	'CCACAG'	=>	'1,1',
	'CCACAT'	=>	'1,1',
	'CCACCA'	=>	'0,0',
	'CCACCC'	=>	'0,1',
	'CCACCG'	=>	'0,1',
	'CCACCT'	=>	'0,1',
	'CCACGA'	=>	'1,0',
	'CCACGC'	=>	'1,1',
	'CCACGG'	=>	'1,1',
	'CCACGT'	=>	'1,1',
	'CCACTA'	=>	'1,0',
	'CCACTC'	=>	'1,1',
	'CCACTG'	=>	'1,1',
	'CCACTT'	=>	'1,1',
	'CCAGAA'	=>	'2,0',
	'CCAGAC'	=>	'2,1',
	'CCAGAG'	=>	'2,1',
	'CCAGAT'	=>	'2,1',
	'CCAGCA'	=>	'1,0',
	'CCAGCC'	=>	'1,1',
	'CCAGCG'	=>	'1,1',
	'CCAGCT'	=>	'1,1',
	'CCAGGA'	=>	'2,0',
	'CCAGGC'	=>	'2,1',
	'CCAGGG'	=>	'2,1',
	'CCAGGT'	=>	'2,1',
	'CCAGTA'	=>	'2,0',
	'CCAGTC'	=>	'2,1',
	'CCAGTG'	=>	'2,1',
	'CCAGTT'	=>	'2,1',
	'CCATAC'	=>	'2,1',
	'CCATAT'	=>	'2,1',
	'CCATCA'	=>	'1,0',
	'CCATCC'	=>	'1,1',
	'CCATCG'	=>	'1,1',
	'CCATCT'	=>	'1,1',
	'CCATGC'	=>	'2,1',
	'CCATGG'	=>	'2,1',
	'CCATGT'	=>	'2,1',
	'CCATTA'	=>	'1,1',
	'CCATTC'	=>	'2,1',
	'CCATTG'	=>	'1,2',
	'CCATTT'	=>	'2,1',
	'CCCAAA'	=>	'2,1',
	'CCCAAC'	=>	'2,0',
	'CCCAAG'	=>	'2,1',
	'CCCAAT'	=>	'2,1',
	'CCCACA'	=>	'1,1',
	'CCCACC'	=>	'1,0',
	'CCCACG'	=>	'1,1',
	'CCCACT'	=>	'1,1',
	'CCCAGA'	=>	'1,2',
	'CCCAGC'	=>	'2,0',
	'CCCAGG'	=>	'1,2',
	'CCCAGT'	=>	'2,1',
	'CCCATA'	=>	'2,1',
	'CCCATC'	=>	'2,0',
	'CCCATG'	=>	'2,1',
	'CCCATT'	=>	'2,1',
	'CCCCAA'	=>	'1,1',
	'CCCCAC'	=>	'1,0',
	'CCCCAG'	=>	'1,1',
	'CCCCAT'	=>	'1,1',
	'CCCCCA'	=>	'0,1',
	'CCCCCC'	=>	'0,0',
	'CCCCCG'	=>	'0,1',
	'CCCCCT'	=>	'0,1',
	'CCCCGA'	=>	'1,1',
	'CCCCGC'	=>	'1,0',
	'CCCCGG'	=>	'1,1',
	'CCCCGT'	=>	'1,1',
	'CCCCTA'	=>	'1,1',
	'CCCCTC'	=>	'1,0',
	'CCCCTG'	=>	'1,1',
	'CCCCTT'	=>	'1,1',
	'CCCGAA'	=>	'2,1',
	'CCCGAC'	=>	'2,0',
	'CCCGAG'	=>	'2,1',
	'CCCGAT'	=>	'2,1',
	'CCCGCA'	=>	'1,1',
	'CCCGCC'	=>	'1,0',
	'CCCGCG'	=>	'1,1',
	'CCCGCT'	=>	'1,1',
	'CCCGGA'	=>	'2,1',
	'CCCGGC'	=>	'2,0',
	'CCCGGG'	=>	'2,1',
	'CCCGGT'	=>	'2,1',
	'CCCGTA'	=>	'2,1',
	'CCCGTC'	=>	'2,0',
	'CCCGTG'	=>	'2,1',
	'CCCGTT'	=>	'2,1',
	'CCCTAC'	=>	'2,0',
	'CCCTAT'	=>	'2,1',
	'CCCTCA'	=>	'1,1',
	'CCCTCC'	=>	'1,0',
	'CCCTCG'	=>	'1,1',
	'CCCTCT'	=>	'1,1',
	'CCCTGC'	=>	'2,0',
	'CCCTGG'	=>	'2,1',
	'CCCTGT'	=>	'2,1',
	'CCCTTA'	=>	'1,2',
	'CCCTTC'	=>	'2,0',
	'CCCTTG'	=>	'1,2',
	'CCCTTT'	=>	'2,1',
	'CCGAAA'	=>	'2,1',
	'CCGAAC'	=>	'2,1',
	'CCGAAG'	=>	'2,0',
	'CCGAAT'	=>	'2,1',
	'CCGACA'	=>	'1,1',
	'CCGACC'	=>	'1,1',
	'CCGACG'	=>	'1,0',
	'CCGACT'	=>	'1,1',
	'CCGAGA'	=>	'1,2',
	'CCGAGC'	=>	'2,1',
	'CCGAGG'	=>	'1,1',
	'CCGAGT'	=>	'2,1',
	'CCGATA'	=>	'2,1',
	'CCGATC'	=>	'2,1',
	'CCGATG'	=>	'2,0',
	'CCGATT'	=>	'2,1',
	'CCGCAA'	=>	'1,1',
	'CCGCAC'	=>	'1,1',
	'CCGCAG'	=>	'1,0',
	'CCGCAT'	=>	'1,1',
	'CCGCCA'	=>	'0,1',
	'CCGCCC'	=>	'0,1',
	'CCGCCG'	=>	'0,0',
	'CCGCCT'	=>	'0,1',
	'CCGCGA'	=>	'1,1',
	'CCGCGC'	=>	'1,1',
	'CCGCGG'	=>	'1,0',
	'CCGCGT'	=>	'1,1',
	'CCGCTA'	=>	'1,1',
	'CCGCTC'	=>	'1,1',
	'CCGCTG'	=>	'1,0',
	'CCGCTT'	=>	'1,1',
	'CCGGAA'	=>	'2,1',
	'CCGGAC'	=>	'2,1',
	'CCGGAG'	=>	'2,0',
	'CCGGAT'	=>	'2,1',
	'CCGGCA'	=>	'1,1',
	'CCGGCC'	=>	'1,1',
	'CCGGCG'	=>	'1,0',
	'CCGGCT'	=>	'1,1',
	'CCGGGA'	=>	'2,1',
	'CCGGGC'	=>	'2,1',
	'CCGGGG'	=>	'2,0',
	'CCGGGT'	=>	'2,1',
	'CCGGTA'	=>	'2,1',
	'CCGGTC'	=>	'2,1',
	'CCGGTG'	=>	'2,0',
	'CCGGTT'	=>	'2,1',
	'CCGTAC'	=>	'2,1',
	'CCGTAT'	=>	'2,1',
	'CCGTCA'	=>	'1,1',
	'CCGTCC'	=>	'1,1',
	'CCGTCG'	=>	'1,0',
	'CCGTCT'	=>	'1,1',
	'CCGTGC'	=>	'2,1',
	'CCGTGG'	=>	'2,0',
	'CCGTGT'	=>	'2,1',
	'CCGTTA'	=>	'1,2',
	'CCGTTC'	=>	'2,1',
	'CCGTTG'	=>	'1,1',
	'CCGTTT'	=>	'2,1',
	'CCTAAA'	=>	'2,1',
	'CCTAAC'	=>	'2,1',
	'CCTAAG'	=>	'2,1',
	'CCTAAT'	=>	'2,0',
	'CCTACA'	=>	'1,1',
	'CCTACC'	=>	'1,1',
	'CCTACG'	=>	'1,1',
	'CCTACT'	=>	'1,0',
	'CCTAGA'	=>	'1,2',
	'CCTAGC'	=>	'2,1',
	'CCTAGG'	=>	'1,2',
	'CCTAGT'	=>	'2,0',
	'CCTATA'	=>	'2,1',
	'CCTATC'	=>	'2,1',
	'CCTATG'	=>	'2,1',
	'CCTATT'	=>	'2,0',
	'CCTCAA'	=>	'1,1',
	'CCTCAC'	=>	'1,1',
	'CCTCAG'	=>	'1,1',
	'CCTCAT'	=>	'1,0',
	'CCTCCA'	=>	'0,1',
	'CCTCCC'	=>	'0,1',
	'CCTCCG'	=>	'0,1',
	'CCTCCT'	=>	'0,0',
	'CCTCGA'	=>	'1,1',
	'CCTCGC'	=>	'1,1',
	'CCTCGG'	=>	'1,1',
	'CCTCGT'	=>	'1,0',
	'CCTCTA'	=>	'1,1',
	'CCTCTC'	=>	'1,1',
	'CCTCTG'	=>	'1,1',
	'CCTCTT'	=>	'1,0',
	'CCTGAA'	=>	'2,1',
	'CCTGAC'	=>	'2,1',
	'CCTGAG'	=>	'2,1',
	'CCTGAT'	=>	'2,0',
	'CCTGCA'	=>	'1,1',
	'CCTGCC'	=>	'1,1',
	'CCTGCG'	=>	'1,1',
	'CCTGCT'	=>	'1,0',
	'CCTGGA'	=>	'2,1',
	'CCTGGC'	=>	'2,1',
	'CCTGGG'	=>	'2,1',
	'CCTGGT'	=>	'2,0',
	'CCTGTA'	=>	'2,1',
	'CCTGTC'	=>	'2,1',
	'CCTGTG'	=>	'2,1',
	'CCTGTT'	=>	'2,0',
	'CCTTAC'	=>	'2,1',
	'CCTTAT'	=>	'2,0',
	'CCTTCA'	=>	'1,1',
	'CCTTCC'	=>	'1,1',
	'CCTTCG'	=>	'1,1',
	'CCTTCT'	=>	'1,0',
	'CCTTGC'	=>	'2,1',
	'CCTTGG'	=>	'2,1',
	'CCTTGT'	=>	'2,0',
	'CCTTTA'	=>	'1,2',
	'CCTTTC'	=>	'2,1',
	'CCTTTG'	=>	'1,2',
	'CCTTTT'	=>	'2,0',
	'CGAAAA'	=>	'1,1',
	'CGAAAC'	=>	'2,1',
	'CGAAAG'	=>	'1,2',
	'CGAAAT'	=>	'2,1',
	'CGAACA'	=>	'1,1',
	'CGAACC'	=>	'1,2',
	'CGAACG'	=>	'1,2',
	'CGAACT'	=>	'1,2',
	'CGAAGA'	=>	'0,1',
	'CGAAGC'	=>	'1,1',
	'CGAAGG'	=>	'0,2',
	'CGAAGT'	=>	'1,1',
	'CGAATA'	=>	'1,1',
	'CGAATC'	=>	'1,2',
	'CGAATG'	=>	'1,2',
	'CGAATT'	=>	'1,2',
	'CGACAA'	=>	'1,0',
	'CGACAC'	=>	'1,1',
	'CGACAG'	=>	'1,1',
	'CGACAT'	=>	'1,1',
	'CGACCA'	=>	'1,0',
	'CGACCC'	=>	'1,1',
	'CGACCG'	=>	'1,1',
	'CGACCT'	=>	'1,1',
	'CGACGA'	=>	'0,0',
	'CGACGC'	=>	'0,1',
	'CGACGG'	=>	'0,1',
	'CGACGT'	=>	'0,1',
	'CGACTA'	=>	'1,0',
	'CGACTC'	=>	'1,1',
	'CGACTG'	=>	'1,1',
	'CGACTT'	=>	'1,1',
	'CGAGAA'	=>	'2,0',
	'CGAGAC'	=>	'2,1',
	'CGAGAG'	=>	'2,1',
	'CGAGAT'	=>	'2,1',
	'CGAGCA'	=>	'2,0',
	'CGAGCC'	=>	'2,1',
	'CGAGCG'	=>	'2,1',
	'CGAGCT'	=>	'2,1',
	'CGAGGA'	=>	'1,0',
	'CGAGGC'	=>	'1,1',
	'CGAGGG'	=>	'1,1',
	'CGAGGT'	=>	'1,1',
	'CGAGTA'	=>	'2,0',
	'CGAGTC'	=>	'2,1',
	'CGAGTG'	=>	'2,1',
	'CGAGTT'	=>	'2,1',
	'CGATAC'	=>	'2,1',
	'CGATAT'	=>	'2,1',
	'CGATCA'	=>	'2,0',
	'CGATCC'	=>	'2,1',
	'CGATCG'	=>	'2,1',
	'CGATCT'	=>	'2,1',
	'CGATGC'	=>	'1,1',
	'CGATGG'	=>	'1,1',
	'CGATGT'	=>	'1,1',
	'CGATTA'	=>	'1,1',
	'CGATTC'	=>	'2,1',
	'CGATTG'	=>	'1,2',
	'CGATTT'	=>	'2,1',
	'CGCAAA'	=>	'1,2',
	'CGCAAC'	=>	'2,0',
	'CGCAAG'	=>	'1,2',
	'CGCAAT'	=>	'2,1',
	'CGCACA'	=>	'1,2',
	'CGCACC'	=>	'2,0',
	'CGCACG'	=>	'1,2',
	'CGCACT'	=>	'2,1',
	'CGCAGA'	=>	'0,2',
	'CGCAGC'	=>	'1,0',
	'CGCAGG'	=>	'0,2',
	'CGCAGT'	=>	'1,1',
	'CGCATA'	=>	'1,2',
	'CGCATC'	=>	'2,0',
	'CGCATG'	=>	'1,2',
	'CGCATT'	=>	'2,1',
	'CGCCAA'	=>	'1,1',
	'CGCCAC'	=>	'1,0',
	'CGCCAG'	=>	'1,1',
	'CGCCAT'	=>	'1,1',
	'CGCCCA'	=>	'1,1',
	'CGCCCC'	=>	'1,0',
	'CGCCCG'	=>	'1,1',
	'CGCCCT'	=>	'1,1',
	'CGCCGA'	=>	'0,1',
	'CGCCGC'	=>	'0,0',
	'CGCCGG'	=>	'0,1',
	'CGCCGT'	=>	'0,1',
	'CGCCTA'	=>	'1,1',
	'CGCCTC'	=>	'1,0',
	'CGCCTG'	=>	'1,1',
	'CGCCTT'	=>	'1,1',
	'CGCGAA'	=>	'2,1',
	'CGCGAC'	=>	'2,0',
	'CGCGAG'	=>	'2,1',
	'CGCGAT'	=>	'2,1',
	'CGCGCA'	=>	'2,1',
	'CGCGCC'	=>	'2,0',
	'CGCGCG'	=>	'2,1',
	'CGCGCT'	=>	'2,1',
	'CGCGGA'	=>	'1,1',
	'CGCGGC'	=>	'1,0',
	'CGCGGG'	=>	'1,1',
	'CGCGGT'	=>	'1,1',
	'CGCGTA'	=>	'2,1',
	'CGCGTC'	=>	'2,0',
	'CGCGTG'	=>	'2,1',
	'CGCGTT'	=>	'2,1',
	'CGCTAC'	=>	'2,0',
	'CGCTAT'	=>	'2,1',
	'CGCTCA'	=>	'2,1',
	'CGCTCC'	=>	'2,0',
	'CGCTCG'	=>	'2,1',
	'CGCTCT'	=>	'2,1',
	'CGCTGC'	=>	'1,0',
	'CGCTGG'	=>	'1,1',
	'CGCTGT'	=>	'1,1',
	'CGCTTA'	=>	'1,2',
	'CGCTTC'	=>	'2,0',
	'CGCTTG'	=>	'1,2',
	'CGCTTT'	=>	'2,1',
	'CGGAAA'	=>	'1,2',
	'CGGAAC'	=>	'2,1',
	'CGGAAG'	=>	'1,1',
	'CGGAAT'	=>	'2,1',
	'CGGACA'	=>	'1,2',
	'CGGACC'	=>	'1,2',
	'CGGACG'	=>	'1,1',
	'CGGACT'	=>	'1,2',
	'CGGAGA'	=>	'0,2',
	'CGGAGC'	=>	'1,1',
	'CGGAGG'	=>	'0,1',
	'CGGAGT'	=>	'1,1',
	'CGGATA'	=>	'1,2',
	'CGGATC'	=>	'2,1',
	'CGGATG'	=>	'1,1',
	'CGGATT'	=>	'2,1',
	'CGGCAA'	=>	'1,1',
	'CGGCAC'	=>	'1,1',
	'CGGCAG'	=>	'1,0',
	'CGGCAT'	=>	'1,1',
	'CGGCCA'	=>	'1,1',
	'CGGCCC'	=>	'1,1',
	'CGGCCG'	=>	'1,0',
	'CGGCCT'	=>	'1,1',
	'CGGCGA'	=>	'0,1',
	'CGGCGC'	=>	'0,1',
	'CGGCGG'	=>	'0,0',
	'CGGCGT'	=>	'0,1',
	'CGGCTA'	=>	'1,1',
	'CGGCTC'	=>	'1,1',
	'CGGCTG'	=>	'1,0',
	'CGGCTT'	=>	'1,1',
	'CGGGAA'	=>	'2,1',
	'CGGGAC'	=>	'2,1',
	'CGGGAG'	=>	'2,0',
	'CGGGAT'	=>	'2,1',
	'CGGGCA'	=>	'2,1',
	'CGGGCC'	=>	'2,1',
	'CGGGCG'	=>	'2,0',
	'CGGGCT'	=>	'2,1',
	'CGGGGA'	=>	'1,1',
	'CGGGGC'	=>	'1,1',
	'CGGGGG'	=>	'1,0',
	'CGGGGT'	=>	'1,1',
	'CGGGTA'	=>	'2,1',
	'CGGGTC'	=>	'2,1',
	'CGGGTG'	=>	'2,0',
	'CGGGTT'	=>	'2,1',
	'CGGTAC'	=>	'2,1',
	'CGGTAT'	=>	'2,1',
	'CGGTCA'	=>	'2,1',
	'CGGTCC'	=>	'2,1',
	'CGGTCG'	=>	'2,0',
	'CGGTCT'	=>	'2,1',
	'CGGTGC'	=>	'1,1',
	'CGGTGG'	=>	'1,0',
	'CGGTGT'	=>	'1,1',
	'CGGTTA'	=>	'1,2',
	'CGGTTC'	=>	'2,1',
	'CGGTTG'	=>	'1,1',
	'CGGTTT'	=>	'2,1',
	'CGTAAA'	=>	'1,2',
	'CGTAAC'	=>	'2,1',
	'CGTAAG'	=>	'1,2',
	'CGTAAT'	=>	'2,0',
	'CGTACA'	=>	'1,2',
	'CGTACC'	=>	'2,1',
	'CGTACG'	=>	'1,2',
	'CGTACT'	=>	'2,0',
	'CGTAGA'	=>	'0,2',
	'CGTAGC'	=>	'1,1',
	'CGTAGG'	=>	'0,2',
	'CGTAGT'	=>	'1,0',
	'CGTATA'	=>	'1,2',
	'CGTATC'	=>	'2,1',
	'CGTATG'	=>	'1,2',
	'CGTATT'	=>	'2,0',
	'CGTCAA'	=>	'1,1',
	'CGTCAC'	=>	'1,1',
	'CGTCAG'	=>	'1,1',
	'CGTCAT'	=>	'1,0',
	'CGTCCA'	=>	'1,1',
	'CGTCCC'	=>	'1,1',
	'CGTCCG'	=>	'1,1',
	'CGTCCT'	=>	'1,0',
	'CGTCGA'	=>	'0,1',
	'CGTCGC'	=>	'0,1',
	'CGTCGG'	=>	'0,1',
	'CGTCGT'	=>	'0,0',
	'CGTCTA'	=>	'1,1',
	'CGTCTC'	=>	'1,1',
	'CGTCTG'	=>	'1,1',
	'CGTCTT'	=>	'1,0',
	'CGTGAA'	=>	'2,1',
	'CGTGAC'	=>	'2,1',
	'CGTGAG'	=>	'2,1',
	'CGTGAT'	=>	'2,0',
	'CGTGCA'	=>	'2,1',
	'CGTGCC'	=>	'2,1',
	'CGTGCG'	=>	'2,1',
	'CGTGCT'	=>	'2,0',
	'CGTGGA'	=>	'1,1',
	'CGTGGC'	=>	'1,1',
	'CGTGGG'	=>	'1,1',
	'CGTGGT'	=>	'1,0',
	'CGTGTA'	=>	'2,1',
	'CGTGTC'	=>	'2,1',
	'CGTGTG'	=>	'2,1',
	'CGTGTT'	=>	'2,0',
	'CGTTAC'	=>	'2,1',
	'CGTTAT'	=>	'2,0',
	'CGTTCA'	=>	'2,1',
	'CGTTCC'	=>	'2,1',
	'CGTTCG'	=>	'2,1',
	'CGTTCT'	=>	'2,0',
	'CGTTGC'	=>	'1,1',
	'CGTTGG'	=>	'1,1',
	'CGTTGT'	=>	'1,0',
	'CGTTTA'	=>	'1,2',
	'CGTTTC'	=>	'2,1',
	'CGTTTG'	=>	'1,2',
	'CGTTTT'	=>	'2,0',
	'CTAAAA'	=>	'2,0',
	'CTAAAC'	=>	'2,1',
	'CTAAAG'	=>	'2,1',
	'CTAAAT'	=>	'2,1',
	'CTAACA'	=>	'2,0',
	'CTAACC'	=>	'2,1',
	'CTAACG'	=>	'2,1',
	'CTAACT'	=>	'2,1',
	'CTAAGA'	=>	'1,1',
	'CTAAGC'	=>	'2,1',
	'CTAAGG'	=>	'1,2',
	'CTAAGT'	=>	'2,1',
	'CTAATA'	=>	'1,0',
	'CTAATC'	=>	'1,1',
	'CTAATG'	=>	'1,1',
	'CTAATT'	=>	'1,1',
	'CTACAA'	=>	'1,0',
	'CTACAC'	=>	'1,1',
	'CTACAG'	=>	'1,1',
	'CTACAT'	=>	'1,1',
	'CTACCA'	=>	'1,0',
	'CTACCC'	=>	'1,1',
	'CTACCG'	=>	'1,1',
	'CTACCT'	=>	'1,1',
	'CTACGA'	=>	'1,0',
	'CTACGC'	=>	'1,1',
	'CTACGG'	=>	'1,1',
	'CTACGT'	=>	'1,1',
	'CTACTA'	=>	'0,0',
	'CTACTC'	=>	'0,1',
	'CTACTG'	=>	'0,1',
	'CTACTT'	=>	'0,1',
	'CTAGAA'	=>	'2,0',
	'CTAGAC'	=>	'2,1',
	'CTAGAG'	=>	'2,1',
	'CTAGAT'	=>	'2,1',
	'CTAGCA'	=>	'2,0',
	'CTAGCC'	=>	'2,1',
	'CTAGCG'	=>	'2,1',
	'CTAGCT'	=>	'2,1',
	'CTAGGA'	=>	'2,0',
	'CTAGGC'	=>	'2,1',
	'CTAGGG'	=>	'2,1',
	'CTAGGT'	=>	'2,1',
	'CTAGTA'	=>	'1,0',
	'CTAGTC'	=>	'1,1',
	'CTAGTG'	=>	'1,1',
	'CTAGTT'	=>	'1,1',
	'CTATAC'	=>	'2,1',
	'CTATAT'	=>	'2,1',
	'CTATCA'	=>	'1,1',
	'CTATCC'	=>	'1,2',
	'CTATCG'	=>	'1,2',
	'CTATCT'	=>	'1,2',
	'CTATGC'	=>	'2,1',
	'CTATGG'	=>	'1,2',
	'CTATGT'	=>	'2,1',
	'CTATTA'	=>	'0,1',
	'CTATTC'	=>	'1,1',
	'CTATTG'	=>	'0,2',
	'CTATTT'	=>	'1,1',
	'CTCAAA'	=>	'2,1',
	'CTCAAC'	=>	'2,0',
	'CTCAAG'	=>	'2,1',
	'CTCAAT'	=>	'2,1',
	'CTCACA'	=>	'2,1',
	'CTCACC'	=>	'2,0',
	'CTCACG'	=>	'2,1',
	'CTCACT'	=>	'2,1',
	'CTCAGA'	=>	'1,2',
	'CTCAGC'	=>	'2,0',
	'CTCAGG'	=>	'1,2',
	'CTCAGT'	=>	'2,1',
	'CTCATA'	=>	'1,1',
	'CTCATC'	=>	'1,0',
	'CTCATG'	=>	'1,1',
	'CTCATT'	=>	'1,1',
	'CTCCAA'	=>	'1,1',
	'CTCCAC'	=>	'1,0',
	'CTCCAG'	=>	'1,1',
	'CTCCAT'	=>	'1,1',
	'CTCCCA'	=>	'1,1',
	'CTCCCC'	=>	'1,0',
	'CTCCCG'	=>	'1,1',
	'CTCCCT'	=>	'1,1',
	'CTCCGA'	=>	'1,1',
	'CTCCGC'	=>	'1,0',
	'CTCCGG'	=>	'1,1',
	'CTCCGT'	=>	'1,1',
	'CTCCTA'	=>	'0,1',
	'CTCCTC'	=>	'0,0',
	'CTCCTG'	=>	'0,1',
	'CTCCTT'	=>	'0,1',
	'CTCGAA'	=>	'2,1',
	'CTCGAC'	=>	'2,0',
	'CTCGAG'	=>	'2,1',
	'CTCGAT'	=>	'2,1',
	'CTCGCA'	=>	'2,1',
	'CTCGCC'	=>	'2,0',
	'CTCGCG'	=>	'2,1',
	'CTCGCT'	=>	'2,1',
	'CTCGGA'	=>	'2,1',
	'CTCGGC'	=>	'2,0',
	'CTCGGG'	=>	'2,1',
	'CTCGGT'	=>	'2,1',
	'CTCGTA'	=>	'1,1',
	'CTCGTC'	=>	'1,0',
	'CTCGTG'	=>	'1,1',
	'CTCGTT'	=>	'1,1',
	'CTCTAC'	=>	'2,0',
	'CTCTAT'	=>	'2,1',
	'CTCTCA'	=>	'1,2',
	'CTCTCC'	=>	'2,0',
	'CTCTCG'	=>	'1,2',
	'CTCTCT'	=>	'2,1',
	'CTCTGC'	=>	'2,0',
	'CTCTGG'	=>	'1,2',
	'CTCTGT'	=>	'2,1',
	'CTCTTA'	=>	'0,2',
	'CTCTTC'	=>	'1,0',
	'CTCTTG'	=>	'0,2',
	'CTCTTT'	=>	'1,1',
	'CTGAAA'	=>	'2,1',
	'CTGAAC'	=>	'2,1',
	'CTGAAG'	=>	'2,0',
	'CTGAAT'	=>	'2,1',
	'CTGACA'	=>	'2,1',
	'CTGACC'	=>	'2,1',
	'CTGACG'	=>	'2,0',
	'CTGACT'	=>	'2,1',
	'CTGAGA'	=>	'1,2',
	'CTGAGC'	=>	'2,1',
	'CTGAGG'	=>	'1,1',
	'CTGAGT'	=>	'2,1',
	'CTGATA'	=>	'1,1',
	'CTGATC'	=>	'1,1',
	'CTGATG'	=>	'1,0',
	'CTGATT'	=>	'1,1',
	'CTGCAA'	=>	'1,1',
	'CTGCAC'	=>	'1,1',
	'CTGCAG'	=>	'1,0',
	'CTGCAT'	=>	'1,1',
	'CTGCCA'	=>	'1,1',
	'CTGCCC'	=>	'1,1',
	'CTGCCG'	=>	'1,0',
	'CTGCCT'	=>	'1,1',
	'CTGCGA'	=>	'1,1',
	'CTGCGC'	=>	'1,1',
	'CTGCGG'	=>	'1,0',
	'CTGCGT'	=>	'1,1',
	'CTGCTA'	=>	'0,1',
	'CTGCTC'	=>	'0,1',
	'CTGCTG'	=>	'0,0',
	'CTGCTT'	=>	'0,1',
	'CTGGAA'	=>	'2,1',
	'CTGGAC'	=>	'2,1',
	'CTGGAG'	=>	'2,0',
	'CTGGAT'	=>	'2,1',
	'CTGGCA'	=>	'2,1',
	'CTGGCC'	=>	'2,1',
	'CTGGCG'	=>	'2,0',
	'CTGGCT'	=>	'2,1',
	'CTGGGA'	=>	'2,1',
	'CTGGGC'	=>	'2,1',
	'CTGGGG'	=>	'2,0',
	'CTGGGT'	=>	'2,1',
	'CTGGTA'	=>	'1,1',
	'CTGGTC'	=>	'1,1',
	'CTGGTG'	=>	'1,0',
	'CTGGTT'	=>	'1,1',
	'CTGTAC'	=>	'2,1',
	'CTGTAT'	=>	'2,1',
	'CTGTCA'	=>	'1,2',
	'CTGTCC'	=>	'1,2',
	'CTGTCG'	=>	'1,1',
	'CTGTCT'	=>	'1,2',
	'CTGTGC'	=>	'2,1',
	'CTGTGG'	=>	'1,1',
	'CTGTGT'	=>	'2,1',
	'CTGTTA'	=>	'0,2',
	'CTGTTC'	=>	'1,1',
	'CTGTTG'	=>	'0,1',
	'CTGTTT'	=>	'1,1',
	'CTTAAA'	=>	'2,1',
	'CTTAAC'	=>	'2,1',
	'CTTAAG'	=>	'2,1',
	'CTTAAT'	=>	'2,0',
	'CTTACA'	=>	'2,1',
	'CTTACC'	=>	'2,1',
	'CTTACG'	=>	'2,1',
	'CTTACT'	=>	'2,0',
	'CTTAGA'	=>	'1,2',
	'CTTAGC'	=>	'2,1',
	'CTTAGG'	=>	'1,2',
	'CTTAGT'	=>	'2,0',
	'CTTATA'	=>	'1,1',
	'CTTATC'	=>	'1,1',
	'CTTATG'	=>	'1,1',
	'CTTATT'	=>	'1,0',
	'CTTCAA'	=>	'1,1',
	'CTTCAC'	=>	'1,1',
	'CTTCAG'	=>	'1,1',
	'CTTCAT'	=>	'1,0',
	'CTTCCA'	=>	'1,1',
	'CTTCCC'	=>	'1,1',
	'CTTCCG'	=>	'1,1',
	'CTTCCT'	=>	'1,0',
	'CTTCGA'	=>	'1,1',
	'CTTCGC'	=>	'1,1',
	'CTTCGG'	=>	'1,1',
	'CTTCGT'	=>	'1,0',
	'CTTCTA'	=>	'0,1',
	'CTTCTC'	=>	'0,1',
	'CTTCTG'	=>	'0,1',
	'CTTCTT'	=>	'0,0',
	'CTTGAA'	=>	'2,1',
	'CTTGAC'	=>	'2,1',
	'CTTGAG'	=>	'2,1',
	'CTTGAT'	=>	'2,0',
	'CTTGCA'	=>	'2,1',
	'CTTGCC'	=>	'2,1',
	'CTTGCG'	=>	'2,1',
	'CTTGCT'	=>	'2,0',
	'CTTGGA'	=>	'2,1',
	'CTTGGC'	=>	'2,1',
	'CTTGGG'	=>	'2,1',
	'CTTGGT'	=>	'2,0',
	'CTTGTA'	=>	'1,1',
	'CTTGTC'	=>	'1,1',
	'CTTGTG'	=>	'1,1',
	'CTTGTT'	=>	'1,0',
	'CTTTAC'	=>	'2,1',
	'CTTTAT'	=>	'2,0',
	'CTTTCA'	=>	'1,2',
	'CTTTCC'	=>	'2,1',
	'CTTTCG'	=>	'1,2',
	'CTTTCT'	=>	'2,0',
	'CTTTGC'	=>	'2,1',
	'CTTTGG'	=>	'1,2',
	'CTTTGT'	=>	'2,0',
	'CTTTTA'	=>	'0,2',
	'CTTTTC'	=>	'1,1',
	'CTTTTG'	=>	'0,2',
	'CTTTTT'	=>	'1,0',
	'GAAAAA'	=>	'1,0',
	'GAAAAC'	=>	'2,0',
	'GAAAAG'	=>	'1,1',
	'GAAAAT'	=>	'2,0',
	'GAAACA'	=>	'2,0',
	'GAAACC'	=>	'2,1',
	'GAAACG'	=>	'2,1',
	'GAAACT'	=>	'2,1',
	'GAAAGA'	=>	'2,0',
	'GAAAGC'	=>	'2,1',
	'GAAAGG'	=>	'2,1',
	'GAAAGT'	=>	'2,1',
	'GAAATA'	=>	'2,0',
	'GAAATC'	=>	'2,1',
	'GAAATG'	=>	'2,1',
	'GAAATT'	=>	'2,1',
	'GAACAA'	=>	'1,0',
	'GAACAC'	=>	'2,0',
	'GAACAG'	=>	'1,1',
	'GAACAT'	=>	'2,0',
	'GAACCA'	=>	'2,0',
	'GAACCC'	=>	'2,1',
	'GAACCG'	=>	'2,1',
	'GAACCT'	=>	'2,1',
	'GAACGA'	=>	'2,0',
	'GAACGC'	=>	'2,1',
	'GAACGG'	=>	'2,1',
	'GAACGT'	=>	'2,1',
	'GAACTA'	=>	'2,0',
	'GAACTC'	=>	'2,1',
	'GAACTG'	=>	'2,1',
	'GAACTT'	=>	'2,1',
	'GAAGAA'	=>	'0,0',
	'GAAGAC'	=>	'1,0',
	'GAAGAG'	=>	'0,1',
	'GAAGAT'	=>	'1,0',
	'GAAGCA'	=>	'1,0',
	'GAAGCC'	=>	'1,1',
	'GAAGCG'	=>	'1,1',
	'GAAGCT'	=>	'1,1',
	'GAAGGA'	=>	'1,0',
	'GAAGGC'	=>	'1,1',
	'GAAGGG'	=>	'1,1',
	'GAAGGT'	=>	'1,1',
	'GAAGTA'	=>	'1,0',
	'GAAGTC'	=>	'1,1',
	'GAAGTG'	=>	'1,1',
	'GAAGTT'	=>	'1,1',
	'GAATAC'	=>	'2,0',
	'GAATAT'	=>	'2,0',
	'GAATCA'	=>	'2,0',
	'GAATCC'	=>	'2,1',
	'GAATCG'	=>	'2,1',
	'GAATCT'	=>	'2,1',
	'GAATGC'	=>	'2,1',
	'GAATGG'	=>	'2,1',
	'GAATGT'	=>	'2,1',
	'GAATTA'	=>	'2,0',
	'GAATTC'	=>	'2,1',
	'GAATTG'	=>	'2,1',
	'GAATTT'	=>	'2,1',
	'GACAAA'	=>	'2,0',
	'GACAAC'	=>	'1,0',
	'GACAAG'	=>	'2,0',
	'GACAAT'	=>	'1,1',
	'GACACA'	=>	'2,1',
	'GACACC'	=>	'2,0',
	'GACACG'	=>	'2,1',
	'GACACT'	=>	'2,1',
	'GACAGA'	=>	'2,1',
	'GACAGC'	=>	'2,0',
	'GACAGG'	=>	'2,1',
	'GACAGT'	=>	'2,1',
	'GACATA'	=>	'2,1',
	'GACATC'	=>	'2,0',
	'GACATG'	=>	'2,1',
	'GACATT'	=>	'2,1',
	'GACCAA'	=>	'2,0',
	'GACCAC'	=>	'1,0',
	'GACCAG'	=>	'2,0',
	'GACCAT'	=>	'1,1',
	'GACCCA'	=>	'2,1',
	'GACCCC'	=>	'2,0',
	'GACCCG'	=>	'2,1',
	'GACCCT'	=>	'2,1',
	'GACCGA'	=>	'2,1',
	'GACCGC'	=>	'2,0',
	'GACCGG'	=>	'2,1',
	'GACCGT'	=>	'2,1',
	'GACCTA'	=>	'2,1',
	'GACCTC'	=>	'2,0',
	'GACCTG'	=>	'2,1',
	'GACCTT'	=>	'2,1',
	'GACGAA'	=>	'1,0',
	'GACGAC'	=>	'0,0',
	'GACGAG'	=>	'1,0',
	'GACGAT'	=>	'0,1',
	'GACGCA'	=>	'1,1',
	'GACGCC'	=>	'1,0',
	'GACGCG'	=>	'1,1',
	'GACGCT'	=>	'1,1',
	'GACGGA'	=>	'1,1',
	'GACGGC'	=>	'1,0',
	'GACGGG'	=>	'1,1',
	'GACGGT'	=>	'1,1',
	'GACGTA'	=>	'1,1',
	'GACGTC'	=>	'1,0',
	'GACGTG'	=>	'1,1',
	'GACGTT'	=>	'1,1',
	'GACTAC'	=>	'1,0',
	'GACTAT'	=>	'1,1',
	'GACTCA'	=>	'2,1',
	'GACTCC'	=>	'2,0',
	'GACTCG'	=>	'2,1',
	'GACTCT'	=>	'2,1',
	'GACTGC'	=>	'2,0',
	'GACTGG'	=>	'2,1',
	'GACTGT'	=>	'2,1',
	'GACTTA'	=>	'2,1',
	'GACTTC'	=>	'2,0',
	'GACTTG'	=>	'2,1',
	'GACTTT'	=>	'2,1',
	'GAGAAA'	=>	'1,1',
	'GAGAAC'	=>	'2,0',
	'GAGAAG'	=>	'1,0',
	'GAGAAT'	=>	'2,0',
	'GAGACA'	=>	'2,1',
	'GAGACC'	=>	'2,1',
	'GAGACG'	=>	'2,0',
	'GAGACT'	=>	'2,1',
	'GAGAGA'	=>	'2,1',
	'GAGAGC'	=>	'2,1',
	'GAGAGG'	=>	'2,0',
	'GAGAGT'	=>	'2,1',
	'GAGATA'	=>	'2,1',
	'GAGATC'	=>	'2,1',
	'GAGATG'	=>	'2,0',
	'GAGATT'	=>	'2,1',
	'GAGCAA'	=>	'1,1',
	'GAGCAC'	=>	'2,0',
	'GAGCAG'	=>	'1,0',
	'GAGCAT'	=>	'2,0',
	'GAGCCA'	=>	'2,1',
	'GAGCCC'	=>	'2,1',
	'GAGCCG'	=>	'2,0',
	'GAGCCT'	=>	'2,1',
	'GAGCGA'	=>	'2,1',
	'GAGCGC'	=>	'2,1',
	'GAGCGG'	=>	'2,0',
	'GAGCGT'	=>	'2,1',
	'GAGCTA'	=>	'2,1',
	'GAGCTC'	=>	'2,1',
	'GAGCTG'	=>	'2,0',
	'GAGCTT'	=>	'2,1',
	'GAGGAA'	=>	'0,1',
	'GAGGAC'	=>	'1,0',
	'GAGGAG'	=>	'0,0',
	'GAGGAT'	=>	'1,0',
	'GAGGCA'	=>	'1,1',
	'GAGGCC'	=>	'1,1',
	'GAGGCG'	=>	'1,0',
	'GAGGCT'	=>	'1,1',
	'GAGGGA'	=>	'1,1',
	'GAGGGC'	=>	'1,1',
	'GAGGGG'	=>	'1,0',
	'GAGGGT'	=>	'1,1',
	'GAGGTA'	=>	'1,1',
	'GAGGTC'	=>	'1,1',
	'GAGGTG'	=>	'1,0',
	'GAGGTT'	=>	'1,1',
	'GAGTAC'	=>	'2,0',
	'GAGTAT'	=>	'2,0',
	'GAGTCA'	=>	'2,1',
	'GAGTCC'	=>	'2,1',
	'GAGTCG'	=>	'2,0',
	'GAGTCT'	=>	'2,1',
	'GAGTGC'	=>	'2,1',
	'GAGTGG'	=>	'2,0',
	'GAGTGT'	=>	'2,1',
	'GAGTTA'	=>	'2,1',
	'GAGTTC'	=>	'2,1',
	'GAGTTG'	=>	'2,0',
	'GAGTTT'	=>	'2,1',
	'GATAAA'	=>	'2,0',
	'GATAAC'	=>	'1,1',
	'GATAAG'	=>	'2,0',
	'GATAAT'	=>	'1,0',
	'GATACA'	=>	'2,1',
	'GATACC'	=>	'2,1',
	'GATACG'	=>	'2,1',
	'GATACT'	=>	'2,0',
	'GATAGA'	=>	'2,1',
	'GATAGC'	=>	'2,1',
	'GATAGG'	=>	'2,1',
	'GATAGT'	=>	'2,0',
	'GATATA'	=>	'2,1',
	'GATATC'	=>	'2,1',
	'GATATG'	=>	'2,1',
	'GATATT'	=>	'2,0',
	'GATCAA'	=>	'2,0',
	'GATCAC'	=>	'1,1',
	'GATCAG'	=>	'2,0',
	'GATCAT'	=>	'1,0',
	'GATCCA'	=>	'2,1',
	'GATCCC'	=>	'2,1',
	'GATCCG'	=>	'2,1',
	'GATCCT'	=>	'2,0',
	'GATCGA'	=>	'2,1',
	'GATCGC'	=>	'2,1',
	'GATCGG'	=>	'2,1',
	'GATCGT'	=>	'2,0',
	'GATCTA'	=>	'2,1',
	'GATCTC'	=>	'2,1',
	'GATCTG'	=>	'2,1',
	'GATCTT'	=>	'2,0',
	'GATGAA'	=>	'1,0',
	'GATGAC'	=>	'0,1',
	'GATGAG'	=>	'1,0',
	'GATGAT'	=>	'0,0',
	'GATGCA'	=>	'1,1',
	'GATGCC'	=>	'1,1',
	'GATGCG'	=>	'1,1',
	'GATGCT'	=>	'1,0',
	'GATGGA'	=>	'1,1',
	'GATGGC'	=>	'1,1',
	'GATGGG'	=>	'1,1',
	'GATGGT'	=>	'1,0',
	'GATGTA'	=>	'1,1',
	'GATGTC'	=>	'1,1',
	'GATGTG'	=>	'1,1',
	'GATGTT'	=>	'1,0',
	'GATTAC'	=>	'1,1',
	'GATTAT'	=>	'1,0',
	'GATTCA'	=>	'2,1',
	'GATTCC'	=>	'2,1',
	'GATTCG'	=>	'2,1',
	'GATTCT'	=>	'2,0',
	'GATTGC'	=>	'2,1',
	'GATTGG'	=>	'2,1',
	'GATTGT'	=>	'2,0',
	'GATTTA'	=>	'2,1',
	'GATTTC'	=>	'2,1',
	'GATTTG'	=>	'2,1',
	'GATTTT'	=>	'2,0',
	'GCAAAA'	=>	'2,0',
	'GCAAAC'	=>	'2,1',
	'GCAAAG'	=>	'2,1',
	'GCAAAT'	=>	'2,1',
	'GCAACA'	=>	'1,0',
	'GCAACC'	=>	'1,1',
	'GCAACG'	=>	'1,1',
	'GCAACT'	=>	'1,1',
	'GCAAGA'	=>	'2,0',
	'GCAAGC'	=>	'2,1',
	'GCAAGG'	=>	'2,1',
	'GCAAGT'	=>	'2,1',
	'GCAATA'	=>	'2,0',
	'GCAATC'	=>	'2,1',
	'GCAATG'	=>	'2,1',
	'GCAATT'	=>	'2,1',
	'GCACAA'	=>	'2,0',
	'GCACAC'	=>	'2,1',
	'GCACAG'	=>	'2,1',
	'GCACAT'	=>	'2,1',
	'GCACCA'	=>	'1,0',
	'GCACCC'	=>	'1,1',
	'GCACCG'	=>	'1,1',
	'GCACCT'	=>	'1,1',
	'GCACGA'	=>	'2,0',
	'GCACGC'	=>	'2,1',
	'GCACGG'	=>	'2,1',
	'GCACGT'	=>	'2,1',
	'GCACTA'	=>	'2,0',
	'GCACTC'	=>	'2,1',
	'GCACTG'	=>	'2,1',
	'GCACTT'	=>	'2,1',
	'GCAGAA'	=>	'1,0',
	'GCAGAC'	=>	'1,1',
	'GCAGAG'	=>	'1,1',
	'GCAGAT'	=>	'1,1',
	'GCAGCA'	=>	'0,0',
	'GCAGCC'	=>	'0,1',
	'GCAGCG'	=>	'0,1',
	'GCAGCT'	=>	'0,1',
	'GCAGGA'	=>	'1,0',
	'GCAGGC'	=>	'1,1',
	'GCAGGG'	=>	'1,1',
	'GCAGGT'	=>	'1,1',
	'GCAGTA'	=>	'1,0',
	'GCAGTC'	=>	'1,1',
	'GCAGTG'	=>	'1,1',
	'GCAGTT'	=>	'1,1',
	'GCATAC'	=>	'2,1',
	'GCATAT'	=>	'2,1',
	'GCATCA'	=>	'1,0',
	'GCATCC'	=>	'1,1',
	'GCATCG'	=>	'1,1',
	'GCATCT'	=>	'1,1',
	'GCATGC'	=>	'2,1',
	'GCATGG'	=>	'2,1',
	'GCATGT'	=>	'2,1',
	'GCATTA'	=>	'2,0',
	'GCATTC'	=>	'2,1',
	'GCATTG'	=>	'2,1',
	'GCATTT'	=>	'2,1',
	'GCCAAA'	=>	'2,1',
	'GCCAAC'	=>	'2,0',
	'GCCAAG'	=>	'2,1',
	'GCCAAT'	=>	'2,1',
	'GCCACA'	=>	'1,1',
	'GCCACC'	=>	'1,0',
	'GCCACG'	=>	'1,1',
	'GCCACT'	=>	'1,1',
	'GCCAGA'	=>	'2,1',
	'GCCAGC'	=>	'2,0',
	'GCCAGG'	=>	'2,1',
	'GCCAGT'	=>	'2,1',
	'GCCATA'	=>	'2,1',
	'GCCATC'	=>	'2,0',
	'GCCATG'	=>	'2,1',
	'GCCATT'	=>	'2,1',
	'GCCCAA'	=>	'2,1',
	'GCCCAC'	=>	'2,0',
	'GCCCAG'	=>	'2,1',
	'GCCCAT'	=>	'2,1',
	'GCCCCA'	=>	'1,1',
	'GCCCCC'	=>	'1,0',
	'GCCCCG'	=>	'1,1',
	'GCCCCT'	=>	'1,1',
	'GCCCGA'	=>	'2,1',
	'GCCCGC'	=>	'2,0',
	'GCCCGG'	=>	'2,1',
	'GCCCGT'	=>	'2,1',
	'GCCCTA'	=>	'2,1',
	'GCCCTC'	=>	'2,0',
	'GCCCTG'	=>	'2,1',
	'GCCCTT'	=>	'2,1',
	'GCCGAA'	=>	'1,1',
	'GCCGAC'	=>	'1,0',
	'GCCGAG'	=>	'1,1',
	'GCCGAT'	=>	'1,1',
	'GCCGCA'	=>	'0,1',
	'GCCGCC'	=>	'0,0',
	'GCCGCG'	=>	'0,1',
	'GCCGCT'	=>	'0,1',
	'GCCGGA'	=>	'1,1',
	'GCCGGC'	=>	'1,0',
	'GCCGGG'	=>	'1,1',
	'GCCGGT'	=>	'1,1',
	'GCCGTA'	=>	'1,1',
	'GCCGTC'	=>	'1,0',
	'GCCGTG'	=>	'1,1',
	'GCCGTT'	=>	'1,1',
	'GCCTAC'	=>	'2,0',
	'GCCTAT'	=>	'2,1',
	'GCCTCA'	=>	'1,1',
	'GCCTCC'	=>	'1,0',
	'GCCTCG'	=>	'1,1',
	'GCCTCT'	=>	'1,1',
	'GCCTGC'	=>	'2,0',
	'GCCTGG'	=>	'2,1',
	'GCCTGT'	=>	'2,1',
	'GCCTTA'	=>	'2,1',
	'GCCTTC'	=>	'2,0',
	'GCCTTG'	=>	'2,1',
	'GCCTTT'	=>	'2,1',
	'GCGAAA'	=>	'2,1',
	'GCGAAC'	=>	'2,1',
	'GCGAAG'	=>	'2,0',
	'GCGAAT'	=>	'2,1',
	'GCGACA'	=>	'1,1',
	'GCGACC'	=>	'1,1',
	'GCGACG'	=>	'1,0',
	'GCGACT'	=>	'1,1',
	'GCGAGA'	=>	'2,1',
	'GCGAGC'	=>	'2,1',
	'GCGAGG'	=>	'2,0',
	'GCGAGT'	=>	'2,1',
	'GCGATA'	=>	'2,1',
	'GCGATC'	=>	'2,1',
	'GCGATG'	=>	'2,0',
	'GCGATT'	=>	'2,1',
	'GCGCAA'	=>	'2,1',
	'GCGCAC'	=>	'2,1',
	'GCGCAG'	=>	'2,0',
	'GCGCAT'	=>	'2,1',
	'GCGCCA'	=>	'1,1',
	'GCGCCC'	=>	'1,1',
	'GCGCCG'	=>	'1,0',
	'GCGCCT'	=>	'1,1',
	'GCGCGA'	=>	'2,1',
	'GCGCGC'	=>	'2,1',
	'GCGCGG'	=>	'2,0',
	'GCGCGT'	=>	'2,1',
	'GCGCTA'	=>	'2,1',
	'GCGCTC'	=>	'2,1',
	'GCGCTG'	=>	'2,0',
	'GCGCTT'	=>	'2,1',
	'GCGGAA'	=>	'1,1',
	'GCGGAC'	=>	'1,1',
	'GCGGAG'	=>	'1,0',
	'GCGGAT'	=>	'1,1',
	'GCGGCA'	=>	'0,1',
	'GCGGCC'	=>	'0,1',
	'GCGGCG'	=>	'0,0',
	'GCGGCT'	=>	'0,1',
	'GCGGGA'	=>	'1,1',
	'GCGGGC'	=>	'1,1',
	'GCGGGG'	=>	'1,0',
	'GCGGGT'	=>	'1,1',
	'GCGGTA'	=>	'1,1',
	'GCGGTC'	=>	'1,1',
	'GCGGTG'	=>	'1,0',
	'GCGGTT'	=>	'1,1',
	'GCGTAC'	=>	'2,1',
	'GCGTAT'	=>	'2,1',
	'GCGTCA'	=>	'1,1',
	'GCGTCC'	=>	'1,1',
	'GCGTCG'	=>	'1,0',
	'GCGTCT'	=>	'1,1',
	'GCGTGC'	=>	'2,1',
	'GCGTGG'	=>	'2,0',
	'GCGTGT'	=>	'2,1',
	'GCGTTA'	=>	'2,1',
	'GCGTTC'	=>	'2,1',
	'GCGTTG'	=>	'2,0',
	'GCGTTT'	=>	'2,1',
	'GCTAAA'	=>	'2,1',
	'GCTAAC'	=>	'2,1',
	'GCTAAG'	=>	'2,1',
	'GCTAAT'	=>	'2,0',
	'GCTACA'	=>	'1,1',
	'GCTACC'	=>	'1,1',
	'GCTACG'	=>	'1,1',
	'GCTACT'	=>	'1,0',
	'GCTAGA'	=>	'2,1',
	'GCTAGC'	=>	'2,1',
	'GCTAGG'	=>	'2,1',
	'GCTAGT'	=>	'2,0',
	'GCTATA'	=>	'2,1',
	'GCTATC'	=>	'2,1',
	'GCTATG'	=>	'2,1',
	'GCTATT'	=>	'2,0',
	'GCTCAA'	=>	'2,1',
	'GCTCAC'	=>	'2,1',
	'GCTCAG'	=>	'2,1',
	'GCTCAT'	=>	'2,0',
	'GCTCCA'	=>	'1,1',
	'GCTCCC'	=>	'1,1',
	'GCTCCG'	=>	'1,1',
	'GCTCCT'	=>	'1,0',
	'GCTCGA'	=>	'2,1',
	'GCTCGC'	=>	'2,1',
	'GCTCGG'	=>	'2,1',
	'GCTCGT'	=>	'2,0',
	'GCTCTA'	=>	'2,1',
	'GCTCTC'	=>	'2,1',
	'GCTCTG'	=>	'2,1',
	'GCTCTT'	=>	'2,0',
	'GCTGAA'	=>	'1,1',
	'GCTGAC'	=>	'1,1',
	'GCTGAG'	=>	'1,1',
	'GCTGAT'	=>	'1,0',
	'GCTGCA'	=>	'0,1',
	'GCTGCC'	=>	'0,1',
	'GCTGCG'	=>	'0,1',
	'GCTGCT'	=>	'0,0',
	'GCTGGA'	=>	'1,1',
	'GCTGGC'	=>	'1,1',
	'GCTGGG'	=>	'1,1',
	'GCTGGT'	=>	'1,0',
	'GCTGTA'	=>	'1,1',
	'GCTGTC'	=>	'1,1',
	'GCTGTG'	=>	'1,1',
	'GCTGTT'	=>	'1,0',
	'GCTTAC'	=>	'2,1',
	'GCTTAT'	=>	'2,0',
	'GCTTCA'	=>	'1,1',
	'GCTTCC'	=>	'1,1',
	'GCTTCG'	=>	'1,1',
	'GCTTCT'	=>	'1,0',
	'GCTTGC'	=>	'2,1',
	'GCTTGG'	=>	'2,1',
	'GCTTGT'	=>	'2,0',
	'GCTTTA'	=>	'2,1',
	'GCTTTC'	=>	'2,1',
	'GCTTTG'	=>	'2,1',
	'GCTTTT'	=>	'2,0',
	'GGAAAA'	=>	'2,0',
	'GGAAAC'	=>	'2,1',
	'GGAAAG'	=>	'2,1',
	'GGAAAT'	=>	'2,1',
	'GGAACA'	=>	'2,0',
	'GGAACC'	=>	'2,1',
	'GGAACG'	=>	'2,1',
	'GGAACT'	=>	'2,1',
	'GGAAGA'	=>	'1,0',
	'GGAAGC'	=>	'1,1',
	'GGAAGG'	=>	'1,1',
	'GGAAGT'	=>	'1,1',
	'GGAATA'	=>	'2,0',
	'GGAATC'	=>	'2,1',
	'GGAATG'	=>	'2,1',
	'GGAATT'	=>	'2,1',
	'GGACAA'	=>	'2,0',
	'GGACAC'	=>	'2,1',
	'GGACAG'	=>	'2,1',
	'GGACAT'	=>	'2,1',
	'GGACCA'	=>	'2,0',
	'GGACCC'	=>	'2,1',
	'GGACCG'	=>	'2,1',
	'GGACCT'	=>	'2,1',
	'GGACGA'	=>	'1,0',
	'GGACGC'	=>	'1,1',
	'GGACGG'	=>	'1,1',
	'GGACGT'	=>	'1,1',
	'GGACTA'	=>	'2,0',
	'GGACTC'	=>	'2,1',
	'GGACTG'	=>	'2,1',
	'GGACTT'	=>	'2,1',
	'GGAGAA'	=>	'1,0',
	'GGAGAC'	=>	'1,1',
	'GGAGAG'	=>	'1,1',
	'GGAGAT'	=>	'1,1',
	'GGAGCA'	=>	'1,0',
	'GGAGCC'	=>	'1,1',
	'GGAGCG'	=>	'1,1',
	'GGAGCT'	=>	'1,1',
	'GGAGGA'	=>	'0,0',
	'GGAGGC'	=>	'0,1',
	'GGAGGG'	=>	'0,1',
	'GGAGGT'	=>	'0,1',
	'GGAGTA'	=>	'1,0',
	'GGAGTC'	=>	'1,1',
	'GGAGTG'	=>	'1,1',
	'GGAGTT'	=>	'1,1',
	'GGATAC'	=>	'2,1',
	'GGATAT'	=>	'2,1',
	'GGATCA'	=>	'2,0',
	'GGATCC'	=>	'2,1',
	'GGATCG'	=>	'2,1',
	'GGATCT'	=>	'2,1',
	'GGATGC'	=>	'1,1',
	'GGATGG'	=>	'1,1',
	'GGATGT'	=>	'1,1',
	'GGATTA'	=>	'2,0',
	'GGATTC'	=>	'2,1',
	'GGATTG'	=>	'2,1',
	'GGATTT'	=>	'2,1',
	'GGCAAA'	=>	'2,1',
	'GGCAAC'	=>	'2,0',
	'GGCAAG'	=>	'2,1',
	'GGCAAT'	=>	'2,1',
	'GGCACA'	=>	'2,1',
	'GGCACC'	=>	'2,0',
	'GGCACG'	=>	'2,1',
	'GGCACT'	=>	'2,1',
	'GGCAGA'	=>	'1,1',
	'GGCAGC'	=>	'1,0',
	'GGCAGG'	=>	'1,1',
	'GGCAGT'	=>	'1,1',
	'GGCATA'	=>	'2,1',
	'GGCATC'	=>	'2,0',
	'GGCATG'	=>	'2,1',
	'GGCATT'	=>	'2,1',
	'GGCCAA'	=>	'2,1',
	'GGCCAC'	=>	'2,0',
	'GGCCAG'	=>	'2,1',
	'GGCCAT'	=>	'2,1',
	'GGCCCA'	=>	'2,1',
	'GGCCCC'	=>	'2,0',
	'GGCCCG'	=>	'2,1',
	'GGCCCT'	=>	'2,1',
	'GGCCGA'	=>	'1,1',
	'GGCCGC'	=>	'1,0',
	'GGCCGG'	=>	'1,1',
	'GGCCGT'	=>	'1,1',
	'GGCCTA'	=>	'2,1',
	'GGCCTC'	=>	'2,0',
	'GGCCTG'	=>	'2,1',
	'GGCCTT'	=>	'2,1',
	'GGCGAA'	=>	'1,1',
	'GGCGAC'	=>	'1,0',
	'GGCGAG'	=>	'1,1',
	'GGCGAT'	=>	'1,1',
	'GGCGCA'	=>	'1,1',
	'GGCGCC'	=>	'1,0',
	'GGCGCG'	=>	'1,1',
	'GGCGCT'	=>	'1,1',
	'GGCGGA'	=>	'0,1',
	'GGCGGC'	=>	'0,0',
	'GGCGGG'	=>	'0,1',
	'GGCGGT'	=>	'0,1',
	'GGCGTA'	=>	'1,1',
	'GGCGTC'	=>	'1,0',
	'GGCGTG'	=>	'1,1',
	'GGCGTT'	=>	'1,1',
	'GGCTAC'	=>	'2,0',
	'GGCTAT'	=>	'2,1',
	'GGCTCA'	=>	'2,1',
	'GGCTCC'	=>	'2,0',
	'GGCTCG'	=>	'2,1',
	'GGCTCT'	=>	'2,1',
	'GGCTGC'	=>	'1,0',
	'GGCTGG'	=>	'1,1',
	'GGCTGT'	=>	'1,1',
	'GGCTTA'	=>	'2,1',
	'GGCTTC'	=>	'2,0',
	'GGCTTG'	=>	'2,1',
	'GGCTTT'	=>	'2,1',
	'GGGAAA'	=>	'2,1',
	'GGGAAC'	=>	'2,1',
	'GGGAAG'	=>	'2,0',
	'GGGAAT'	=>	'2,1',
	'GGGACA'	=>	'2,1',
	'GGGACC'	=>	'2,1',
	'GGGACG'	=>	'2,0',
	'GGGACT'	=>	'2,1',
	'GGGAGA'	=>	'1,1',
	'GGGAGC'	=>	'1,1',
	'GGGAGG'	=>	'1,0',
	'GGGAGT'	=>	'1,1',
	'GGGATA'	=>	'2,1',
	'GGGATC'	=>	'2,1',
	'GGGATG'	=>	'2,0',
	'GGGATT'	=>	'2,1',
	'GGGCAA'	=>	'2,1',
	'GGGCAC'	=>	'2,1',
	'GGGCAG'	=>	'2,0',
	'GGGCAT'	=>	'2,1',
	'GGGCCA'	=>	'2,1',
	'GGGCCC'	=>	'2,1',
	'GGGCCG'	=>	'2,0',
	'GGGCCT'	=>	'2,1',
	'GGGCGA'	=>	'1,1',
	'GGGCGC'	=>	'1,1',
	'GGGCGG'	=>	'1,0',
	'GGGCGT'	=>	'1,1',
	'GGGCTA'	=>	'2,1',
	'GGGCTC'	=>	'2,1',
	'GGGCTG'	=>	'2,0',
	'GGGCTT'	=>	'2,1',
	'GGGGAA'	=>	'1,1',
	'GGGGAC'	=>	'1,1',
	'GGGGAG'	=>	'1,0',
	'GGGGAT'	=>	'1,1',
	'GGGGCA'	=>	'1,1',
	'GGGGCC'	=>	'1,1',
	'GGGGCG'	=>	'1,0',
	'GGGGCT'	=>	'1,1',
	'GGGGGA'	=>	'0,1',
	'GGGGGC'	=>	'0,1',
	'GGGGGG'	=>	'0,0',
	'GGGGGT'	=>	'0,1',
	'GGGGTA'	=>	'1,1',
	'GGGGTC'	=>	'1,1',
	'GGGGTG'	=>	'1,0',
	'GGGGTT'	=>	'1,1',
	'GGGTAC'	=>	'2,1',
	'GGGTAT'	=>	'2,1',
	'GGGTCA'	=>	'2,1',
	'GGGTCC'	=>	'2,1',
	'GGGTCG'	=>	'2,0',
	'GGGTCT'	=>	'2,1',
	'GGGTGC'	=>	'1,1',
	'GGGTGG'	=>	'1,0',
	'GGGTGT'	=>	'1,1',
	'GGGTTA'	=>	'2,1',
	'GGGTTC'	=>	'2,1',
	'GGGTTG'	=>	'2,0',
	'GGGTTT'	=>	'2,1',
	'GGTAAA'	=>	'2,1',
	'GGTAAC'	=>	'2,1',
	'GGTAAG'	=>	'2,1',
	'GGTAAT'	=>	'2,0',
	'GGTACA'	=>	'2,1',
	'GGTACC'	=>	'2,1',
	'GGTACG'	=>	'2,1',
	'GGTACT'	=>	'2,0',
	'GGTAGA'	=>	'1,1',
	'GGTAGC'	=>	'1,1',
	'GGTAGG'	=>	'1,1',
	'GGTAGT'	=>	'1,0',
	'GGTATA'	=>	'2,1',
	'GGTATC'	=>	'2,1',
	'GGTATG'	=>	'2,1',
	'GGTATT'	=>	'2,0',
	'GGTCAA'	=>	'2,1',
	'GGTCAC'	=>	'2,1',
	'GGTCAG'	=>	'2,1',
	'GGTCAT'	=>	'2,0',
	'GGTCCA'	=>	'2,1',
	'GGTCCC'	=>	'2,1',
	'GGTCCG'	=>	'2,1',
	'GGTCCT'	=>	'2,0',
	'GGTCGA'	=>	'1,1',
	'GGTCGC'	=>	'1,1',
	'GGTCGG'	=>	'1,1',
	'GGTCGT'	=>	'1,0',
	'GGTCTA'	=>	'2,1',
	'GGTCTC'	=>	'2,1',
	'GGTCTG'	=>	'2,1',
	'GGTCTT'	=>	'2,0',
	'GGTGAA'	=>	'1,1',
	'GGTGAC'	=>	'1,1',
	'GGTGAG'	=>	'1,1',
	'GGTGAT'	=>	'1,0',
	'GGTGCA'	=>	'1,1',
	'GGTGCC'	=>	'1,1',
	'GGTGCG'	=>	'1,1',
	'GGTGCT'	=>	'1,0',
	'GGTGGA'	=>	'0,1',
	'GGTGGC'	=>	'0,1',
	'GGTGGG'	=>	'0,1',
	'GGTGGT'	=>	'0,0',
	'GGTGTA'	=>	'1,1',
	'GGTGTC'	=>	'1,1',
	'GGTGTG'	=>	'1,1',
	'GGTGTT'	=>	'1,0',
	'GGTTAC'	=>	'2,1',
	'GGTTAT'	=>	'2,0',
	'GGTTCA'	=>	'2,1',
	'GGTTCC'	=>	'2,1',
	'GGTTCG'	=>	'2,1',
	'GGTTCT'	=>	'2,0',
	'GGTTGC'	=>	'1,1',
	'GGTTGG'	=>	'1,1',
	'GGTTGT'	=>	'1,0',
	'GGTTTA'	=>	'2,1',
	'GGTTTC'	=>	'2,1',
	'GGTTTG'	=>	'2,1',
	'GGTTTT'	=>	'2,0',
	'GTAAAA'	=>	'2,0',
	'GTAAAC'	=>	'2,1',
	'GTAAAG'	=>	'2,1',
	'GTAAAT'	=>	'2,1',
	'GTAACA'	=>	'2,0',
	'GTAACC'	=>	'2,1',
	'GTAACG'	=>	'2,1',
	'GTAACT'	=>	'2,1',
	'GTAAGA'	=>	'2,0',
	'GTAAGC'	=>	'2,1',
	'GTAAGG'	=>	'2,1',
	'GTAAGT'	=>	'2,1',
	'GTAATA'	=>	'1,0',
	'GTAATC'	=>	'1,1',
	'GTAATG'	=>	'1,1',
	'GTAATT'	=>	'1,1',
	'GTACAA'	=>	'2,0',
	'GTACAC'	=>	'2,1',
	'GTACAG'	=>	'2,1',
	'GTACAT'	=>	'2,1',
	'GTACCA'	=>	'2,0',
	'GTACCC'	=>	'2,1',
	'GTACCG'	=>	'2,1',
	'GTACCT'	=>	'2,1',
	'GTACGA'	=>	'2,0',
	'GTACGC'	=>	'2,1',
	'GTACGG'	=>	'2,1',
	'GTACGT'	=>	'2,1',
	'GTACTA'	=>	'1,0',
	'GTACTC'	=>	'1,1',
	'GTACTG'	=>	'1,1',
	'GTACTT'	=>	'1,1',
	'GTAGAA'	=>	'1,0',
	'GTAGAC'	=>	'1,1',
	'GTAGAG'	=>	'1,1',
	'GTAGAT'	=>	'1,1',
	'GTAGCA'	=>	'1,0',
	'GTAGCC'	=>	'1,1',
	'GTAGCG'	=>	'1,1',
	'GTAGCT'	=>	'1,1',
	'GTAGGA'	=>	'1,0',
	'GTAGGC'	=>	'1,1',
	'GTAGGG'	=>	'1,1',
	'GTAGGT'	=>	'1,1',
	'GTAGTA'	=>	'0,0',
	'GTAGTC'	=>	'0,1',
	'GTAGTG'	=>	'0,1',
	'GTAGTT'	=>	'0,1',
	'GTATAC'	=>	'2,1',
	'GTATAT'	=>	'2,1',
	'GTATCA'	=>	'2,0',
	'GTATCC'	=>	'2,1',
	'GTATCG'	=>	'2,1',
	'GTATCT'	=>	'2,1',
	'GTATGC'	=>	'2,1',
	'GTATGG'	=>	'2,1',
	'GTATGT'	=>	'2,1',
	'GTATTA'	=>	'1,0',
	'GTATTC'	=>	'1,1',
	'GTATTG'	=>	'1,1',
	'GTATTT'	=>	'1,1',
	'GTCAAA'	=>	'2,1',
	'GTCAAC'	=>	'2,0',
	'GTCAAG'	=>	'2,1',
	'GTCAAT'	=>	'2,1',
	'GTCACA'	=>	'2,1',
	'GTCACC'	=>	'2,0',
	'GTCACG'	=>	'2,1',
	'GTCACT'	=>	'2,1',
	'GTCAGA'	=>	'2,1',
	'GTCAGC'	=>	'2,0',
	'GTCAGG'	=>	'2,1',
	'GTCAGT'	=>	'2,1',
	'GTCATA'	=>	'1,1',
	'GTCATC'	=>	'1,0',
	'GTCATG'	=>	'1,1',
	'GTCATT'	=>	'1,1',
	'GTCCAA'	=>	'2,1',
	'GTCCAC'	=>	'2,0',
	'GTCCAG'	=>	'2,1',
	'GTCCAT'	=>	'2,1',
	'GTCCCA'	=>	'2,1',
	'GTCCCC'	=>	'2,0',
	'GTCCCG'	=>	'2,1',
	'GTCCCT'	=>	'2,1',
	'GTCCGA'	=>	'2,1',
	'GTCCGC'	=>	'2,0',
	'GTCCGG'	=>	'2,1',
	'GTCCGT'	=>	'2,1',
	'GTCCTA'	=>	'1,1',
	'GTCCTC'	=>	'1,0',
	'GTCCTG'	=>	'1,1',
	'GTCCTT'	=>	'1,1',
	'GTCGAA'	=>	'1,1',
	'GTCGAC'	=>	'1,0',
	'GTCGAG'	=>	'1,1',
	'GTCGAT'	=>	'1,1',
	'GTCGCA'	=>	'1,1',
	'GTCGCC'	=>	'1,0',
	'GTCGCG'	=>	'1,1',
	'GTCGCT'	=>	'1,1',
	'GTCGGA'	=>	'1,1',
	'GTCGGC'	=>	'1,0',
	'GTCGGG'	=>	'1,1',
	'GTCGGT'	=>	'1,1',
	'GTCGTA'	=>	'0,1',
	'GTCGTC'	=>	'0,0',
	'GTCGTG'	=>	'0,1',
	'GTCGTT'	=>	'0,1',
	'GTCTAC'	=>	'2,0',
	'GTCTAT'	=>	'2,1',
	'GTCTCA'	=>	'2,1',
	'GTCTCC'	=>	'2,0',
	'GTCTCG'	=>	'2,1',
	'GTCTCT'	=>	'2,1',
	'GTCTGC'	=>	'2,0',
	'GTCTGG'	=>	'2,1',
	'GTCTGT'	=>	'2,1',
	'GTCTTA'	=>	'1,1',
	'GTCTTC'	=>	'1,0',
	'GTCTTG'	=>	'1,1',
	'GTCTTT'	=>	'1,1',
	'GTGAAA'	=>	'2,1',
	'GTGAAC'	=>	'2,1',
	'GTGAAG'	=>	'2,0',
	'GTGAAT'	=>	'2,1',
	'GTGACA'	=>	'2,1',
	'GTGACC'	=>	'2,1',
	'GTGACG'	=>	'2,0',
	'GTGACT'	=>	'2,1',
	'GTGAGA'	=>	'2,1',
	'GTGAGC'	=>	'2,1',
	'GTGAGG'	=>	'2,0',
	'GTGAGT'	=>	'2,1',
	'GTGATA'	=>	'1,1',
	'GTGATC'	=>	'1,1',
	'GTGATG'	=>	'1,0',
	'GTGATT'	=>	'1,1',
	'GTGCAA'	=>	'2,1',
	'GTGCAC'	=>	'2,1',
	'GTGCAG'	=>	'2,0',
	'GTGCAT'	=>	'2,1',
	'GTGCCA'	=>	'2,1',
	'GTGCCC'	=>	'2,1',
	'GTGCCG'	=>	'2,0',
	'GTGCCT'	=>	'2,1',
	'GTGCGA'	=>	'2,1',
	'GTGCGC'	=>	'2,1',
	'GTGCGG'	=>	'2,0',
	'GTGCGT'	=>	'2,1',
	'GTGCTA'	=>	'1,1',
	'GTGCTC'	=>	'1,1',
	'GTGCTG'	=>	'1,0',
	'GTGCTT'	=>	'1,1',
	'GTGGAA'	=>	'1,1',
	'GTGGAC'	=>	'1,1',
	'GTGGAG'	=>	'1,0',
	'GTGGAT'	=>	'1,1',
	'GTGGCA'	=>	'1,1',
	'GTGGCC'	=>	'1,1',
	'GTGGCG'	=>	'1,0',
	'GTGGCT'	=>	'1,1',
	'GTGGGA'	=>	'1,1',
	'GTGGGC'	=>	'1,1',
	'GTGGGG'	=>	'1,0',
	'GTGGGT'	=>	'1,1',
	'GTGGTA'	=>	'0,1',
	'GTGGTC'	=>	'0,1',
	'GTGGTG'	=>	'0,0',
	'GTGGTT'	=>	'0,1',
	'GTGTAC'	=>	'2,1',
	'GTGTAT'	=>	'2,1',
	'GTGTCA'	=>	'2,1',
	'GTGTCC'	=>	'2,1',
	'GTGTCG'	=>	'2,0',
	'GTGTCT'	=>	'2,1',
	'GTGTGC'	=>	'2,1',
	'GTGTGG'	=>	'2,0',
	'GTGTGT'	=>	'2,1',
	'GTGTTA'	=>	'1,1',
	'GTGTTC'	=>	'1,1',
	'GTGTTG'	=>	'1,0',
	'GTGTTT'	=>	'1,1',
	'GTTAAA'	=>	'2,1',
	'GTTAAC'	=>	'2,1',
	'GTTAAG'	=>	'2,1',
	'GTTAAT'	=>	'2,0',
	'GTTACA'	=>	'2,1',
	'GTTACC'	=>	'2,1',
	'GTTACG'	=>	'2,1',
	'GTTACT'	=>	'2,0',
	'GTTAGA'	=>	'2,1',
	'GTTAGC'	=>	'2,1',
	'GTTAGG'	=>	'2,1',
	'GTTAGT'	=>	'2,0',
	'GTTATA'	=>	'1,1',
	'GTTATC'	=>	'1,1',
	'GTTATG'	=>	'1,1',
	'GTTATT'	=>	'1,0',
	'GTTCAA'	=>	'2,1',
	'GTTCAC'	=>	'2,1',
	'GTTCAG'	=>	'2,1',
	'GTTCAT'	=>	'2,0',
	'GTTCCA'	=>	'2,1',
	'GTTCCC'	=>	'2,1',
	'GTTCCG'	=>	'2,1',
	'GTTCCT'	=>	'2,0',
	'GTTCGA'	=>	'2,1',
	'GTTCGC'	=>	'2,1',
	'GTTCGG'	=>	'2,1',
	'GTTCGT'	=>	'2,0',
	'GTTCTA'	=>	'1,1',
	'GTTCTC'	=>	'1,1',
	'GTTCTG'	=>	'1,1',
	'GTTCTT'	=>	'1,0',
	'GTTGAA'	=>	'1,1',
	'GTTGAC'	=>	'1,1',
	'GTTGAG'	=>	'1,1',
	'GTTGAT'	=>	'1,0',
	'GTTGCA'	=>	'1,1',
	'GTTGCC'	=>	'1,1',
	'GTTGCG'	=>	'1,1',
	'GTTGCT'	=>	'1,0',
	'GTTGGA'	=>	'1,1',
	'GTTGGC'	=>	'1,1',
	'GTTGGG'	=>	'1,1',
	'GTTGGT'	=>	'1,0',
	'GTTGTA'	=>	'0,1',
	'GTTGTC'	=>	'0,1',
	'GTTGTG'	=>	'0,1',
	'GTTGTT'	=>	'0,0',
	'GTTTAC'	=>	'2,1',
	'GTTTAT'	=>	'2,0',
	'GTTTCA'	=>	'2,1',
	'GTTTCC'	=>	'2,1',
	'GTTTCG'	=>	'2,1',
	'GTTTCT'	=>	'2,0',
	'GTTTGC'	=>	'2,1',
	'GTTTGG'	=>	'2,1',
	'GTTTGT'	=>	'2,0',
	'GTTTTA'	=>	'1,1',
	'GTTTTC'	=>	'1,1',
	'GTTTTG'	=>	'1,1',
	'GTTTTT'	=>	'1,0',
	'TACAAA'	=>	'2,0',
	'TACAAC'	=>	'1,0',
	'TACAAG'	=>	'2,0',
	'TACAAT'	=>	'1,1',
	'TACACA'	=>	'2,1',
	'TACACC'	=>	'2,0',
	'TACACG'	=>	'2,1',
	'TACACT'	=>	'2,1',
	'TACAGA'	=>	'3,0',
	'TACAGC'	=>	'2,0',
	'TACAGG'	=>	'3,0',
	'TACAGT'	=>	'2,1',
	'TACATA'	=>	'2,1',
	'TACATC'	=>	'2,0',
	'TACATG'	=>	'3,0',
	'TACATT'	=>	'2,1',
	'TACCAA'	=>	'2,0',
	'TACCAC'	=>	'1,0',
	'TACCAG'	=>	'2,0',
	'TACCAT'	=>	'1,1',
	'TACCCA'	=>	'2,1',
	'TACCCC'	=>	'2,0',
	'TACCCG'	=>	'2,1',
	'TACCCT'	=>	'2,1',
	'TACCGA'	=>	'2,1',
	'TACCGC'	=>	'2,0',
	'TACCGG'	=>	'2,1',
	'TACCGT'	=>	'2,1',
	'TACCTA'	=>	'2,1',
	'TACCTC'	=>	'2,0',
	'TACCTG'	=>	'2,1',
	'TACCTT'	=>	'2,1',
	'TACGAA'	=>	'2,0',
	'TACGAC'	=>	'1,0',
	'TACGAG'	=>	'2,0',
	'TACGAT'	=>	'1,1',
	'TACGCA'	=>	'2,1',
	'TACGCC'	=>	'2,0',
	'TACGCG'	=>	'2,1',
	'TACGCT'	=>	'2,1',
	'TACGGA'	=>	'2,1',
	'TACGGC'	=>	'2,0',
	'TACGGG'	=>	'2,1',
	'TACGGT'	=>	'2,1',
	'TACGTA'	=>	'2,1',
	'TACGTC'	=>	'2,0',
	'TACGTG'	=>	'2,1',
	'TACGTT'	=>	'2,1',
	'TACTAC'	=>	'0,0',
	'TACTAT'	=>	'0,1',
	'TACTCA'	=>	'1,1',
	'TACTCC'	=>	'1,0',
	'TACTCG'	=>	'1,1',
	'TACTCT'	=>	'1,1',
	'TACTGC'	=>	'1,0',
	'TACTGG'	=>	'2,0',
	'TACTGT'	=>	'1,1',
	'TACTTA'	=>	'2,0',
	'TACTTC'	=>	'1,0',
	'TACTTG'	=>	'2,0',
	'TACTTT'	=>	'1,1',
	'TATAAA'	=>	'2,0',
	'TATAAC'	=>	'1,1',
	'TATAAG'	=>	'2,0',
	'TATAAT'	=>	'1,0',
	'TATACA'	=>	'2,1',
	'TATACC'	=>	'2,1',
	'TATACG'	=>	'2,1',
	'TATACT'	=>	'2,0',
	'TATAGA'	=>	'3,0',
	'TATAGC'	=>	'2,1',
	'TATAGG'	=>	'3,0',
	'TATAGT'	=>	'2,0',
	'TATATA'	=>	'2,1',
	'TATATC'	=>	'2,1',
	'TATATG'	=>	'3,0',
	'TATATT'	=>	'2,0',
	'TATCAA'	=>	'2,0',
	'TATCAC'	=>	'1,1',
	'TATCAG'	=>	'2,0',
	'TATCAT'	=>	'1,0',
	'TATCCA'	=>	'2,1',
	'TATCCC'	=>	'2,1',
	'TATCCG'	=>	'2,1',
	'TATCCT'	=>	'2,0',
	'TATCGA'	=>	'2,1',
	'TATCGC'	=>	'2,1',
	'TATCGG'	=>	'2,1',
	'TATCGT'	=>	'2,0',
	'TATCTA'	=>	'2,1',
	'TATCTC'	=>	'2,1',
	'TATCTG'	=>	'2,1',
	'TATCTT'	=>	'2,0',
	'TATGAA'	=>	'2,0',
	'TATGAC'	=>	'1,1',
	'TATGAG'	=>	'2,0',
	'TATGAT'	=>	'1,0',
	'TATGCA'	=>	'2,1',
	'TATGCC'	=>	'2,1',
	'TATGCG'	=>	'2,1',
	'TATGCT'	=>	'2,0',
	'TATGGA'	=>	'2,1',
	'TATGGC'	=>	'2,1',
	'TATGGG'	=>	'2,1',
	'TATGGT'	=>	'2,0',
	'TATGTA'	=>	'2,1',
	'TATGTC'	=>	'2,1',
	'TATGTG'	=>	'2,1',
	'TATGTT'	=>	'2,0',
	'TATTAC'	=>	'0,1',
	'TATTAT'	=>	'0,0',
	'TATTCA'	=>	'1,1',
	'TATTCC'	=>	'1,1',
	'TATTCG'	=>	'1,1',
	'TATTCT'	=>	'1,0',
	'TATTGC'	=>	'1,1',
	'TATTGG'	=>	'2,0',
	'TATTGT'	=>	'1,0',
	'TATTTA'	=>	'2,0',
	'TATTTC'	=>	'1,1',
	'TATTTG'	=>	'2,0',
	'TATTTT'	=>	'1,0',
	'TCAAAA'	=>	'2,0',
	'TCAAAC'	=>	'2,1',
	'TCAAAG'	=>	'2,1',
	'TCAAAT'	=>	'2,1',
	'TCAACA'	=>	'1,0',
	'TCAACC'	=>	'1,1',
	'TCAACG'	=>	'1,1',
	'TCAACT'	=>	'1,1',
	'TCAAGA'	=>	'2,0',
	'TCAAGC'	=>	'2,1',
	'TCAAGG'	=>	'2,1',
	'TCAAGT'	=>	'2,1',
	'TCAATA'	=>	'2,0',
	'TCAATC'	=>	'2,1',
	'TCAATG'	=>	'2,1',
	'TCAATT'	=>	'2,1',
	'TCACAA'	=>	'2,0',
	'TCACAC'	=>	'2,1',
	'TCACAG'	=>	'2,1',
	'TCACAT'	=>	'2,1',
	'TCACCA'	=>	'1,0',
	'TCACCC'	=>	'1,1',
	'TCACCG'	=>	'1,1',
	'TCACCT'	=>	'1,1',
	'TCACGA'	=>	'2,0',
	'TCACGC'	=>	'2,1',
	'TCACGG'	=>	'2,1',
	'TCACGT'	=>	'2,1',
	'TCACTA'	=>	'1,1',
	'TCACTC'	=>	'1,2',
	'TCACTG'	=>	'1,2',
	'TCACTT'	=>	'1,2',
	'TCAGAA'	=>	'2,0',
	'TCAGAC'	=>	'2,1',
	'TCAGAG'	=>	'2,1',
	'TCAGAT'	=>	'2,1',
	'TCAGCA'	=>	'1,0',
	'TCAGCC'	=>	'1,1',
	'TCAGCG'	=>	'1,1',
	'TCAGCT'	=>	'1,1',
	'TCAGGA'	=>	'2,0',
	'TCAGGC'	=>	'2,1',
	'TCAGGG'	=>	'2,1',
	'TCAGGT'	=>	'2,1',
	'TCAGTA'	=>	'2,0',
	'TCAGTC'	=>	'2,1',
	'TCAGTG'	=>	'2,1',
	'TCAGTT'	=>	'2,1',
	'TCATAC'	=>	'1,1',
	'TCATAT'	=>	'1,1',
	'TCATCA'	=>	'0,0',
	'TCATCC'	=>	'0,1',
	'TCATCG'	=>	'0,1',
	'TCATCT'	=>	'0,1',
	'TCATGC'	=>	'1,1',
	'TCATGG'	=>	'1,1',
	'TCATGT'	=>	'1,1',
	'TCATTA'	=>	'1,0',
	'TCATTC'	=>	'1,1',
	'TCATTG'	=>	'1,1',
	'TCATTT'	=>	'1,1',
	'TCCAAA'	=>	'2,1',
	'TCCAAC'	=>	'2,0',
	'TCCAAG'	=>	'2,1',
	'TCCAAT'	=>	'2,1',
	'TCCACA'	=>	'1,1',
	'TCCACC'	=>	'1,0',
	'TCCACG'	=>	'1,1',
	'TCCACT'	=>	'1,1',
	'TCCAGA'	=>	'2,1',
	'TCCAGC'	=>	'2,0',
	'TCCAGG'	=>	'2,1',
	'TCCAGT'	=>	'2,1',
	'TCCATA'	=>	'2,1',
	'TCCATC'	=>	'2,0',
	'TCCATG'	=>	'2,1',
	'TCCATT'	=>	'2,1',
	'TCCCAA'	=>	'2,1',
	'TCCCAC'	=>	'2,0',
	'TCCCAG'	=>	'2,1',
	'TCCCAT'	=>	'2,1',
	'TCCCCA'	=>	'1,1',
	'TCCCCC'	=>	'1,0',
	'TCCCCG'	=>	'1,1',
	'TCCCCT'	=>	'1,1',
	'TCCCGA'	=>	'2,1',
	'TCCCGC'	=>	'2,0',
	'TCCCGG'	=>	'2,1',
	'TCCCGT'	=>	'2,1',
	'TCCCTA'	=>	'1,2',
	'TCCCTC'	=>	'2,0',
	'TCCCTG'	=>	'1,2',
	'TCCCTT'	=>	'2,1',
	'TCCGAA'	=>	'2,1',
	'TCCGAC'	=>	'2,0',
	'TCCGAG'	=>	'2,1',
	'TCCGAT'	=>	'2,1',
	'TCCGCA'	=>	'1,1',
	'TCCGCC'	=>	'1,0',
	'TCCGCG'	=>	'1,1',
	'TCCGCT'	=>	'1,1',
	'TCCGGA'	=>	'2,1',
	'TCCGGC'	=>	'2,0',
	'TCCGGG'	=>	'2,1',
	'TCCGGT'	=>	'2,1',
	'TCCGTA'	=>	'2,1',
	'TCCGTC'	=>	'2,0',
	'TCCGTG'	=>	'2,1',
	'TCCGTT'	=>	'2,1',
	'TCCTAC'	=>	'1,0',
	'TCCTAT'	=>	'1,1',
	'TCCTCA'	=>	'0,1',
	'TCCTCC'	=>	'0,0',
	'TCCTCG'	=>	'0,1',
	'TCCTCT'	=>	'0,1',
	'TCCTGC'	=>	'1,0',
	'TCCTGG'	=>	'1,1',
	'TCCTGT'	=>	'1,1',
	'TCCTTA'	=>	'1,1',
	'TCCTTC'	=>	'1,0',
	'TCCTTG'	=>	'1,1',
	'TCCTTT'	=>	'1,1',
	'TCGAAA'	=>	'2,1',
	'TCGAAC'	=>	'2,1',
	'TCGAAG'	=>	'2,0',
	'TCGAAT'	=>	'2,1',
	'TCGACA'	=>	'1,1',
	'TCGACC'	=>	'1,1',
	'TCGACG'	=>	'1,0',
	'TCGACT'	=>	'1,1',
	'TCGAGA'	=>	'2,1',
	'TCGAGC'	=>	'2,1',
	'TCGAGG'	=>	'2,0',
	'TCGAGT'	=>	'2,1',
	'TCGATA'	=>	'2,1',
	'TCGATC'	=>	'2,1',
	'TCGATG'	=>	'2,0',
	'TCGATT'	=>	'2,1',
	'TCGCAA'	=>	'2,1',
	'TCGCAC'	=>	'2,1',
	'TCGCAG'	=>	'2,0',
	'TCGCAT'	=>	'2,1',
	'TCGCCA'	=>	'1,1',
	'TCGCCC'	=>	'1,1',
	'TCGCCG'	=>	'1,0',
	'TCGCCT'	=>	'1,1',
	'TCGCGA'	=>	'2,1',
	'TCGCGC'	=>	'2,1',
	'TCGCGG'	=>	'2,0',
	'TCGCGT'	=>	'2,1',
	'TCGCTA'	=>	'1,2',
	'TCGCTC'	=>	'1,2',
	'TCGCTG'	=>	'1,1',
	'TCGCTT'	=>	'1,2',
	'TCGGAA'	=>	'2,1',
	'TCGGAC'	=>	'2,1',
	'TCGGAG'	=>	'2,0',
	'TCGGAT'	=>	'2,1',
	'TCGGCA'	=>	'1,1',
	'TCGGCC'	=>	'1,1',
	'TCGGCG'	=>	'1,0',
	'TCGGCT'	=>	'1,1',
	'TCGGGA'	=>	'2,1',
	'TCGGGC'	=>	'2,1',
	'TCGGGG'	=>	'2,0',
	'TCGGGT'	=>	'2,1',
	'TCGGTA'	=>	'2,1',
	'TCGGTC'	=>	'2,1',
	'TCGGTG'	=>	'2,0',
	'TCGGTT'	=>	'2,1',
	'TCGTAC'	=>	'1,1',
	'TCGTAT'	=>	'1,1',
	'TCGTCA'	=>	'0,1',
	'TCGTCC'	=>	'0,1',
	'TCGTCG'	=>	'0,0',
	'TCGTCT'	=>	'0,1',
	'TCGTGC'	=>	'1,1',
	'TCGTGG'	=>	'1,0',
	'TCGTGT'	=>	'1,1',
	'TCGTTA'	=>	'1,1',
	'TCGTTC'	=>	'1,1',
	'TCGTTG'	=>	'1,0',
	'TCGTTT'	=>	'1,1',
	'TCTAAA'	=>	'2,1',
	'TCTAAC'	=>	'2,1',
	'TCTAAG'	=>	'2,1',
	'TCTAAT'	=>	'2,0',
	'TCTACA'	=>	'1,1',
	'TCTACC'	=>	'1,1',
	'TCTACG'	=>	'1,1',
	'TCTACT'	=>	'1,0',
	'TCTAGA'	=>	'2,1',
	'TCTAGC'	=>	'2,1',
	'TCTAGG'	=>	'2,1',
	'TCTAGT'	=>	'2,0',
	'TCTATA'	=>	'2,1',
	'TCTATC'	=>	'2,1',
	'TCTATG'	=>	'2,1',
	'TCTATT'	=>	'2,0',
	'TCTCAA'	=>	'2,1',
	'TCTCAC'	=>	'2,1',
	'TCTCAG'	=>	'2,1',
	'TCTCAT'	=>	'2,0',
	'TCTCCA'	=>	'1,1',
	'TCTCCC'	=>	'1,1',
	'TCTCCG'	=>	'1,1',
	'TCTCCT'	=>	'1,0',
	'TCTCGA'	=>	'2,1',
	'TCTCGC'	=>	'2,1',
	'TCTCGG'	=>	'2,1',
	'TCTCGT'	=>	'2,0',
	'TCTCTA'	=>	'1,2',
	'TCTCTC'	=>	'2,1',
	'TCTCTG'	=>	'1,2',
	'TCTCTT'	=>	'2,0',
	'TCTGAA'	=>	'2,1',
	'TCTGAC'	=>	'2,1',
	'TCTGAG'	=>	'2,1',
	'TCTGAT'	=>	'2,0',
	'TCTGCA'	=>	'1,1',
	'TCTGCC'	=>	'1,1',
	'TCTGCG'	=>	'1,1',
	'TCTGCT'	=>	'1,0',
	'TCTGGA'	=>	'2,1',
	'TCTGGC'	=>	'2,1',
	'TCTGGG'	=>	'2,1',
	'TCTGGT'	=>	'2,0',
	'TCTGTA'	=>	'2,1',
	'TCTGTC'	=>	'2,1',
	'TCTGTG'	=>	'2,1',
	'TCTGTT'	=>	'2,0',
	'TCTTAC'	=>	'1,1',
	'TCTTAT'	=>	'1,0',
	'TCTTCA'	=>	'0,1',
	'TCTTCC'	=>	'0,1',
	'TCTTCG'	=>	'0,1',
	'TCTTCT'	=>	'0,0',
	'TCTTGC'	=>	'1,1',
	'TCTTGG'	=>	'1,1',
	'TCTTGT'	=>	'1,0',
	'TCTTTA'	=>	'1,1',
	'TCTTTC'	=>	'1,1',
	'TCTTTG'	=>	'1,1',
	'TCTTTT'	=>	'1,0',
	'TGCAAA'	=>	'3,0',
	'TGCAAC'	=>	'2,0',
	'TGCAAG'	=>	'3,0',
	'TGCAAT'	=>	'2,1',
	'TGCACA'	=>	'2,1',
	'TGCACC'	=>	'2,0',
	'TGCACG'	=>	'2,1',
	'TGCACT'	=>	'2,1',
	'TGCAGA'	=>	'2,0',
	'TGCAGC'	=>	'1,0',
	'TGCAGG'	=>	'2,0',
	'TGCAGT'	=>	'1,1',
	'TGCATA'	=>	'2,1',
	'TGCATC'	=>	'2,0',
	'TGCATG'	=>	'3,0',
	'TGCATT'	=>	'2,1',
	'TGCCAA'	=>	'2,1',
	'TGCCAC'	=>	'2,0',
	'TGCCAG'	=>	'2,1',
	'TGCCAT'	=>	'2,1',
	'TGCCCA'	=>	'2,1',
	'TGCCCC'	=>	'2,0',
	'TGCCCG'	=>	'2,1',
	'TGCCCT'	=>	'2,1',
	'TGCCGA'	=>	'1,1',
	'TGCCGC'	=>	'1,0',
	'TGCCGG'	=>	'1,1',
	'TGCCGT'	=>	'1,1',
	'TGCCTA'	=>	'2,1',
	'TGCCTC'	=>	'2,0',
	'TGCCTG'	=>	'2,1',
	'TGCCTT'	=>	'2,1',
	'TGCGAA'	=>	'2,1',
	'TGCGAC'	=>	'2,0',
	'TGCGAG'	=>	'2,1',
	'TGCGAT'	=>	'2,1',
	'TGCGCA'	=>	'2,1',
	'TGCGCC'	=>	'2,0',
	'TGCGCG'	=>	'2,1',
	'TGCGCT'	=>	'2,1',
	'TGCGGA'	=>	'1,1',
	'TGCGGC'	=>	'1,0',
	'TGCGGG'	=>	'1,1',
	'TGCGGT'	=>	'1,1',
	'TGCGTA'	=>	'2,1',
	'TGCGTC'	=>	'2,0',
	'TGCGTG'	=>	'2,1',
	'TGCGTT'	=>	'2,1',
	'TGCTAC'	=>	'1,0',
	'TGCTAT'	=>	'1,1',
	'TGCTCA'	=>	'1,1',
	'TGCTCC'	=>	'1,0',
	'TGCTCG'	=>	'1,1',
	'TGCTCT'	=>	'1,1',
	'TGCTGC'	=>	'0,0',
	'TGCTGG'	=>	'1,0',
	'TGCTGT'	=>	'0,1',
	'TGCTTA'	=>	'2,0',
	'TGCTTC'	=>	'1,0',
	'TGCTTG'	=>	'2,0',
	'TGCTTT'	=>	'1,1',
	'TGGAAA'	=>	'2,1',
	'TGGAAC'	=>	'3,0',
	'TGGAAG'	=>	'2,0',
	'TGGAAT'	=>	'3,0',
	'TGGACA'	=>	'2,1',
	'TGGACC'	=>	'2,1',
	'TGGACG'	=>	'2,0',
	'TGGACT'	=>	'2,1',
	'TGGAGA'	=>	'1,1',
	'TGGAGC'	=>	'2,0',
	'TGGAGG'	=>	'1,0',
	'TGGAGT'	=>	'2,0',
	'TGGATA'	=>	'2,1',
	'TGGATC'	=>	'3,0',
	'TGGATG'	=>	'2,0',
	'TGGATT'	=>	'3,0',
	'TGGCAA'	=>	'2,1',
	'TGGCAC'	=>	'2,1',
	'TGGCAG'	=>	'2,0',
	'TGGCAT'	=>	'2,1',
	'TGGCCA'	=>	'2,1',
	'TGGCCC'	=>	'2,1',
	'TGGCCG'	=>	'2,0',
	'TGGCCT'	=>	'2,1',
	'TGGCGA'	=>	'1,1',
	'TGGCGC'	=>	'1,1',
	'TGGCGG'	=>	'1,0',
	'TGGCGT'	=>	'1,1',
	'TGGCTA'	=>	'1,2',
	'TGGCTC'	=>	'1,2',
	'TGGCTG'	=>	'1,1',
	'TGGCTT'	=>	'1,2',
	'TGGGAA'	=>	'2,1',
	'TGGGAC'	=>	'2,1',
	'TGGGAG'	=>	'2,0',
	'TGGGAT'	=>	'2,1',
	'TGGGCA'	=>	'2,1',
	'TGGGCC'	=>	'2,1',
	'TGGGCG'	=>	'2,0',
	'TGGGCT'	=>	'2,1',
	'TGGGGA'	=>	'1,1',
	'TGGGGC'	=>	'1,1',
	'TGGGGG'	=>	'1,0',
	'TGGGGT'	=>	'1,1',
	'TGGGTA'	=>	'2,1',
	'TGGGTC'	=>	'2,1',
	'TGGGTG'	=>	'2,0',
	'TGGGTT'	=>	'2,1',
	'TGGTAC'	=>	'2,0',
	'TGGTAT'	=>	'2,0',
	'TGGTCA'	=>	'1,1',
	'TGGTCC'	=>	'1,1',
	'TGGTCG'	=>	'1,0',
	'TGGTCT'	=>	'1,1',
	'TGGTGC'	=>	'1,0',
	'TGGTGG'	=>	'0,0',
	'TGGTGT'	=>	'1,0',
	'TGGTTA'	=>	'1,1',
	'TGGTTC'	=>	'2,0',
	'TGGTTG'	=>	'1,0',
	'TGGTTT'	=>	'2,0',
	'TGTAAA'	=>	'3,0',
	'TGTAAC'	=>	'2,1',
	'TGTAAG'	=>	'3,0',
	'TGTAAT'	=>	'2,0',
	'TGTACA'	=>	'2,1',
	'TGTACC'	=>	'2,1',
	'TGTACG'	=>	'2,1',
	'TGTACT'	=>	'2,0',
	'TGTAGA'	=>	'2,0',
	'TGTAGC'	=>	'1,1',
	'TGTAGG'	=>	'2,0',
	'TGTAGT'	=>	'1,0',
	'TGTATA'	=>	'2,1',
	'TGTATC'	=>	'2,1',
	'TGTATG'	=>	'3,0',
	'TGTATT'	=>	'2,0',
	'TGTCAA'	=>	'2,1',
	'TGTCAC'	=>	'2,1',
	'TGTCAG'	=>	'2,1',
	'TGTCAT'	=>	'2,0',
	'TGTCCA'	=>	'2,1',
	'TGTCCC'	=>	'2,1',
	'TGTCCG'	=>	'2,1',
	'TGTCCT'	=>	'2,0',
	'TGTCGA'	=>	'1,1',
	'TGTCGC'	=>	'1,1',
	'TGTCGG'	=>	'1,1',
	'TGTCGT'	=>	'1,0',
	'TGTCTA'	=>	'2,1',
	'TGTCTC'	=>	'2,1',
	'TGTCTG'	=>	'2,1',
	'TGTCTT'	=>	'2,0',
	'TGTGAA'	=>	'2,1',
	'TGTGAC'	=>	'2,1',
	'TGTGAG'	=>	'2,1',
	'TGTGAT'	=>	'2,0',
	'TGTGCA'	=>	'2,1',
	'TGTGCC'	=>	'2,1',
	'TGTGCG'	=>	'2,1',
	'TGTGCT'	=>	'2,0',
	'TGTGGA'	=>	'1,1',
	'TGTGGC'	=>	'1,1',
	'TGTGGG'	=>	'1,1',
	'TGTGGT'	=>	'1,0',
	'TGTGTA'	=>	'2,1',
	'TGTGTC'	=>	'2,1',
	'TGTGTG'	=>	'2,1',
	'TGTGTT'	=>	'2,0',
	'TGTTAC'	=>	'1,1',
	'TGTTAT'	=>	'1,0',
	'TGTTCA'	=>	'1,1',
	'TGTTCC'	=>	'1,1',
	'TGTTCG'	=>	'1,1',
	'TGTTCT'	=>	'1,0',
	'TGTTGC'	=>	'0,1',
	'TGTTGG'	=>	'1,0',
	'TGTTGT'	=>	'0,0',
	'TGTTTA'	=>	'2,0',
	'TGTTTC'	=>	'1,1',
	'TGTTTG'	=>	'2,0',
	'TGTTTT'	=>	'1,0',
	'TTAAAA'	=>	'2,0',
	'TTAAAC'	=>	'2,1',
	'TTAAAG'	=>	'2,1',
	'TTAAAT'	=>	'2,1',
	'TTAACA'	=>	'2,0',
	'TTAACC'	=>	'2,1',
	'TTAACG'	=>	'2,1',
	'TTAACT'	=>	'2,1',
	'TTAAGA'	=>	'2,0',
	'TTAAGC'	=>	'2,1',
	'TTAAGG'	=>	'2,1',
	'TTAAGT'	=>	'2,1',
	'TTAATA'	=>	'1,0',
	'TTAATC'	=>	'1,1',
	'TTAATG'	=>	'1,1',
	'TTAATT'	=>	'1,1',
	'TTACAA'	=>	'1,1',
	'TTACAC'	=>	'1,2',
	'TTACAG'	=>	'1,2',
	'TTACAT'	=>	'1,2',
	'TTACCA'	=>	'1,1',
	'TTACCC'	=>	'1,2',
	'TTACCG'	=>	'1,2',
	'TTACCT'	=>	'1,2',
	'TTACGA'	=>	'1,1',
	'TTACGC'	=>	'1,2',
	'TTACGG'	=>	'1,2',
	'TTACGT'	=>	'1,2',
	'TTACTA'	=>	'0,1',
	'TTACTC'	=>	'0,2',
	'TTACTG'	=>	'0,2',
	'TTACTT'	=>	'0,2',
	'TTAGAA'	=>	'2,0',
	'TTAGAC'	=>	'2,1',
	'TTAGAG'	=>	'2,1',
	'TTAGAT'	=>	'2,1',
	'TTAGCA'	=>	'2,0',
	'TTAGCC'	=>	'2,1',
	'TTAGCG'	=>	'2,1',
	'TTAGCT'	=>	'2,1',
	'TTAGGA'	=>	'2,0',
	'TTAGGC'	=>	'2,1',
	'TTAGGG'	=>	'2,1',
	'TTAGGT'	=>	'2,1',
	'TTAGTA'	=>	'1,0',
	'TTAGTC'	=>	'1,1',
	'TTAGTG'	=>	'1,1',
	'TTAGTT'	=>	'1,1',
	'TTATAC'	=>	'2,0',
	'TTATAT'	=>	'2,0',
	'TTATCA'	=>	'1,0',
	'TTATCC'	=>	'1,1',
	'TTATCG'	=>	'1,1',
	'TTATCT'	=>	'1,1',
	'TTATGC'	=>	'2,0',
	'TTATGG'	=>	'1,1',
	'TTATGT'	=>	'2,0',
	'TTATTA'	=>	'0,0',
	'TTATTC'	=>	'1,0',
	'TTATTG'	=>	'0,1',
	'TTATTT'	=>	'1,0',
	'TTCAAA'	=>	'2,1',
	'TTCAAC'	=>	'2,0',
	'TTCAAG'	=>	'3,0',
	'TTCAAT'	=>	'2,1',
	'TTCACA'	=>	'2,1',
	'TTCACC'	=>	'2,0',
	'TTCACG'	=>	'2,1',
	'TTCACT'	=>	'2,1',
	'TTCAGA'	=>	'2,1',
	'TTCAGC'	=>	'2,0',
	'TTCAGG'	=>	'3,0',
	'TTCAGT'	=>	'2,1',
	'TTCATA'	=>	'1,1',
	'TTCATC'	=>	'1,0',
	'TTCATG'	=>	'2,0',
	'TTCATT'	=>	'1,1',
	'TTCCAA'	=>	'2,1',
	'TTCCAC'	=>	'2,0',
	'TTCCAG'	=>	'2,1',
	'TTCCAT'	=>	'2,1',
	'TTCCCA'	=>	'2,1',
	'TTCCCC'	=>	'2,0',
	'TTCCCG'	=>	'2,1',
	'TTCCCT'	=>	'2,1',
	'TTCCGA'	=>	'2,1',
	'TTCCGC'	=>	'2,0',
	'TTCCGG'	=>	'2,1',
	'TTCCGT'	=>	'2,1',
	'TTCCTA'	=>	'1,1',
	'TTCCTC'	=>	'1,0',
	'TTCCTG'	=>	'1,1',
	'TTCCTT'	=>	'1,1',
	'TTCGAA'	=>	'2,1',
	'TTCGAC'	=>	'2,0',
	'TTCGAG'	=>	'2,1',
	'TTCGAT'	=>	'2,1',
	'TTCGCA'	=>	'2,1',
	'TTCGCC'	=>	'2,0',
	'TTCGCG'	=>	'2,1',
	'TTCGCT'	=>	'2,1',
	'TTCGGA'	=>	'2,1',
	'TTCGGC'	=>	'2,0',
	'TTCGGG'	=>	'2,1',
	'TTCGGT'	=>	'2,1',
	'TTCGTA'	=>	'1,1',
	'TTCGTC'	=>	'1,0',
	'TTCGTG'	=>	'1,1',
	'TTCGTT'	=>	'1,1',
	'TTCTAC'	=>	'1,0',
	'TTCTAT'	=>	'1,1',
	'TTCTCA'	=>	'1,1',
	'TTCTCC'	=>	'1,0',
	'TTCTCG'	=>	'1,1',
	'TTCTCT'	=>	'1,1',
	'TTCTGC'	=>	'1,0',
	'TTCTGG'	=>	'2,0',
	'TTCTGT'	=>	'1,1',
	'TTCTTA'	=>	'1,0',
	'TTCTTC'	=>	'0,0',
	'TTCTTG'	=>	'1,0',
	'TTCTTT'	=>	'0,1',
	'TTGAAA'	=>	'2,1',
	'TTGAAC'	=>	'3,0',
	'TTGAAG'	=>	'2,0',
	'TTGAAT'	=>	'3,0',
	'TTGACA'	=>	'2,1',
	'TTGACC'	=>	'2,1',
	'TTGACG'	=>	'2,0',
	'TTGACT'	=>	'2,1',
	'TTGAGA'	=>	'2,1',
	'TTGAGC'	=>	'3,0',
	'TTGAGG'	=>	'2,0',
	'TTGAGT'	=>	'3,0',
	'TTGATA'	=>	'1,1',
	'TTGATC'	=>	'2,0',
	'TTGATG'	=>	'1,0',
	'TTGATT'	=>	'2,0',
	'TTGCAA'	=>	'1,2',
	'TTGCAC'	=>	'1,2',
	'TTGCAG'	=>	'1,1',
	'TTGCAT'	=>	'1,2',
	'TTGCCA'	=>	'1,2',
	'TTGCCC'	=>	'1,2',
	'TTGCCG'	=>	'1,1',
	'TTGCCT'	=>	'1,2',
	'TTGCGA'	=>	'1,2',
	'TTGCGC'	=>	'1,2',
	'TTGCGG'	=>	'1,1',
	'TTGCGT'	=>	'1,2',
	'TTGCTA'	=>	'0,2',
	'TTGCTC'	=>	'0,2',
	'TTGCTG'	=>	'0,1',
	'TTGCTT'	=>	'0,2',
	'TTGGAA'	=>	'2,1',
	'TTGGAC'	=>	'2,1',
	'TTGGAG'	=>	'2,0',
	'TTGGAT'	=>	'2,1',
	'TTGGCA'	=>	'2,1',
	'TTGGCC'	=>	'2,1',
	'TTGGCG'	=>	'2,0',
	'TTGGCT'	=>	'2,1',
	'TTGGGA'	=>	'2,1',
	'TTGGGC'	=>	'2,1',
	'TTGGGG'	=>	'2,0',
	'TTGGGT'	=>	'2,1',
	'TTGGTA'	=>	'1,1',
	'TTGGTC'	=>	'1,1',
	'TTGGTG'	=>	'1,0',
	'TTGGTT'	=>	'1,1',
	'TTGTAC'	=>	'2,0',
	'TTGTAT'	=>	'2,0',
	'TTGTCA'	=>	'1,1',
	'TTGTCC'	=>	'1,1',
	'TTGTCG'	=>	'1,0',
	'TTGTCT'	=>	'1,1',
	'TTGTGC'	=>	'2,0',
	'TTGTGG'	=>	'1,0',
	'TTGTGT'	=>	'2,0',
	'TTGTTA'	=>	'0,1',
	'TTGTTC'	=>	'1,0',
	'TTGTTG'	=>	'0,0',
	'TTGTTT'	=>	'1,0',
	'TTTAAA'	=>	'2,1',
	'TTTAAC'	=>	'2,1',
	'TTTAAG'	=>	'3,0',
	'TTTAAT'	=>	'2,0',
	'TTTACA'	=>	'2,1',
	'TTTACC'	=>	'2,1',
	'TTTACG'	=>	'2,1',
	'TTTACT'	=>	'2,0',
	'TTTAGA'	=>	'2,1',
	'TTTAGC'	=>	'2,1',
	'TTTAGG'	=>	'3,0',
	'TTTAGT'	=>	'2,0',
	'TTTATA'	=>	'1,1',
	'TTTATC'	=>	'1,1',
	'TTTATG'	=>	'2,0',
	'TTTATT'	=>	'1,0',
	'TTTCAA'	=>	'2,1',
	'TTTCAC'	=>	'2,1',
	'TTTCAG'	=>	'2,1',
	'TTTCAT'	=>	'2,0',
	'TTTCCA'	=>	'2,1',
	'TTTCCC'	=>	'2,1',
	'TTTCCG'	=>	'2,1',
	'TTTCCT'	=>	'2,0',
	'TTTCGA'	=>	'2,1',
	'TTTCGC'	=>	'2,1',
	'TTTCGG'	=>	'2,1',
	'TTTCGT'	=>	'2,0',
	'TTTCTA'	=>	'1,1',
	'TTTCTC'	=>	'1,1',
	'TTTCTG'	=>	'1,1',
	'TTTCTT'	=>	'1,0',
	'TTTGAA'	=>	'2,1',
	'TTTGAC'	=>	'2,1',
	'TTTGAG'	=>	'2,1',
	'TTTGAT'	=>	'2,0',
	'TTTGCA'	=>	'2,1',
	'TTTGCC'	=>	'2,1',
	'TTTGCG'	=>	'2,1',
	'TTTGCT'	=>	'2,0',
	'TTTGGA'	=>	'2,1',
	'TTTGGC'	=>	'2,1',
	'TTTGGG'	=>	'2,1',
	'TTTGGT'	=>	'2,0',
	'TTTGTA'	=>	'1,1',
	'TTTGTC'	=>	'1,1',
	'TTTGTG'	=>	'1,1',
	'TTTGTT'	=>	'1,0',
	'TTTTAC'	=>	'1,1',
	'TTTTAT'	=>	'1,0',
	'TTTTCA'	=>	'1,1',
	'TTTTCC'	=>	'1,1',
	'TTTTCG'	=>	'1,1',
	'TTTTCT'	=>	'1,0',
	'TTTTGC'	=>	'1,1',
	'TTTTGG'	=>	'2,0',
	'TTTTGT'	=>	'1,0',
	'TTTTTA'	=>	'1,0',
	'TTTTTC'	=>	'0,1',
	'TTTTTG'	=>	'1,0',
	'TTTTTT'	=>	'0,0',

    );
return \%m;
}
