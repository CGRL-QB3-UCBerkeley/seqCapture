#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;
use Getopt::Std;
use Getopt::Long;
use List::Util qw[min max];
use List::Util qw(sum);
use List::MoreUtils qw(part);
use List::MoreUtils qw/ uniq /;
use File::Temp;
use List::Util 'shuffle';


die(qq/

Usage: seqCapture generateoutfiles [options]

Basic options:

-aln             DIR                  [FULL PATH] A folder with all individual alignments (locus1.aln, 
                                      locus2.aln, locus3.aln...)

-codingPos       DIR                  [FULL PATH] A folder with all coding end pisitions (locus1_coding_end_position.txt, 
                                      locus2_coding_end_position.txt, locus3_coding_end_position.txt...) [null]

-resDir          DIR                  [FULL PATH] Results folder

-out             INT                  Options for output files (only run one option each run!!!!!): 
                                      1 - phylip for raxml
                                      2 - SNAPP
                                      3 - Splitstree
                                      4 - STRUCTURE
                                      5 - SmartPCA
                                      6 - GPhocs
                                      7 - BPP3  ## need to run TransExonCapPhylo Combine if exon+flanking is desired
                                      8 - Nexus ## need to run -out 1 first!                                                 
                                      9 - produce raxml alignment for each marker 
                                      10 - TESS                                                                               
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
additional options:
-informative                          For out = 2-5,10        only keep informative sites [null]
-random                               For out = 2-5,10        random pick one snp per locus [null]
-start          INT                   For out = 1             user defined start position for partition file [0]
-codon                                For out = 1             if you want to parition the data by 1st+2nd,3rd [null] 
-codingcodon                          For out = 1             if you want to parition the data by 1st+2nd,3rd and noncoding [null]
-coding                               For out = 1             if you want to parition the data by coding\/noncoding [null]
-numSample      INT                   For out = 6             number of alignments for GPhocs; 
                                                              if numSample = 0, use all alignments [0]  
-subsamples     FILE                  For out = 1-5,7,10      a subset of samples to be included in the analyses (format see blow).
                                                              Notice that this file has to be modified based on "sampleID.txt" produced by
                                                              TransExonCapPhylo Alignment!!! [null]
-toploci        INT                   For out = 1-5,7,10      If used 0 then all loci are inculded in the downstream analyses. If INT (>0)
                                                              then only keep the INT high quality loci (ranked by percent missing data) 
                                                              NOTE: -toploci and -topmissing should be used in combination![0]
                                                              [ranking loci by length first and pick those containing least % missing data]
-topmissing     FLOAT                 For out = 1-5,7,10      A threshold (0 - 1) defines ratio of missing data in any loci in individuals that will 
                                                              be counted as low quality loci. For example, if 0.6 is given it means that for a 
                                                              given alignment, if any individual sequence contains more than 60% missing data (Ns) 
                                                              then this entire alignment is removed from the dataset 
                                                              If given 0 then no missing data is allowed in any alignments.
                                                              If given 1 then this filter is turned off, all alignments will be used for toploci selection
                                                              NOTE: -toploci and -topmissing should be used in combination![1]
-topsubloci     INT                   For out = 1             divide the number of alignments sampled in toploci by INT to generate sub-datasets. 
                                                              For example if you sample 260 loci in toploci and use -topsubloci 5, it will generate 
                                                              5 phylip and partition files each with 52 loci (randomly subsampling without replacement). 
                                                              0 means no subsampling [0] 
-blocks         DIR                   For out = 8             a folder with raxml outputs. If you have both exons (with or without flanking) 
                                                              and UCEs, name them as exons_Partition.txt, exons.phylip, uces_Partition.txt,
                                                              uces.phylip and put them in the same folder. Make sure use -start when 
                                                              generating UCE raxml inputs.
                                                              if you only have exons  (with or without flanking), just name the files as 
                                                              exons_Partition.txt and exons.phylip and put them in the same folder. The same
                                                              applies to uces markers only 
-totlen         INT                   For out = 8             total length of the alignment. You can get this number from the last line of Partition.txt.
                                                              if both exons and uces are used, then you can get this number from the last line of 
                                                              uces_Partition.txt (make sure that -start is used when you generate raxml input for UCE)
-popInfo        FILE                  For out = 4             Population information (sampleID popID). popID must be numerical
                                                              Note: if the -subsamples (above) is used then the popInfo must be modified 
                                                              based on the subsamples file, not the "sampleID.txt"!!! [null]
-missing        FLOAT                 For out = 2-5,10        missing data allowed to keep a SNP [0.2]
-het            FLOAT                 For out = 2-5,10        maximum proportion of shared 
                                                              polymorphic sites in a locus   [0.2]

Format for subsamples: must based on the order in sampleID.txt (1 means to include the samples; 0 means to exclude):
sample1 0
sample2 1
sample3 1
sample4 0
sample5 1
sample6 1
sample7 0
sample8 1

\n\n/) unless (@ARGV);
  my %opts = (missing => [0.2], het=> [0.2], numSample => [0], start => [0], toploci => [0], topmissing=>[1], topsubloci=>[0]);
  my ($aln, $out,  $res, $snp, $geno, $resDir, $popInfo,$codingPos, $subsamples, $blocks, $totlen ) = (undef, undef, undef, undef,undef,undef,undef, undef,undef,undef, undef);
  my $informative;
  my $random;
  my $nondi;
  my $codon;
  my $codingcodon;
  my $coding;
  

  GetOptions('out=s@{1,}' => \$out, 'aln=s@{1,1}' => \$aln, 'snp=s@{1,1}' => \$snp, 'geno=s@{1,1}' => \$geno, 'popInfo=s@{1,1}' => \$popInfo, 'resDir=s@{1,1}' => \$resDir, 'codingPos=s@{1,1}' => \$codingPos, 'informative'  => \$informative, 'random' => \$random, 'missing=s@{1,1}' => \$opts{missing}, 'het=s@{1,1}' => \$opts{het}, 'coding'  => \$coding,  'codon'  => \$codon, 'codingcodon' => \$codingcodon,  'start=s@{1,1}' => \$opts{start}, 'numSample=s@{1,1}' => \$opts{numSample}, 'subsamples=s@{1,1}' => \$subsamples, 'blocks=s@{1,1}' => \$blocks, 'totlen=s@{1,1}' => \$totlen , 'toploci=s@{1,1}' => \$opts{toploci}, 'topmissing=s@{1,1}' => \$opts{topmissing}, 'topsubloci=s@{1,1}' => \$opts{topsubloci});
    
  my $alndir = dir ($aln);
  my $snpdir = $1 . "Individual_SNPs/" if $alndir =~ /(\S+)Individual_ALNs/; 
  my $genodir = $1 . "Individual_GENOs/" if $alndir =~ /(\S+)Individual_ALNs/; 

  my $resdir = dir ($resDir);
  my $codingdir = dir ($codingPos) if $codingPos;

  my $blockdir = dir ($blocks) if $blocks;
  
  mkdir $resdir unless -e $resdir;
  my $topmissing = $opts{topmissing}[0];
 
  
  my $topsubloci = $opts{topsubloci}[0];
  
  my $species = 0;   
  my $samplelist = $alndir . "sampleID.txt";
  open (my $list, "<", $samplelist);
  while (<$list>) {
    $species ++ unless $_ =~ /^$/; 
  }
  close $list;
  
  my $mis = $opts{missing}[0];
  my $hets = $opts{het}[0];
  my $pop = @{$popInfo}[0];
  my $numsam = $opts{numSample}[0];
  my $start = $opts{start}[0];
  my $ln = @{$totlen}[0] if $totlen;
  my $toploci = $opts{toploci}[0];
  
### if only focus on a subset of the samples for structure and SNAPP...
  my $subfile =  @{$subsamples}[0] if  @{$subsamples}[0];


  if ($subfile) {
   
    my %sub;
    my $include = 1;
    my $string = "-f";
    
    my $subsnp = $snpdir . "subsample/";
    mkdir $subsnp unless -e $subsnp;

    my $subgeno = $genodir . "subsample/";
    mkdir $subgeno unless -e $subgeno;

    my $subaln = $alndir . "subsample/";
    mkdir $subaln unless -e $subaln;
    
    
    my $newsampleID = $subsnp . "sampleID.txt";
    my $newsampleID2 = $subgeno . "sampleID.txt";
    my $newsampleID3 = $subaln . "sampleID.txt";
    
    open (INFILE, "<", $subfile);
    open (OUTFILE, ">", $newsampleID);
    
    while (<INFILE>) {
      chomp (my @line = split/\s+/,$_);
      if ($line[1] == 1) {
	$sub{$include}++;
	print OUTFILE $line[0], "\n"; 
      }
      $include ++;
    }
    close INFILE;
    close OUTFILE;
     
    system ("cp $newsampleID  $newsampleID2");
    system ("cp $newsampleID  $newsampleID3");
    
    my $dds = 1;  
    foreach my $id (sort {$a <=> $b} keys %sub) {
      $string .= $id. "," if $dds < scalar keys %sub;
      $string .= $id if $dds == scalar keys %sub;
      $dds++;
    }

    my @SNP = <$snpdir*_SNP>;

    foreach (@SNP) {
      my $file = $_;
      my $new = $subsnp . basename($file);     
      system ("cut $string $file > $new");
    }
   
    my @GENO = <$genodir*_geno>;
    foreach (@GENO) {
      my $file = $_;
      my $new = $subgeno . basename($file);     
      system ("cut $string $file > $new");
    }


    my %miss;
     open (IN2, "<", $subfile);
      while (<IN2>) {
      chomp (my @line = split/\s+/,$_);
      if ($line[1] == 0) {
	$miss{$line[0]}++;
      }
    }
    close IN2;
    

    my @ALN = <$alndir*.aln>;
    foreach (@ALN) {
      my $file = $_;
      my $new = $subaln . basename($file);
      open (IN, "<", $file);
      open (OUT, ">", $new);

      while (<IN>) {
	chomp (my $line = $_);
	if ($line =~ m/^>(\S+)/) {
    
	  my $id = $1;
	 
	  if ($miss{$id}) {
	    chomp (my $seq = <IN>);
	  }
	  else {
	    print OUT ">", $id, "\n";
	    chomp (my $seq1 = <IN>);
	    print OUT $seq1, "\n";
	  }
	  
	}
	
      }
      close IN;
      close OUT;
    }

   
    
    $snpdir = $subsnp;
    $genodir = $subgeno;
    $alndir = $subaln;

    
    $species = scalar keys %sub;
    #print $species, " samples are included in the analyses!\n";
    #print $snpdir , "\n";
    #exit;
    #my $nondir = $1 . 'Individual_Non_diallelic/' if $res =~ /(\S+)Individual_SNPs/; ##
  }
  

    if ($toploci > 0) {

      
      my $subsnp = $snpdir . "topsample/";
      system ("rm -r $subsnp") if -e $subsnp;
      mkdir $subsnp unless -e $subsnp;
      
      my $subgeno = $genodir . "topsample/";
      system ("rm -r $subgeno") if -e $subgeno;
      mkdir $subgeno unless -e $subgeno;
      
      my $subaln = $alndir . "topsample/";
      system ("rm -r $subaln") if -e $subaln;
      mkdir $subaln unless -e $subaln;
    
   
  
  
    my %totalmissing;
    my @ALN = <$alndir*.aln>;
    foreach (@ALN) {
      
      my $file = $_;       
      my $locus = $1 if basename($file) =~ /(Contig\d+)\.aln/; 
      my $countd = 0;

      open (IN, "<", $file);
      while (<IN>) {
	chomp (my $line = $_);
	if ($line =~ m/^>(\S+)/) {
	
	  $countd++;
	  my $id = $1;
	  chomp (my $seq = <IN>);
	  my $seq1 = $seq;
	  
	  
	  my $missing = ($seq1 =~ s/[-Nn]//ig);
	  
	  $missing = 0 if ! $missing;
	  my $length = length $seq;

	  my $efflen = $length - $missing;
	  
	  my $ratio = $missing/$length;
	  #print $line, "\t", $ratio, "\n";
          
	  delete $totalmissing{$locus} if $ratio > $topmissing;
	  
	  last if $ratio > $topmissing;

	  $totalmissing{$locus}{'ratio'} += $ratio;
	  $totalmissing{$locus}{'len'} += $efflen;
	}	
      }
    
      
      $totalmissing{$locus}{'ratio'} = $totalmissing{$locus}{'ratio'}/$countd if $totalmissing{$locus};
      $totalmissing{$locus}{'len'} = $totalmissing{$locus}{'len'}/$countd if $totalmissing{$locus};
      
      close IN;
    }

    
    my %loci1;
    my $loci_to_keep1;
    foreach my $locus (sort {$totalmissing{$b}{'len'} <=> $totalmissing{$a}{'len'}} keys %totalmissing) {
      $loci_to_keep1++;
      if  ($loci_to_keep1 <= $toploci + 200) {
	$loci1{$locus}{'ratio'} = $totalmissing{$locus}{'ratio'};
      }	
    }
    
    
    my $loci_to_keep2;
    my %loci2;
    foreach my $locus (sort {$loci1{$a}{'ratio'} <=> $loci1{$b}{'ratio'}} keys %loci1) {
      $loci_to_keep2++;
      if  ($loci_to_keep2 <= $toploci) {
	$loci2{$locus}++;
      }	
    }
    
    @ALN = <$alndir*.aln>;
    foreach (@ALN) {
      my $file = $_;
      my $locus = $1 if basename($file) =~ /(Contig\d+)\.aln/;
      if ($loci2{$locus}) {
	system ("cp $file  $subaln");
      }
    }
    
     my @SNP = <$snpdir*_SNP>;
     foreach (@SNP) {
      my $file = $_;
      my $locus = $1 if basename($file) =~ /(Contig\d+)_/;
      if ($loci2{$locus}) {
	system ("cp $file  $subsnp");
      }
    }
   
    my @GENO = <$genodir*_geno>;
    foreach (@GENO) {
      my $file = $_;
      my $locus = $1 if basename($file) =~ /(Contig\d+)_/;
      if ($loci2{$locus}) {
	system ("cp $file  $subgeno");
      }
    }
        
    my $newsampleID = $snpdir . "sampleID.txt";
    system ("cp $newsampleID  $subgeno");
    system ("cp $newsampleID  $subsnp");
    system ("cp $newsampleID  $subaln");
    
    $snpdir = $subsnp;
    $genodir = $subgeno;
    $alndir = $subaln;
  
  }

  
  
  print "\nNow producing output files for your various fancy analyses to blow everyone's minds !!!!! \n";

  
  foreach my $k (@{$out}) {


    if ($k == 9) {
      print "\nNow it is producing raxml alignment for each marker!\n";
      
      IndAln ($alndir, $species, $resdir,"0", "0", $start) if (!$coding && !$codon && !$codingcodon);
      IndAln ($alndir, $species, $resdir,"1",$codingdir,$start) if ($coding && !$codon && !$codingcodon);
      IndAln ($alndir, $species, $resdir,"2", "0",$start) if (!$coding && $codon && !$codingcodon);
      IndAln ($alndir, $species, $resdir,"3", $codingdir,$start) if (!$coding && !$codon && $codingcodon);

      print "\n";
      print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
      print "\n";
    }
     
    if ($k == 8) {
      print "\nNow it is producing Nexus files!\n";
      Bayes ($blockdir,  $species, $resdir, $ln);
      print "\n";
      print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
      print "\n";
    }
     if ($k == 10) {
      print "\nNow it is producing nexus input files for TESS!\n";
      if ($random && $informative) {
	print "\n\nYou are keeping one random informative SNP per locus for TESS!\n\n";
	tess ($snpdir, $species, $mis, "ri", $hets, $resdir);
	
      }
      if ($random && !$informative) {
	print "\n\nYou are keeping one random SNP per locus for TESS even it is non-informative!\n\n";
	tess ($snpdir, $species, $mis, "r",$hets, $resdir);
      }
      if (!$random && $informative) {
	print "\n\nYou are keeping all informative snps for TESS!\n\n";
	tess ($snpdir, $species, $mis,"i", $hets, $resdir);
      }
      if (!$random && !$informative) {
	print "\n\nYou are keeping all snps for TESS?!\n\n";
	tess ($snpdir, $species, $mis,"a",$hets, $resdir);
      }	 
    }
    
    
     if ($k == 6) {
      print "\nNow it is producing input files for GPhocs!\n";       print "All samples are included!\n";
      GPhocs ($alndir, $resdir, $numsam, $species);
      print "\n";
      print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
      print "\n";
    }

     if ($k == 7) {
      print "\nNow it is producing input files for BPP3!\n";       
      BPP3 ($alndir, $species, $resdir);
      print "\n";
      print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
      print "\n";
    }
    if ($k == 1) {
      print "\nNow it is producing phylip input files for RaxMl!\n";    
      raxml ($alndir, $species, $resdir,"0", "0", $start,$topsubloci) if (!$coding && !$codon && !$codingcodon);
      raxml ($alndir, $species, $resdir,"1",$codingdir,$start,$topsubloci) if ($coding && !$codon && !$codingcodon);
      raxml ($alndir, $species, $resdir,"2", "0",$start,$topsubloci) if (!$coding && $codon && !$codingcodon);
      raxml ($alndir, $species, $resdir,"3", $codingdir,$start,$topsubloci) if (!$coding && !$codon && $codingcodon);
    }
    if ($k == 2) {
      print "\nNow it is producing nexus input files for SNAPP!\n";
      if ($random && $informative) {
	print "\n\nYou are keeping one random informative SNP per locus for SNAPP!\n\n";
	SNAPP ($snpdir, $species, $mis, "ri", $hets, $resdir);
	
      }
      if ($random && !$informative) {
	print "\n\nYou are keeping one random SNP per locus for SNAPP even it is non-informative!\n\n";
	SNAPP ($snpdir, $species, $mis, "r",$hets, $resdir);
      }
      if (!$random && $informative) {
	print "\n\nYou are keeping all informative snps for SNAPP!\n\n";
	SNAPP ($snpdir, $species, $mis,"i", $hets, $resdir);
      }
      if (!$random && !$informative) {
	print "\n\nYou are keeping all snps for SNAPP?!\n\n";
	SNAPP ($snpdir, $species, $mis,"a",$hets, $resdir);
      }	 
    }
    if ($k == 3) {
      print "\nNow it is producing input files for Splits tree!\n";
      if ($random && $informative) {
	print "\n\nYou are keeping one random informative SNP per locus for Splits tree!\n\n";
	Adegenet ($snpdir, $species, $mis, "ri", $hets, $resdir);
      }
      if ($random && !$informative) {
	print "\n\nYou are keeping one random SNP per locus for Splits tree even it is non-informative!\n\n";
	Adegenet ($snpdir, $species, $mis, "r",$hets, $resdir);
      }
      if (!$random && $informative) {
	print "\n\nYou are keeping all informative snps for Splits tree!\n\n";
	Adegenet ($snpdir, $species, $mis,"i", $hets, $resdir); 
      }
      if (!$random && !$informative) {
	print "\n\nYou are keeping all snps for Splits tree?!\n\n";
	Adegenet ($snpdir, $species, $mis,"a",$hets, $resdir);
      }
    }
    if ($k == 4) {
        print "\nNow it is producing input files for Structure!\n";
      if ($random && $informative) {
	print "\n\nYou are keeping one random informative SNP per locus for Structure!\n\n";
	Structure ($snpdir,$genodir, $species, $mis, "ri", $hets, $resdir, $pop);
      }
      if ($random && !$informative) {
	print "\n\nYou are keeping one random SNP per locus for Structure even it is non-informative!\n\n";
	Structure ($snpdir, $genodir, $species, $mis, "r",$hets, $resdir,$pop);
      }
      if (!$random && $informative) {
	print "\n\nYou are keeping all informative snps for Structure!\n\n";
	Structure ($snpdir, $genodir, $species, $mis,"i", $hets, $resdir,$pop);
      }
      if (!$random && !$informative) {
	print "\n\nYou are keeping all snps for Structure?!\n\n";
	Structure ($snpdir, $genodir, $species, $mis,"a",$hets, $resdir,$pop);
      }
    }
    if ($k == 5) {
      print "\nNow it is producing input files for smartPCA!\n";
       if ($random && $informative) {
	print "\n\nYou are keeping one random informative SNP per locus for smartPCA!\n\n";
	PCA ($snpdir, $species, $mis, "ri", $hets, $resdir);
      }
      if ($random && !$informative) {
	print "\n\nYou are keeping one random SNP per locus for smartPCA even it is non-informative!\n\n";
	PCA ($snpdir, $species, $mis, "r",$hets, $resdir);
      }
      if (!$random && $informative) {
	print "\n\nYou are keeping all informative snps for smartPCA!\n\n";
	PCA ($snpdir, $species, $mis,"i", $hets, $resdir);
      }
      if (!$random && !$informative) {
	print "\n\nYou are keeping all snps for smartPCA?!\n\n";
	PCA ($snpdir, $species, $mis,"a",$hets, $resdir);
      }
    }
  }
  


sub tess  {
  my ($res, $number, $missing, $flag, $het, $res2) = @_;
  
  my $resdir = $res2 . "TESSInput/";
  mkdir $resdir unless -e $resdir;

  my $nondir = $1 . 'Individual_Non_diallelic/' if $res =~ /(\S+)Individual/; ##
  my %non;
  my @all = <$nondir*Contig*_Non_diallelic_SNPID.txt>;
  foreach (@all) {
    my $file = $_;
    my $lib = $1 if basename ($file) =~ /(Contig\d+)_Non_diallelic_SNPID.txt/;
    open (NON, "<", $file);
    while (<NON>) {
      chomp (my @line = split /\s+/,$_);
      $non{$lib}{$line[0]}++;
    }
    close NON;  
  }

  
  my @SNP = <$res*_SNP>;
  my $sample = $res. "sampleID.txt";
  my $out =  $resdir . "Final_TESS";
  $out = $out . "_random_informative" if $flag eq 'ri';
  $out = $out . "_random" if $flag eq 'r';
  $out = $out . "_informative" if $flag eq 'i';
  $out = $out . "_all" if $flag eq 'a';

  
  my $file = $resdir . "Final_geno2.txt";
  open (OUT, ">", $file);
  
  foreach (@SNP) {
    my $file = $_;
    my $locus = $1 if basename($file) =~ /(Contig\S+)_SNP/;  
    open (IN, "<", $file);
    my $d = 1;  
    while (<IN>) {
      chomp (my $line = $_);
      print OUT $locus, "\t", $d, "\t", $line, "\n";
      $d++;
    }
    close IN;
  }
  close OUT;
 
 
  open (IN, "<", $file);
  my %hash;
  my @ind;  
  while (<IN>) {
    chomp (my @line = split /\s+/,$_);
    my $geno =join('',@line[2 .. $#line]);
    $geno =~ s/-1/9/g;
    my @miss = ($geno =~ m/9/g);
    my @het = ($geno =~ m/1/g);
    my @homR =  ($geno =~ m/0/g);
    my @homA =  ($geno =~ m/2/g);
    if ($flag eq 'ri' || $flag eq "i") {  
      if (scalar(@homA) + scalar(@het) > 1 && scalar(@homR) > 0 ) {
	if (scalar(@het)/$number <= $het && scalar(@miss)/$number <= $missing && !$non{$line[0]}{$line[1]}) {
	    $hash{$line[0]}{$line[1]} = $geno;	  
	}
      }
    }
    
    if ($flag eq 'r' || $flag eq "a") {
      if ( scalar(@homR) > 0 ) {
	if (scalar(@homA) > 0 || scalar(@het) > 0) {
	  if (scalar(@het)/$number <= $het && scalar(@miss)/$number <= $missing && !$non{$line[0]}{$line[1]} ) {
	    $hash{$line[0]}{$line[1]} = $geno;
	  }
	}
      } 
    }
  }
  close IN;
 
  open (ID, "<", $sample);  
  while (<ID>) {
    chomp (my $line = $_);
    push (@ind, $line);
  }
  close ID;
  
  if ($flag eq 'r' || $flag eq 'ri') {
    my $selected = $resdir . "selected_one_SNP_per_Contig_for_TESS.txt";
   
    my $selected1 = random (\%hash, $resdir, $selected);
    system ("cut -f3 $selected1 > $out ");
    
  }
  if ($flag eq 'a' || $flag eq 'i') {
    my $selected1 = $resdir . "selected_ALL_SNP_chosen_for_TESS.txt";
    open (O, ">", $selected1);
    for my $locus (sort {$a cmp $b} keys %hash){
      for my $pos (sort {$a <=> $b} keys %{$hash{$locus}}){
	print O $locus, "\t", $pos, "\t", $hash{$locus}{$pos},"\n";
      }
    }
    close O;
    system ("cut -f3 $selected1 > $out ");
    
  }
  
}  



sub IndAln {
  my ($dir, $s, $res2, $tag, $coding, $start) = @_;

  my $resdir = $res2 . "IndAln/";
  mkdir $resdir unless -e $resdir;
 

  if ($tag == 3) {
  my @aln1 = <$dir*.aln>;
  my @site =  <$coding*_coding_end_position.txt>; 

  my %pos;
  foreach (@site) {
    my $file = $_;
    
    my $locus = $1 if basename($file) =~ /(Contig\S+)_coding_end_position.txt/;
    open (IN, "<", $file);
    while (<IN>) {
      chomp (my $line = $_);
      $pos{$locus} = $line;   
      
    }
    close IN;
  } 
 
 
  foreach (@aln1) {
    
    my %aln;  
    my $file = $_;
  
    my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
    open (IN, "<", $file);
    my $ids;
    my $sample; 
    while (<IN>) {
      chomp (my $line = $_);    
      if ($line =~ m/^>(\S+)/ ) {
	$ids = $1;
	
	chomp (my $seq = <IN>);
	if ($seq =~ /A|T|G|C/) {
	  $sample++;
	  $aln{$ids}{$locus} = $seq;
	}
      }
    }
    close IN;
  
    ##raxml ouput:
    my $out = $resdir . $locus. ".phylip";	
    my $par =  $resdir . $locus.  "_Partition.txt";
 

  my $all;
  open (OUT, ">", $out);
  open (PAR, ">", $par);
  
  my $c;
  $c = $start-1 if $start > 0;
  $c = 0 if $start == 0;
  foreach my $s (sort {$a cmp $b} keys %aln) {
    foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
      if ($pos{$contig} > 3) {
      
      print PAR  "DNA,", $contig,"c12", "=";
      print PAR $c+1, "-", $c+$pos{$contig}, "\\","3,", $c+2, "-", $c+$pos{$contig}, "\\","3", "\n";
       
      print PAR  "DNA,", $contig,"c3", "=";
      print PAR $c+3, "-", $c+$pos{$contig}, "\\3","\n";

      $c += $pos{$contig};
     if ( $pos{$contig} < length ($aln{$s}{$contig}) ) {

      print PAR  "DNA,", $contig, "_f", "=";
      
      print PAR $c+1, "-", $c + length ($aln{$s}{$contig})  - $pos{$contig}, "\n";
     }
     
      $c =  length ($aln{$s}{$contig})  - $pos{$contig} + $c;
     
      $all += length ($aln{$s}{$contig});
      }
     else {
      print PAR  "DNA,", $contig, "_f", "=";
      print PAR $c+1, "-", $c+length ($aln{$s}{$contig}), "\n";
      $c += length ($aln{$s}{$contig});
      $all += length ($aln{$s}{$contig});
     }
    }
   last;
  } ## foreach my $s (sort {$a cmp $b} keys %aln) {
  close PAR;
 # print "\nThe total length of the alignment is: ", $all, "bp!", "\n"; 
  
  print OUT $sample, " ", $all,"\n";
  foreach my $s (sort {$a cmp $b} keys %aln) {
    my $total;    
    print OUT $s, " ";
    #my $name_space = length ($sample) + 1;
    foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
      $total .= $aln{$s}{$contig};
    } 
    print OUT " ", $total, "\n";
  } ##foreach my $s (sort {$a cmp $b} keys %aln) {
  close OUT;
  
  #print "\nthe total number of loci in the nexus matrix is: ", scalar (@aln1), "!", "\n";
  
  #print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   } #foreach (@aln1) { 
 } ## if $tag == 3



  if ($tag == 2) {
  my @aln1 = <$dir*.aln>;

  
  
  foreach (@aln1) {
    my %aln;  
    my $file = $_;
    my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
    open (IN, "<", $file);
    my $ids;
    my $sample;
    while (<IN>) {
      chomp (my $line = $_);    
      if ($line =~ m/^>(\S+)/ ) {
	$ids = $1;
	chomp (my $seq = <IN>);
	if ($seq =~ /A|T|G|C/) {
	  $sample++;
	  $aln{$ids}{$locus} = $seq;
	}
      }
    }
    close IN;
  
    ##raxml ouput:
    my $out = $resdir . $locus. ".phylip";	
    my $par =  $resdir . $locus.  "_Partition.txt";
 
  my $all;
  open (OUT, ">", $out);
  open (PAR, ">", $par);
  
  my $c;
  $c = $start-1 if $start > 0;
  $c = 0 if $start == 0;
  
  foreach my $s (sort {$a cmp $b} keys %aln) {
    foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
      print PAR  "DNA,", $contig,"c12", "=";
      print PAR $c+1, "-", $c+length ($aln{$s}{$contig}), "\\3,", $c+2, "-", $c+length ($aln{$s}{$contig}), "\\3", "\n";
       
      print PAR  "DNA,", $contig,"c3", "=";
      print PAR $c+3, "-", $c+length ($aln{$s}{$contig}), "\\3","\n";
      $c += length ($aln{$s}{$contig});
      $all += length ($aln{$s}{$contig});
    }
    last;
  } ## foreach my $s (sort {$a cmp $b} keys %aln) {
  close PAR;
  #print "\nThe total length of the alignment is: ", $all, "bp!", "\n"; 
  
  print OUT $sample, " ", $all,"\n";
  foreach my $s (sort {$a cmp $b} keys %aln) {
    my $total;    
    print OUT $s, " ";
    #my $name_space = length ($sample) + 1;
    foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
      $total .= $aln{$s}{$contig};
    } 
    print OUT " ", $total, "\n";
  } ##foreach my $s (sort {$a cmp $b} keys %aln) {
  close OUT;
  
  #print "\nthe total number of loci in the nexus matrix is: ", scalar (@aln1), "!", "\n";
  
  #print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  }
 } ## if $tag == 2


  if ($tag == 1) {
  my @aln1 = <$dir*.aln>;
  my @site =  <$coding*_coding_end_position.txt>; 

  my %pos;
    foreach (@site) {
    
    my $file = $_;  
    my $locus = $1 if basename($file) =~ /(Contig\d+)_coding_end_position.txt/;
    open (IN, "<", $file);
    while (<IN>) {
    chomp (my $line = $_);
    $pos{$locus} = $line;   
    
    }
    close IN;
  }
  
  
  foreach (@aln1) {
    my %aln;  
    my $file = $_;
    my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
    open (IN, "<", $file);
    my $ids;
    my $sample;
    while (<IN>) {
      chomp (my $line = $_);    
      if ($line =~ m/^>(\S+)/ ) {
	$ids = $1;
	chomp (my $seq = <IN>);
	if ($seq =~ /A|T|G|C/) {
	  $sample++;
	  $aln{$ids}{$locus} = $seq;
	}
      }
    }
    close IN;
  
    ##raxml ouput:
    my $out = $resdir . $locus. ".phylip";	
    my $par =  $resdir . $locus.  "_Partition.txt";
 

  my $all;
  open (OUT, ">", $out);
  open (PAR, ">", $par);
  
  my $c;
  $c = $start-1 if $start > 0;
  $c = 0 if $start == 0;
  
  foreach my $s (sort {$a cmp $b} keys %aln) {
    foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
    if ($pos{$contig} > 3) {
      print PAR  "DNA,", $contig, "_c", "=";    
      print PAR $c+1, "-", $c+ $pos{$contig}, "\n";     
      $c += $pos{$contig};
 
     if ( $pos{$contig} < length ($aln{$s}{$contig}) ) {

      print PAR  "DNA,", $contig, "_f", "=";
      
      print PAR $c+1, "-", $c + length ($aln{$s}{$contig})  - $pos{$contig}, "\n";
     }
     
      $c =  length ($aln{$s}{$contig})  - $pos{$contig} + $c;
     
      $all += length ($aln{$s}{$contig});
      }
     else {
      print PAR  "DNA,", $contig, "_f", "=";
      print PAR $c+1, "-", $c+length ($aln{$s}{$contig}), "\n";
      $c += length ($aln{$s}{$contig});
      $all += length ($aln{$s}{$contig});
     }



    }
    last;
  } ## foreach my $s (sort {$a cmp $b} keys %aln) {
  close PAR;
 # print "\nThe total length of the alignment is: ", $all, "bp!", "\n"; 
  
  print OUT $sample, " ", $all,"\n";
  foreach my $s (sort {$a cmp $b} keys %aln) {
    my $total;    
    print OUT $s, " ";
    #my $name_space = length ($sample) + 1;
    foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
      $total .= $aln{$s}{$contig};
    } 
    print OUT " ", $total, "\n";
  } ##foreach my $s (sort {$a cmp $b} keys %aln) {
  close OUT;
  
  #print "\nthe total number of loci in the nexus matrix is: ", scalar (@aln1), "!", "\n";
  
  #print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  }
 } ## if $tag == 1


  if ($tag == 0) {
  my @aln1 = <$dir*.aln>;
  
  
  foreach (@aln1) {
    my %aln;  
    my $file = $_;
    my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
    open (IN, "<", $file);
    my $ids;
    my $sample;
    while (<IN>) {
      chomp (my $line = $_);    
      if ($line =~ m/^>(\S+)/ ) {
	$ids = $1;
	chomp (my $seq = <IN>);
	if ($seq =~ /A|T|G|C/) {
	  $sample++;
	  $aln{$ids}{$locus} = $seq;
	}
      }
    }
    close IN;
  
    ##raxml ouput:
    my $out = $resdir . $locus. ".phylip";	
    my $par =  $resdir . $locus.  "_Partition.txt";
 
  my $all;
  open (OUT, ">", $out);
  open (PAR, ">", $par);
  
  my $c;
  $c = $start-1 if $start > 0;
  $c = 0 if $start == 0;
  
  foreach my $s (sort {$a cmp $b} keys %aln) {
    foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
      print PAR  "DNA,", $contig, "=";
      print PAR $c+1, "-", $c+length ($aln{$s}{$contig}), "\n";
      $c += length ($aln{$s}{$contig});
      $all += length ($aln{$s}{$contig});
    }
    last;
  } ## foreach my $s (sort {$a cmp $b} keys %aln) {
  close PAR;
  #print "\nThe total length of the alignment is: ", $all, "bp!", "\n"; 
  
  print OUT $sample, " ", $all,"\n";
  foreach my $s (sort {$a cmp $b} keys %aln) {
    my $total;    
    print OUT $s, " ";
    #my $name_space = length ($sample) + 1;
    foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
      $total .= $aln{$s}{$contig};
    } 
    print OUT " ", $total, "\n";
  } ##foreach my $s (sort {$a cmp $b} keys %aln) {
  close OUT;
  
  #print "\nthe total number of loci in the nexus matrix is: ", scalar (@aln1), "!", "\n";
  
 # print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  }
 } ## if $tag == 0
 
} ##sub IndAln



sub Bayes {
  my ($blockdir,  $species, $res2, $len) = @_;
  
  my $resdir = $res2 . "NexusInput/";
  mkdir $resdir unless -e $resdir;
  my $nexusout =  $resdir . "nexus.input";
  open (OUT, ">", $nexusout);
  print OUT "#NEXUS", "\n";
  print OUT "\n";
  print OUT "BEGIN DATA;\n";
  print OUT "\t\tDIMENSIONS\t";
  print OUT "NTAX=", $species, "\t";
  print OUT "NCHAR=", $len, ";\n";
  print OUT "\t\tFORMAT\tDATATYPE=DNA\tinterleave\tMISSING=N\tGAP=-;", "\n";
  print OUT "MATRIX\n\n";
  

  my @alignment = <$blockdir*phylip>;
  my @partition = <$blockdir*Partition.txt>;


  foreach (@alignment) {
    my $file = $_;
    my $lib = $1 if basename($file) =~ m/(\S+)\.phylip/;
    open (IN, "<", $file);
    my $first = <IN>;
    while (<IN>) {
      chomp (my $line = $_);
      print OUT $line, "\n";
    }
    close IN;
    print OUT "\n";
  }

  print OUT ";","\n";
  print OUT "END;\n\n";
  print OUT "Begin assumptions;", "\n";

  my @contigs;
  my $dd;
  
  foreach (@partition) {
    my $file = $_;
    my $lib = $1 if basename($file) =~ m/(\S+)_Partition\.txt/;
    
    open (IN, "<", $file);
    while (<IN>) {
      chomp (my @line = split /,/, $_);
      my @a = split /=/, $line[1];
      push @contigs, $a[0];
      $dd++;
      print OUT "\t\tcharset\t", $a[0], " = ", $a[1] , ";", "\n";
    }
    close IN;
  }
  print OUT "endblock;\n\n\n";

  print OUT "begin mrbayes;\n";
  @partition = <$blockdir*Partition.txt>;

  foreach (@partition) {
    my $file = $_;
    my $lib = $1 if basename($file) =~ m/(\S+)_Partition\.txt/;
    
    open (IN, "<", $file);
    while (<IN>) {
      chomp (my @line = split /,/, $_);
      my @a = split /=/, $line[1];
      print OUT "\t\tcharset\t", $a[0], " = ", $a[1] , ";", "\n";
    }
    close IN;
  }
  print OUT "\n";
  print OUT " partition by_pos = ", $dd, ":";

  my $d2; 
  foreach (@contigs) {
    $d2++;
    print OUT $_, "," if $d2 < $dd;
    print OUT $_, ";" if $d2 == $dd;

  }
  print OUT "set partition=by_pos;\n";
  for (my $i = 1; $i <= $dd; $i++) {
    print OUT "\t\tlset applyto=($i) nst = 6 rates = invgamma;\n";
  }

  print OUT "\n";
  print OUT "\t\tunlink revmat=(all) shape=(all) pinvar=(all) statefreq=(all);","\n";
  print OUT "\t\tmcmc Nruns=1 Nchains=4 ngen=10000000 printfreq=1000 samplefreq=1000 savebrlens=yes;","\n";
  print OUT "end;","\n";
  

}
  
sub BPP3 {
  my ($alndir, $species, $res2) = @_;
  my $resdir = $res2 . "BPP3Input/";
  mkdir $resdir unless -e $resdir;
  
  my @aln1 = <$alndir*.aln>;
 
  my $bppout =  $resdir . "bpp3.input";
  open (OUT, ">", $bppout);
  
  foreach (@aln1) {
    my $file = $_;
 
    my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
    open (IN, "<", $file);
    my $id = <IN>;
    chomp (my $seq = <IN>);
    my $length = length ($seq);
    seek IN ,0,0;
    print OUT $species, " ", $length, "\n\n";
    while (<IN>) {
      chomp (my $line = $_);    
      if ($line =~ m/^>(\S+)/ ) {
	my $id = $1;
	print OUT "^", $id, "  ";
	chomp (my $seq = <IN>);
	print OUT $seq, "\n\n";
      }
    }
      close IN;
      print OUT "\n";
      
  } ## foreach (@aln1)
  
   close OUT; 
    

} 



sub GPhocs {
  my ($alndir, $res2, $numsam,$species) = @_;
  my $resdir = $res2 . "GPhocsInput/";
  mkdir $resdir unless -e $resdir;

  my $final = $resdir . "GPhocsInput.aln";
  open (FINAL,">", $final);

  
  my @aln1 = <$alndir*.aln>;

  if ($numsam == 0 ) { #use all alignments  
    my $total = scalar (@aln1);
    print FINAL $total, "\n\n";
  }
  else {
    print FINAL $numsam, "\n\n";
  }
  
  foreach (@aln1) {
    my $file = $_;
    my $locus = $1 if basename($file) =~ m/(Contig\d+).aln/;
  
    print FINAL  $locus, "\t", $species, "\t";  
    my $length;
    open (IN, "<", $file);
    while (<IN>) {
      chomp (my $line = $_);
      if ($line !~ m/^>(\S+)/) {
	$length = length ($line);
	last;
      }
      
    }
    seek IN, 0,0;

    print FINAL $length, "\n";
    
    while (<IN>) {
      chomp (my $line = $_);
      if ($line =~ m/^>(\S+)/) {
	print FINAL $1, "\t";
      }
      else {
	my $seq = $line;
	$seq =~ s/-/N/g ;
	print FINAL $seq, "\n";
	
      } 
    }
    close IN;
    print FINAL "\n\n\n"
    
  } ##foreach (@aln1) {
    
  
} #sub GPhocs 

sub raxml {
  my ($dir, $sample, $res2, $tag, $coding, $start,$topsubloci) = @_;

  my $resdir = $res2 . "raxmlInput/";
  mkdir $resdir unless -e $resdir;
  my %sub;
  if ($topsubloci > 0) {
    my @deck = <$dir*.aln>;
    my $num_picks = (scalar @deck) /$topsubloci ;
    
    foreach (my $i = 1; $i <= $topsubloci; $i ++) {
      
      my @shuffled_indexes = shuffle(0..$#deck);
      
      
      my @pick_indexes = @shuffled_indexes[ 0 .. $num_picks - 1 ];  
      my @non = @shuffled_indexes[ $num_picks .. $#deck];  
      # Pick cards from @deck
      my @picks = @deck[ @pick_indexes ];
      @deck = @deck[ @non  ];
      push @{$sub{'a'}}, \@picks;
    }
    
    
    my $d = 0;
    foreach (@{$sub{'a'}}) {
      $d++;
	if ($tag == 3) {
	  my @aln1 = @{$_};
	  my @site =  <$coding*_coding_end_position.txt>; 
	  
	  my %aln;  
	  foreach (@aln1) {
	    my $file = $_;
	    my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
	    open (IN, "<", $file);
	    my $ids;
	    
	    while (<IN>) {
	      chomp (my $line = $_);    
	      if ($line =~ m/^>(\S+)/ ) {
		$ids = $1;	
	      }      
	      else {		
		$aln{$ids}{$locus} .= $line;
	      }
	    }
	    close IN;
	  } ## foreach (@aln1)
	  
	  my %pos;
	  foreach (@site) {
	    my $file = $_;
	    
	    my $locus = $1 if basename($file) =~ /(Contig\S+)_coding_end_position.txt/;
	    open (IN, "<", $file);
	    while (<IN>) {
	      chomp (my $line = $_);
	      $pos{$locus} = $line;   
	      
	    }
	    close IN;
	  } 
	  
	  
	  ##raxml ouput:
	  my $out = $resdir . "Final_alignment_subset$d.phylip";
	  my $par =  $resdir . "Final_alignment_Partition_subset$d.txt";
	  
	  my $all;
	  open (OUT, ">", $out);
	  open (PAR, ">", $par);
	  
	  my $c;
	  $c = $start-1 if $start > 0;
	  $c = 0 if $start == 0;
	  foreach my $s (sort {$a cmp $b} keys %aln) {
	    foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	      if ($pos{$contig} > 3) {
		
		print PAR  "DNA,", $contig,"c12", "=";
		print PAR $c+1, "-", $c+$pos{$contig}, "\\","3,", $c+2, "-", $c+$pos{$contig}, "\\","3", "\n";
		
		print PAR  "DNA,", $contig,"c3", "=";
	    print PAR $c+3, "-", $c+$pos{$contig}, "\\3","\n";
	    
	    $c += $pos{$contig};
	    if ( $pos{$contig} < length ($aln{$s}{$contig}) ) {
	      
	      print PAR  "DNA,", $contig, "_f", "=";
	      
	      print PAR $c+1, "-", $c + length ($aln{$s}{$contig})  - $pos{$contig}, "\n";
	    }
	    
	    $c =  length ($aln{$s}{$contig})  - $pos{$contig} + $c;
	    
	    $all += length ($aln{$s}{$contig});
	  }
	  else {
	    print PAR  "DNA,", $contig, "_f", "=";
	    print PAR $c+1, "-", $c+length ($aln{$s}{$contig}), "\n";
	    $c += length ($aln{$s}{$contig});
	    $all += length ($aln{$s}{$contig});
	  }
	}
	last;
      } ## foreach my $s (sort {$a cmp $b} keys %aln) {
      close PAR;
      print "\nThe total length of the alignment in Final_alignment_subset$d.phylip is: ", $all, "bp!", "\n"; 
      
      print OUT $sample, " ", $all,"\n";
      foreach my $s (sort {$a cmp $b} keys %aln) {
	my $total;    
	print OUT $s, " ";
	#my $name_space = length ($sample) + 1;
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  $total .= $aln{$s}{$contig};
	} 
	print OUT " ", $total, "\n";
      } ##foreach my $s (sort {$a cmp $b} keys %aln) {
      close OUT;
      
      print "\nthe total number of loci in the nexus matrix is: ", scalar (@aln1), "!", "\n";
      
      print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
      
  
    } ## if $tag == 3
    
    
    
    if ($tag == 2) {
      my @aln1 = @{$_};
      
      
      my %aln;  
      foreach (@aln1) {
	my $file = $_;
	my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
	open (IN, "<", $file);
	my $ids;
	
	while (<IN>) {
	  chomp (my $line = $_);    
	  if ($line =~ m/^>(\S+)/ ) {
	    $ids = $1;	
	  }      
	  else {		
	    $aln{$ids}{$locus} .= $line;
	  }
	}
	close IN;
      } ## foreach (@aln1)
      
      
      ##raxml ouput:
      my $out = $resdir . "Final_alignment_subset$d.phylip";
      my $par =  $resdir . "Final_alignment_Partition_subset$d.txt";
      
      my $all;
      open (OUT, ">", $out);
      open (PAR, ">", $par);
      
      my $c;
      $c = $start-1 if $start > 0;
      $c = 0 if $start == 0;
      
      foreach my $s (sort {$a cmp $b} keys %aln) {
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  print PAR  "DNA,", $contig,"c12", "=";
	  print PAR $c+1, "-", $c+length ($aln{$s}{$contig}), "\\3,", $c+2, "-", $c+length ($aln{$s}{$contig}), "\\3", "\n";
	  
	  print PAR  "DNA,", $contig,"c3", "=";
	  print PAR $c+3, "-", $c+length ($aln{$s}{$contig}), "\\3","\n";
	  $c += length ($aln{$s}{$contig});
	  $all += length ($aln{$s}{$contig});
	}
	last;
      } ## foreach my $s (sort {$a cmp $b} keys %aln) {
      close PAR;
      print "\nThe total length of the alignment in Final_alignment_subset$d.phylip is: ", $all, "bp!", "\n"; 
      
      print OUT $sample, " ", $all,"\n";
      foreach my $s (sort {$a cmp $b} keys %aln) {
	my $total;    
	print OUT $s, " ";
	#my $name_space = length ($sample) + 1;
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  $total .= $aln{$s}{$contig};
	} 
	print OUT " ", $total, "\n";
      } ##foreach my $s (sort {$a cmp $b} keys %aln) {
      close OUT;
      
      print "\nthe total number of loci in the nexus matrix is: ", scalar (@aln1), "!", "\n";
      
      print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
      
    } ## if $tag == 2
    
    
    if ($tag == 1) {
      my @aln1 = @{$_};
      my @site =  <$coding*_coding_end_position.txt>; 
      
      my %aln;  
      foreach (@aln1) {
	my $file = $_;
	my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
	open (IN, "<", $file);
	my $ids;
	
	while (<IN>) {
	  chomp (my $line = $_);    
	  if ($line =~ m/^>(\S+)/ ) {
	    $ids = $1;	
	  }      
	  else {		
	    $aln{$ids}{$locus} .= $line;
	    
	  }
	}
	close IN;
      } ## foreach (@aln1)
      
      my %pos;
      foreach (@site) {
	my $file = $_;
	my $locus = $1 if basename($file) =~ /(Contig\S+)_coding_end_position.txt/;
	open (IN, "<", $file);
	while (<IN>) {
	  chomp (my $line = $_);
	  $pos{$locus} = $line;        
	}
	close IN;
      } 
      
      
      
      ##raxml ouput:
      my $out = $resdir . "Final_alignment_subset$d.phylip";
      my $par =  $resdir . "Final_alignment_Partition_subset$d.txt";
      
      my $all;
      open (OUT, ">", $out);
      open (PAR, ">", $par);
      
      my $c;
      $c = $start-1 if $start > 0;
      $c = 0 if $start == 0;
      
      foreach my $s (sort {$a cmp $b} keys %aln) {
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  if ($pos{$contig} > 3) {
	    print PAR  "DNA,", $contig, "_c", "=";    
	    print PAR $c+1, "-", $c+ $pos{$contig}, "\n";     
	    $c += $pos{$contig};
	    
	    if ( $pos{$contig} < length ($aln{$s}{$contig}) ) {
	      
	      print PAR  "DNA,", $contig, "_f", "=";
	      
	      print PAR $c+1, "-", $c + length ($aln{$s}{$contig})  - $pos{$contig}, "\n";
	    }
	    
	    $c =  length ($aln{$s}{$contig})  - $pos{$contig} + $c;
	    
	    $all += length ($aln{$s}{$contig});
	  }
	  else {
	    print PAR  "DNA,", $contig, "_f", "=";
	    print PAR $c+1, "-", $c+length ($aln{$s}{$contig}), "\n";
	    $c += length ($aln{$s}{$contig});
	    $all += length ($aln{$s}{$contig});
	  }
	  
	  
	  
	}
	last;
      } ## foreach my $s (sort {$a cmp $b} keys %aln) {
      close PAR;
      print "\nThe total length of the alignment in Final_alignment_subset$d.phylip is: ", $all, "bp!", "\n"; 
      
      print OUT $sample, " ", $all,"\n";
      foreach my $s (sort {$a cmp $b} keys %aln) {
	my $total;    
	print OUT $s, " ";
	#my $name_space = length ($sample) + 1;
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  $total .= $aln{$s}{$contig};
	} 
	print OUT " ", $total, "\n";
      } ##foreach my $s (sort {$a cmp $b} keys %aln) {
      close OUT;
      
      print "\nthe total number of loci in the nexus matrix is: ", scalar (@aln1), "!", "\n";
      
      print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
      
    } ## if $tag == 1
    
    
    if ($tag == 0) {
      my @aln1 = @{$_};
      
      my %aln;  
      foreach (@aln1) {
	my $file = $_;
	my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
	open (IN, "<", $file);
	my $ids;
	
	while (<IN>) {
	  chomp (my $line = $_);    
	  if ($line =~ m/^>(\S+)/ ) {
	    $ids = $1;	
	  }      
	  else {		
	    $aln{$ids}{$locus} .= $line;
	  }
	}
	close IN;
      } ## foreach (@aln1)
      
      
      ##raxml ouput:
      my $out = $resdir . "Final_alignment_subset$d.phylip";
      my $par =  $resdir . "Final_alignment_Partition_subset$d.txt";
      
      my $all;
      open (OUT, ">", $out);
      open (PAR, ">", $par);
      
      my $c;
      $c = $start-1 if $start > 0;
      $c = 0 if $start == 0;
      
      foreach my $s (sort {$a cmp $b} keys %aln) {
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  print PAR  "DNA,", $contig, "=";
	  print PAR $c+1, "-", $c+length ($aln{$s}{$contig}), "\n";
	  $c += length ($aln{$s}{$contig});
	  $all += length ($aln{$s}{$contig});
	}
	last;
      } ## foreach my $s (sort {$a cmp $b} keys %aln) {
      close PAR;
      print "\nThe total length of the alignment in Final_alignment_subset$d.phylip is: ", $all, "bp!", "\n"; 
      
      print OUT $sample, " ", $all,"\n";
      foreach my $s (sort {$a cmp $b} keys %aln) {
	my $total;    
	print OUT $s, " ";
	#my $name_space = length ($sample) + 1;
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  $total .= $aln{$s}{$contig};
	} 
	print OUT " ", $total, "\n";
      } ##foreach my $s (sort {$a cmp $b} keys %aln) {
      close OUT;
      
      print "\nthe total number of loci in the nexus matrix is: ", scalar (@aln1), "!", "\n";
      
      print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
      
    } ## if $tag == 0
	
	
	
      } #foreach (@{$sub{'a'}})
      
      

    
    
  }
  
  
  if ($topsubloci == 0) {
    if ($tag == 3) {
      my @aln1 = <$dir*.aln>;
      my @site =  <$coding*_coding_end_position.txt>; 
      
      my %aln;  
      foreach (@aln1) {
	my $file = $_;
	my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
	open (IN, "<", $file);
	my $ids;
	
	while (<IN>) {
	  chomp (my $line = $_);    
	  if ($line =~ m/^>(\S+)/ ) {
	    $ids = $1;	
	  }      
	  else {		
	    $aln{$ids}{$locus} .= $line;
	  }
	}
	close IN;
      } ## foreach (@aln1)
      
      my %pos;
      foreach (@site) {
	my $file = $_;
	
	my $locus = $1 if basename($file) =~ /(Contig\S+)_coding_end_position.txt/;
	open (IN, "<", $file);
	while (<IN>) {
	  chomp (my $line = $_);
	  $pos{$locus} = $line;   
	  
	}
	close IN;
      } 
      
      
      ##raxml ouput:
      my $out = $resdir . "Final_alignment.phylip";
      my $par =  $resdir . "Final_alignment_Partition.txt";
      
      my $all;
      open (OUT, ">", $out);
      open (PAR, ">", $par);
      
      my $c;
      $c = $start-1 if $start > 0;
      $c = 0 if $start == 0;
      foreach my $s (sort {$a cmp $b} keys %aln) {
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  if ($pos{$contig} > 3) {
	    
	    print PAR  "DNA,", $contig,"c12", "=";
	    print PAR $c+1, "-", $c+$pos{$contig}, "\\","3,", $c+2, "-", $c+$pos{$contig}, "\\","3", "\n";
	    
	    print PAR  "DNA,", $contig,"c3", "=";
	    print PAR $c+3, "-", $c+$pos{$contig}, "\\3","\n";
	    
	    $c += $pos{$contig};
	    if ( $pos{$contig} < length ($aln{$s}{$contig}) ) {
	      
	      print PAR  "DNA,", $contig, "_f", "=";
	      
	      print PAR $c+1, "-", $c + length ($aln{$s}{$contig})  - $pos{$contig}, "\n";
	    }
	    
	    $c =  length ($aln{$s}{$contig})  - $pos{$contig} + $c;
	    
	    $all += length ($aln{$s}{$contig});
	  }
	  else {
	    print PAR  "DNA,", $contig, "_f", "=";
	    print PAR $c+1, "-", $c+length ($aln{$s}{$contig}), "\n";
	    $c += length ($aln{$s}{$contig});
	    $all += length ($aln{$s}{$contig});
	  }
	}
	last;
      } ## foreach my $s (sort {$a cmp $b} keys %aln) {
      close PAR;
      print "\nThe total length of the alignment is: ", $all, "bp!", "\n"; 
      
      print OUT $sample, " ", $all,"\n";
      foreach my $s (sort {$a cmp $b} keys %aln) {
	my $total;    
	print OUT $s, " ";
	#my $name_space = length ($sample) + 1;
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  $total .= $aln{$s}{$contig};
	} 
	print OUT " ", $total, "\n";
      } ##foreach my $s (sort {$a cmp $b} keys %aln) {
      close OUT;
      
      print "\nthe total number of loci in the nexus matrix is: ", scalar (@aln1), "!", "\n";
      
      print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
      
  
    } ## if $tag == 3
    
    
    
    if ($tag == 2) {
      my @aln1 = <$dir*.aln>;
      
      
      my %aln;  
      foreach (@aln1) {
	my $file = $_;
	my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
	open (IN, "<", $file);
	my $ids;
	
	while (<IN>) {
	  chomp (my $line = $_);    
	  if ($line =~ m/^>(\S+)/ ) {
	    $ids = $1;	
	  }      
	  else {		
	    $aln{$ids}{$locus} .= $line;
	  }
	}
	close IN;
      } ## foreach (@aln1)
      
      
      ##raxml ouput:
      my $out = $resdir . "Final_alignment.phylip";
      my $par =  $resdir . "Final_alignment_Partition.txt";
      
      my $all;
      open (OUT, ">", $out);
      open (PAR, ">", $par);
      
      my $c;
      $c = $start-1 if $start > 0;
      $c = 0 if $start == 0;
      
      foreach my $s (sort {$a cmp $b} keys %aln) {
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  print PAR  "DNA,", $contig,"c12", "=";
	  print PAR $c+1, "-", $c+length ($aln{$s}{$contig}), "\\3,", $c+2, "-", $c+length ($aln{$s}{$contig}), "\\3", "\n";
	  
	  print PAR  "DNA,", $contig,"c3", "=";
	  print PAR $c+3, "-", $c+length ($aln{$s}{$contig}), "\\3","\n";
	  $c += length ($aln{$s}{$contig});
	  $all += length ($aln{$s}{$contig});
	}
	last;
      } ## foreach my $s (sort {$a cmp $b} keys %aln) {
      close PAR;
      print "\nThe total length of the alignment is: ", $all, "bp!", "\n"; 
      
      print OUT $sample, " ", $all,"\n";
      foreach my $s (sort {$a cmp $b} keys %aln) {
	my $total;    
	print OUT $s, " ";
	#my $name_space = length ($sample) + 1;
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  $total .= $aln{$s}{$contig};
	} 
	print OUT " ", $total, "\n";
      } ##foreach my $s (sort {$a cmp $b} keys %aln) {
      close OUT;
      
      print "\nthe total number of loci in the nexus matrix is: ", scalar (@aln1), "!", "\n";
      
      print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
      
    } ## if $tag == 2
    
    
    if ($tag == 1) {
      my @aln1 = <$dir*.aln>;
      my @site =  <$coding*_coding_end_position.txt>; 
      
      my %aln;  
      foreach (@aln1) {
	my $file = $_;
	my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
	open (IN, "<", $file);
	my $ids;
	
	while (<IN>) {
	  chomp (my $line = $_);    
	  if ($line =~ m/^>(\S+)/ ) {
	    $ids = $1;	
	  }      
	  else {		
	    $aln{$ids}{$locus} .= $line;
	    
	  }
	}
	close IN;
      } ## foreach (@aln1)
      
      my %pos;
      foreach (@site) {
	my $file = $_;
	my $locus = $1 if basename($file) =~ /(Contig\S+)_coding_end_position.txt/;
	open (IN, "<", $file);
	while (<IN>) {
	  chomp (my $line = $_);
	  $pos{$locus} = $line;        
	}
	close IN;
      } 
      
      
      
      ##raxml ouput:
      my $out = $resdir . "Final_alignment.phylip";
      my $par =  $resdir . "Final_alignment_Partition.txt";
      
      my $all;
      open (OUT, ">", $out);
      open (PAR, ">", $par);
      
      my $c;
      $c = $start-1 if $start > 0;
      $c = 0 if $start == 0;
      
      foreach my $s (sort {$a cmp $b} keys %aln) {
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  if ($pos{$contig} > 3) {
	    print PAR  "DNA,", $contig, "_c", "=";    
	    print PAR $c+1, "-", $c+ $pos{$contig}, "\n";     
	    $c += $pos{$contig};
	    
	    if ( $pos{$contig} < length ($aln{$s}{$contig}) ) {
	      
	      print PAR  "DNA,", $contig, "_f", "=";
	      
	      print PAR $c+1, "-", $c + length ($aln{$s}{$contig})  - $pos{$contig}, "\n";
	    }
	    
	    $c =  length ($aln{$s}{$contig})  - $pos{$contig} + $c;
	    
	    $all += length ($aln{$s}{$contig});
	  }
	  else {
	    print PAR  "DNA,", $contig, "_f", "=";
	    print PAR $c+1, "-", $c+length ($aln{$s}{$contig}), "\n";
	    $c += length ($aln{$s}{$contig});
	    $all += length ($aln{$s}{$contig});
	  }
	  
	  
	  
	}
	last;
      } ## foreach my $s (sort {$a cmp $b} keys %aln) {
      close PAR;
      print "\nThe total length of the alignment is: ", $all, "bp!", "\n"; 
      
      print OUT $sample, " ", $all,"\n";
      foreach my $s (sort {$a cmp $b} keys %aln) {
	my $total;    
	print OUT $s, " ";
	#my $name_space = length ($sample) + 1;
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  $total .= $aln{$s}{$contig};
	} 
	print OUT " ", $total, "\n";
      } ##foreach my $s (sort {$a cmp $b} keys %aln) {
      close OUT;
      
      print "\nthe total number of loci in the nexus matrix is: ", scalar (@aln1), "!", "\n";
      
      print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
      
    } ## if $tag == 1
    
    
    if ($tag == 0) {
      my @aln1 = <$dir*.aln>;
      
      my %aln;  
      foreach (@aln1) {
	my $file = $_;
	my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
	open (IN, "<", $file);
	my $ids;
	
	while (<IN>) {
	  chomp (my $line = $_);    
	  if ($line =~ m/^>(\S+)/ ) {
	    $ids = $1;	
	  }      
	  else {		
	    $aln{$ids}{$locus} .= $line;
	  }
	}
	close IN;
      } ## foreach (@aln1)
      
      
      ##raxml ouput:
      my $out = $resdir . "Final_alignment.phylip";
      my $par =  $resdir . "Final_alignment_Partition.txt";
      
      my $all;
      open (OUT, ">", $out);
      open (PAR, ">", $par);
      
      my $c;
      $c = $start-1 if $start > 0;
      $c = 0 if $start == 0;
      
      foreach my $s (sort {$a cmp $b} keys %aln) {
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  print PAR  "DNA,", $contig, "=";
	  print PAR $c+1, "-", $c+length ($aln{$s}{$contig}), "\n";
	  $c += length ($aln{$s}{$contig});
	  $all += length ($aln{$s}{$contig});
	}
	last;
      } ## foreach my $s (sort {$a cmp $b} keys %aln) {
      close PAR;
      print "\nThe total length of the alignment is: ", $all, "bp!", "\n"; 
      
      print OUT $sample, " ", $all,"\n";
      foreach my $s (sort {$a cmp $b} keys %aln) {
	my $total;    
	print OUT $s, " ";
	#my $name_space = length ($sample) + 1;
	foreach my $contig (sort {$a cmp $b} keys %{$aln{$s}}) {
	  $total .= $aln{$s}{$contig};
	} 
	print OUT " ", $total, "\n";
      } ##foreach my $s (sort {$a cmp $b} keys %aln) {
      close OUT;
      
      print "\nthe total number of loci in the nexus matrix is: ", scalar (@aln1), "!", "\n";
      
      print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
      
    } ## if $tag == 0
  }  ## toploci == 0
} ##sub raxml






sub Dxy {
  my ($res, $species, $mis, $hets, $res2,$pop) = @_;
  my %pop;
  my $d = 1;
  open (POP, "<", $pop);
  while (<POP>) {
    chomp (my @line = split /\s+/, $_);
    push @{$pop{$line[1]}}, $d;
    $d++;
  }
  close POP;
  
  my $resdir = $res2 . "Dxy/";
  mkdir $resdir unless -e $resdir;
  my $nondir = $1 . 'Individual_Non_diallelic/' if $res =~ /(\S+)Individual_SNPs/; ##
  my %non;
  my $allnon = $resdir . "all_nondi";
  system ("cat $nondir*Non_diallelic_SNPID.txt > $allnon");
  open (NON, "<", $allnon);
  while (<NON>) {
    chomp (my @line = split /\s+/,$_);
    $non{$line[0]}{$line[1]}++;
  }
  close NON;
  unlink ($allnon);

}

sub Adegenet {
  my ($res, $number, $missing, $flag, $het, $res2) = @_;
  # ($snpdir, $species, $mis, "ri", $hets, $resdir);

  my %genos = ('9' => ['N','A'], '0' => ['1','1'], '1' => ['1','2'], '2' => ['2','2']);
   
  my $resdir = $res2 . "AdegenetInput/";
  mkdir $resdir unless -e $resdir;
  
  my $nondir = $1 . 'Individual_Non_diallelic/' if $res =~ /(\S+)Individual_SNPs/; ##
  my %non;
  my $allnon = $resdir . "all_nondi";
  system ("cat $nondir*Non_diallelic_SNPID.txt > $allnon");
  open (NON, "<", $allnon);
  while (<NON>) {
    chomp (my @line = split /\s+/,$_);
    $non{$line[0]}{$line[1]}++;
  }
  close NON;
  unlink ($allnon);
   
  my @SNP = <$res*_SNP>;
  my $sample = $res."sampleID.txt";
  my $out = $resdir . "Final_Adegenet";
  $out = $out . "_random_informative" if $flag eq 'ri';
  $out = $out . "_random" if $flag eq 'r';
  $out = $out . "_informative" if $flag eq 'i';
  $out = $out . "_all" if $flag eq 'a';
  
  my $file = $resdir . "Final_geno2.txt";
  open (OUT, ">", $file);
  
  foreach (@SNP) {
    my $file = $_;
    my $locus = $1 if basename($file) =~ /(Contig\S+)_SNP/;  
    open (IN, "<", $file);
    my $d = 1;  
    while (<IN>) {
      chomp (my $line = $_);
      print OUT $locus, "\t", $d, "\t", $line, "\n";
      $d++;
    }
    close IN;
  }
  close OUT;

  open (IN, "<", $file);
  my %hash;
  my @ind;  
  while (<IN>) {
    chomp (my @line = split /\s+/,$_);
    my $geno =join("\t",@line[2 .. $#line]);
    $geno =~ s/-1/9/g; 
    my @miss = ($geno =~ m/9/g);
    my @het = ($geno =~ m/1/g);
    my @homR =  ($geno =~ m/0/g);
    my @homA =  ($geno =~ m/2/g);
    if ($flag eq 'ri' || $flag eq "i") {    
      if (scalar(@homA) + scalar(@het) > 1 && scalar(@homR) > 0 ) {
	if (scalar(@het)/$number <= $het && scalar(@miss)/$number <= $missing && !$non{$line[0]}{$line[1]}) {
	    $hash{$line[0]}{$line[1]} = $geno;	  
	}
      }
    }
    
    if ($flag eq 'r' || $flag eq "a") {
      if ( scalar(@homR) > 0 ) {
	if (scalar(@homA) > 0 || scalar(@het) > 0) {
	  if (scalar(@het)/$number <= $het && scalar(@miss)/$number <= $missing && !$non{$line[0]}{$line[1]}) {
	    $hash{$line[0]}{$line[1]} = $geno;
	  }
	}
      } 
    }
  }
  close IN;

   if ($flag eq 'r' || $flag eq 'ri') {
    my $selected = $resdir . "selected_one_SNP_per_Contig_for_Adegenet.pos";
    my $selected1 = random (\%hash, $resdir, $selected);  
       
    open (IN2, "<", $selected1);   
    my $n = 0;
    while (<IN2>) {       
      $n++;
    }  
    print "The number of sites selected for Splits Tree analyses is: ", $n, "\n";    
    print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    close IN2;

    my $out1 = $resdir . 'cleaned';
    open (IN, "<", $selected1);
    open (OUT,">", $out1);
    foreach (<IN>) {
      chomp (my @line = split /\s+/, $_);    
      foreach my $m (@line[2..$#line]) {
	print OUT $genos{$m}[0] . $genos{$m}[1], "\t";
      }
      print OUT "\n";		
    }
    close IN;
    close OUT;
    my $out2 = $resdir . 'transposed';    
    transpose($out1, $number, $n, $out2);
    adegenet ($out2 , $n, $out);
    unlink ($out1, $out2);
    
    
  }
  if ($flag eq 'a' || $flag eq 'i') {
    my $selected1 = $resdir . "selected_ALL_SNP_chosen_for_Adegenet.pos";
    
    my $n = 0;
    open (O, ">", $selected1);
 
    for my $locus (sort {$a cmp $b} keys %hash){
      for my $pos (sort {$a <=> $b} keys %{$hash{$locus}}){
	$n++;
	print O $locus, "\t", $pos, "\t", $hash{$locus}{$pos},"\n";
      }
    }
    close O;
    print "The number of sites selected for Splits Tree analyses is: ", $n, "\n";   
    print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    
    
    my $out1 = $resdir . 'cleaned';
    open (IN, "<", $selected1) ;   
    open (OUT,">", $out1);
    foreach (<IN>) {
      chomp (my @line = split /\s+/, $_);    
      foreach my $m (@line[2..$#line]) {
	print OUT $genos{$m}[0] . $genos{$m}[1], "\t";
      }
      print OUT "\n";		
    }
    close IN;
    close OUT;
    
     
    my $out2 = $resdir . 'transposed';
   
    transpose($out1, $number, $n, $out2);
    adegenet ($out2 , $n, $out);
    unlink ($out1, $out2);
  
  }
}




sub transpose {
    my ($file, $size, $site1, $out) = @_;
    
    open (OUT1, '>', $out);  
    for (my $i = 1; $i <= $size; $i++) {
      my $count = $site1;
      open (IN1, '<', $file);
      foreach (<IN1>) {
	chomp();	  
	if ($count > 1) {
	  my @line = split /\s+/, $_;
	  print OUT1 $line[$i-1], "\t";
	  $count--;
	}
	elsif ($count == 1) {
	  my @line = split /\s+/, $_;
	  print OUT1 $line[$i-1], "\n";
	}
	else {
	  last;
	}
      }
    }
    close IN1;
    close OUT1;
   # unlink ($file);
  }

sub adegenet {
  my ($in, $sites, $out) = @_;
  open (IN2, "<", $in);
  open (OPTo, '>', $out);
  my $count = 1;
  for ($count = 1; $count <= $sites; $count++) {
    print OPTo "\t", join ("", "Loc". $count);
  }
  print OPTo "\n";
  
  $count = 1; 
  foreach (<IN2>) {
    print OPTo 'Ind'.$count, "\t", join ("", $_);
    $count++;   
  }  
  #unlink ($in); 
  close OPTo;
  close IN2;
}


 sub PCA {
  my ($res, $number, $missing, $flag, $het, $res2) = @_;
  
  my %genos = ('0' => 0 , '1' => 1, '2'=>2, '-1' => 9, '9' => 9);
  
  my $resdir = $res2 . "smartPCAInput/";
  mkdir $resdir unless -e $resdir;

  my $nondir = $1 . 'Individual_Non_diallelic/' if $res =~ /(\S+)Individual/; ##
  my %non;
  my @all = <$nondir*Contig*_Non_diallelic_SNPID.txt>;
  foreach (@all) {
    my $file = $_;
    my $lib = $1 if basename ($file) =~ /(Contig\d+)_Non_diallelic_SNPID.txt/;
    open (NON, "<", $file);
    while (<NON>) {
      chomp (my @line = split /\s+/,$_);
      $non{$lib}{$line[0]}++;
    }
    close NON;  
  }
 
  my @SNP = <$res*_SNP>;
  my $sample = $res."sampleID.txt";
  my $out = "Final_SNAPP";
  $out = $out . "_random_informative" if $flag eq 'ri';
  $out = $out . "_random" if $flag eq 'r';
  $out = $out . "_informative" if $flag eq 'i';
  $out = $out . "_all" if $flag eq 'a';
  
  my $file = $resdir . "Final_geno2.txt";
  open (OUT, ">", $file);
  
  foreach (@SNP) {
    my $file = $_;
    my $locus = $1 if basename($file) =~ /(Contig\S+)_SNP/;  
    open (IN, "<", $file);
    my $d = 1;  
    while (<IN>) {
      chomp (my $line = $_);
      print OUT $locus, "\t", $d, "\t", $line, "\n";
      $d++;
    }
    close IN;
  }
  close OUT;
 
 
  open (IN, "<", $file);
  my %hash;
  my @ind;  
  while (<IN>) {
    chomp (my @line = split /\s+/,$_);
    my $geno =join("\t",@line[2 .. $#line]);
    $geno =~ s/-1/9/g;
    my @miss = ($geno =~ m/9/g);
    my @het = ($geno =~ m/1/g);
    my @homR =  ($geno =~ m/0/g);
    my @homA =  ($geno =~ m/2/g);
    if ($flag eq 'ri' || $flag eq "i") {  
      if (scalar(@homA) + scalar(@het) > 1 && scalar(@homR) > 0 ) {
	if (scalar(@het)/$number <= $het && scalar(@miss)/$number <= $missing && !$non{$line[0]}{$line[1]}) {
	    $hash{$line[0]}{$line[1]} = $geno;	  
	}
      }
    }
    
    if ($flag eq 'r' || $flag eq "a") {
      if ( scalar(@homR) > 0 ) {
	if (scalar(@homA) > 0 || scalar(@het) > 0) {
	  if (scalar(@het)/$number <= $het && scalar(@miss)/$number <= $missing && !$non{$line[0]}{$line[1]}) {
	    $hash{$line[0]}{$line[1]} = $geno;
	  }
	}
      } 
    }
  }
  close IN;
 
  open (ID, "<", $sample);  
  while (<ID>) {
    chomp (my $line = $_);
    push (@ind, $line);
  }
  close ID;
  
  my $snp;
  my $eigen;
  my $samplefile;
  
  if ($flag eq 'r' || $flag eq 'ri') {
    my $selected = $resdir . "selected_one_SNP_per_Contig_for_smartPCA.pos";
    my $selected1 = random (\%hash, $resdir, $selected);  
    $snp =  $resdir . "selected_one_SNP_per_Contig_for_smartPCA.snp";
    
    open (IN2, "<", $selected1);
    open (OUT, ">", $snp);
    my $n = 0;
    while (<IN2>) {
      chomp (my @line = split /\s+/, $_);
      print OUT $line[0] . "pos". $line[1], "\t", "1", "\t","0.0","\t","1","\n";     
      $n++;
    }
    close IN2;
    close OUT;
    
    print "The number of sites selected for PCA analyses is: ", $n, "\n";
    
    print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";


    my $out1 = 'cleaned';
    trim ($number, $selected1, $out1,\%genos);
    #my $out2 = 'transposed';
    $eigen = $resdir . "selected_one_SNP_per_Contig_for_smartPCA.eigen";
    #transpose($out1, $number, $n, $out2);
    #system ("sed 's/\t//g' $out2 > $eigen");
    #unlink ($out1, $out2);
	system ("sed 's/\t//g' $out1 > $eigen");
            unlink ($out1);


    $samplefile = $resdir . "selected_one_SNP_per_Contig_for_smartPCA.sample";
    open (SAMPLE, ">", $samplefile);
    foreach (@ind) {
      print SAMPLE $_, "\t", "U", "\t", "label","\n";
    }
    close SAMPLE;
    
        
  }
  if ($flag eq 'a' || $flag eq 'i') {
    my $selected1 = $resdir . "selected_ALL_SNP_chosen_for_smartPCA.pos";
    $snp =  $resdir . "selected_ALL_SNP_chosen_for_smartPCA.snp";
    my $n = 0;
    open (O, ">", $selected1);
    open (OUT, ">",$snp );
    for my $locus (sort {$a cmp $b} keys %hash){
      for my $pos (sort {$a <=> $b} keys %{$hash{$locus}}){
	$n++;
	print OUT $locus . "pos". $pos, "\t", "1", "\t","0.0","\t","1","\n";
	print O $locus, "\t", $pos, "\t", $hash{$locus}{$pos},"\n";
      }
    }
    close O;
    print "The number of sites selected for PCA analyses is: ", $n, "\n";
    
    print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    
    
    my $out1 = 'cleaned';
    trim ($number, $selected1, $out1,\%genos);
    #my $out2 = 'transposed';
    $eigen = $resdir . "selected_one_SNP_per_Contig_for_smartPCA.eigen";
    #transpose($out1, $number, $n, $out2);
    #system ("sed 's/\t//g' $out2 > $eigen");
    #unlink ($out1, $out2);
	system ("sed 's/\t//g' $out1 > $eigen");
            unlink ($out1);
    
    $samplefile = $resdir . "selected_ALL_SNP_chosen_for_smartPCA.sample";
    open (SAMPLE, ">", $samplefile);
    foreach (@ind) {
      print SAMPLE $_, "\t", "U", "\t", "label","\n";
    }
    close SAMPLE;
  
  }

  my $par = $resdir . "smartPCA.par";
  my $final  = $resdir . "smartPCA_results";
  open (PAR, ">", $par);
  print PAR "#### input files", "\n";
  print PAR "genotypename: $eigen","\n";
  print PAR "snpname: $snp", "\n";
  print PAR "indivname: $samplefile","\n";

  print PAR "#### output files", "\n";
  print PAR "snpweightoutname: $final",".snpeigs","\n";
  print PAR "evecoutname: $final", ".eigs", "\n";
  print PAR "evaloutname: $final",".eval", "\n";
  print PAR "phylipoutname: $final",".fst","\n";
	  
  print PAR "#### run parameters", "\n";
  print PAR "numoutevec: 20","\n";
  print PAR "numoutlieriter: 0", "\n";
  print PAR "outlieroutname: $final", ".out", "\n";
  print PAR "altnormstyle: NO","\n";
  print PAR "missingmode: NO","\n";
  print PAR "nsnpldregress: 0","\n";
  print PAR "noxdata: YES","\n";
  print PAR "nomalexhet: YES", "\n";

  close PAR;

}  

sub Structure {
  my ($snpdir, $res, $number, $missing, $flag, $het, $res2, $pop) = @_;
   my %genos = ('AA' => 0 , 'AC' => 1, 'CA'=>1, 'AG' => 2, 'GA' => 2,  'AT' => 3, 'TA'=>3, 'CC'=> 4, 'CG'=>5, 'GC'=>5, 'CT' => 6, 'TC' => 6, 'GG'=>7, 'GT'=>8, 'TG'=>8, 'TT'=>9, 'NN' => -9);
  # A = 1, C = 2, G = 3, T = 4
  my %genos2 = ('0' => ['1','1'], '1' => ['1','2'], '2' => ['1','3'], '3' => ['1','4'], '4'=>['2','2'], '5'=>['2','3'], '6' => ['2','4'], '7'=>['3','3'], '8'=>['3','4'], '9'=>['4','4'], '-9' => ['-9','-9'] );
  
  my $resdir = $res2 . "StructureInput/";
  mkdir $resdir unless -e $resdir;

  
  my $nondir = $1 . 'Individual_Non_diallelic/' if $res =~ /(\S+)Individual/; ##
  my %non;
  my @all = <$nondir*Contig*_Non_diallelic_SNPID.txt>;
  foreach (@all) {
    my $file = $_;
    my $lib = $1 if basename ($file) =~ /(Contig\d+)_Non_diallelic_SNPID.txt/;
    open (NON, "<", $file);
    while (<NON>) {
      chomp (my @line = split /\s+/,$_);
      $non{$lib}{$line[0]}++;
    }
    close NON;  
  }

  
  my @geno = <$res*_geno>;
  my $sample = $res."sampleID.txt";
  my $out = $resdir ."Final_Structure";
  $out = $out . "_random_informative" if $flag eq 'ri';
  $out = $out . "_random" if $flag eq 'r';
  $out = $out . "_informative" if $flag eq 'i';
  $out = $out . "_all" if $flag eq 'a';

  
  my @snp = <$snpdir*_SNP>; 
  my $file = $resdir . "Final_geno4.txt";
  open (OUT, ">", $file);
  
  foreach (@geno) {
    my $file = $_;
    my $locus = $1 if basename($file) =~ /(Contig\S+)_geno/;  
    open (IN, "<", $file);
    my $d = 1;  
    while (<IN>) {
      chomp (my $line = $_);
      print OUT $locus, "\t", $d, "\t", $line, "\n";
      $d++;
    }
    close IN;
  }
  close OUT;
 
  my %geno_hash;
  open (GENO, "<", $file);
  while (<GENO>) {
    chomp (my @a = split /\s+/, $_);
    $geno_hash{$a[0]}{$a[1]} = join ("\t",@a[2..$#a]); 	
  }
  close GENO;
  

  my $snptmp = 'tmp';
  open (OUT, ">", $snptmp );
  foreach (@snp) {
    my $file = $_;
    my $locus = $1 if basename($file) =~ /(Contig\S+)_SNP/;  
    open (IN, "<", $file);
    my $d = 1;  
    while (<IN>) {
      chomp (my $line = $_);
      print OUT $locus, "\t", $d, "\t", $line, "\n";
      $d++;
    }
    close IN;
  }
  close OUT;
  
 
  open (IN, "<", $snptmp);
  my %hash;
  my @ind;  
  while (<IN>) {
    chomp (my @line = split /\s+/,$_);
    my $geno =join('',@line[2 .. $#line]);
    $geno =~ s/-1/?/g;
    my @miss = ($geno =~ m/\?/g);
    my @het = ($geno =~ m/1/g);
    my @homR =  ($geno =~ m/0/g);
    my @homA =  ($geno =~ m/2/g);
    if ($flag eq 'ri' || $flag eq "i") {  
      if (scalar(@homA) + scalar(@het) > 1 && scalar(@homR) > 0 ) {
	if (scalar(@het)/$number <= $het && scalar(@miss)/$number <= $missing && !$non{$line[0]}{$line[1]}) {
	    $hash{$line[0]}{$line[1]} = $geno_hash{$line[0]}{$line[1]};  
	}
      }
    }
    
    if ($flag eq 'r' || $flag eq "a") {
      if ( scalar(@homR) > 0 ) {
	if (scalar(@homA) > 0 || scalar(@het) > 0) {
	  if (scalar(@het)/$number <= $het && scalar(@miss)/$number <= $missing && !$non{$line[0]}{$line[1]})  {
	    $hash{$line[0]}{$line[1]} = $geno_hash{$line[0]}{$line[1]};  
	  }
	}
      } 
    }
  }
  close IN;

  
  if ($flag eq 'r' || $flag eq 'ri') {
    my $selected = $resdir . "selected_one_SNP_per_Contig_for_structure.txt";
    my $selected1 = random (\%hash, $resdir,$selected);
    
    open (IN2, "<", $selected1);
    my $n = 0;
    while (<IN2>) {
      $n++;
    }
    close IN2;
    
    print "The number of sites selected for structure analyses is: ", $n, "\n";
    
    print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    
    my $out1 = 'cleaned';
    trim($number,$selected1, $out1,\%genos);
    my $out2= 'converted';
    transpose($out1, $number, $n, $out2);
    convert($out2, $out, \%genos2, $pop);
    unlink ($out1, $out2);
  }
  if ($flag eq 'a' || $flag eq 'i') {
    my $selected1 = $resdir . "selected_ALL_SNP_chosen_for_structure.txt";
    open (O, ">", $selected1);
    my $head = $resdir . "head";
    open (OUT, ">", $head);
    print OUT "IndID\tPopID\t";
    
    
    my $n = 0;
    
    for my $locus (sort {$a cmp $b} keys %hash){
      my $dd = -1;     
      for my $pos (sort {$a <=> $b} keys %{$hash{$locus}}){
	print OUT $dd, "\t";
	$dd += 1 if $dd > -1 ;
	$dd += 3 if $dd == -1 ;

	print O $locus, "\t", $pos, "\t", $hash{$locus}{$pos},"\n";

	$n++;
	
      }
    }
    print OUT "\n";
    close O;
    close OUT;
    
    print "The number of sites selected for structure analyses is: ", $n, "\n";
    
    print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    
    my $out1 = 'cleaned';
    trim($number,$selected1, $out1,\%genos);
    my $out2= 'converted';
    transpose($out1, $number, $n, $out2);
    convert($out2, $out, \%genos2, $pop);
    my $tmp = $resdir . "tmp";
    system ("cat $head $out > $tmp");
    system ("mv $tmp $out");
    
    unlink ($out1, $out2,$head);
  }  
}

sub convert {
  my ($file, $final, $genos2,$pop) = @_;
  my %genos2 = %{$genos2};
  open (IN, '<', $file);
  open (OUT,'>', $final);
  
  open (POP, '<',$pop);
  my $d = 1;
  my %pop;
  foreach (<POP>) {
    chomp(my @line = split /\s+/, $_);
    $pop{$d} = $line[1];
    $pop{$d} =~ s/\r//g;
    $d++;
  }
  my $i = 1;
  foreach (<IN>) {
    my @array1;
    my @array2;
    
    chomp (my @line = split /\s+/, $_);
    
    foreach my $gene (@line) {
      push (@array1, "\t", $genos2{$gene}[0]);
    }
    print OUT $i, "\t", $pop{$i},"\t", @array1, "\n";
    foreach my $gene2 (@line) {
      push (@array2, "\t",$genos2{$gene2}[1]);
    } 
    print OUT $i,"\t", $pop{$i},"\t", @array2, "\n";
    $i++;
  }
  #unlink ($file); 
}


sub trim{
  my ($sample, $file, $out, $genos) = @_;
  my %genos = %{$genos};
  open (GENOPROB, '<', $file); 
  open (OUT, '>', $out);
  foreach (<GENOPROB>) {
    my @geno = split /\s+/, $_; 
    for (my $i = 2; $i<=$sample ;$i++) {
      print OUT $genos{$geno[$i]}, "\t";
    }
    print OUT $genos{$geno[$sample+1]},"\n";
  }
  close GENOPROB;
  close OUT;
  #unlink ($file);
}

sub random {
  my ($hash, $dir, $selected) =@_;
  my %hash =%{$hash};
  open (OUT, ">", $selected);
  my @final;
  foreach my $id (sort {$a cmp $b} keys %hash) {    
    my @array = (); 
    foreach my $pos (sort {$a <=> $b} keys %{$hash{$id}}){
      push (@array, $pos);
    } 
    my $randomelement = $array[rand @array];    
    print OUT $id, "\t", $randomelement, "\t", join("\t",$hash{$id}{$randomelement}), "\n" if ($randomelement);  
  } 
  close OUT;
  return $selected;
}


sub dir {
  my ($b) = @_;
  my @a = @{$b};
  my $dir;
  if ($a[0] =~ m/\/$/ ){
    $dir= $a[0]; 
  }
  else {
    $dir = $a[0] . "/";
  }
  return $dir;
}


sub SNAPP {
  my ($res, $number, $missing, $flag, $het, $res2) = @_;
  
  my $resdir = $res2 . "SNAPPInput/";
  mkdir $resdir unless -e $resdir;

   my $nondir = $1 . 'Individual_Non_diallelic/' if $res =~ /(\S+)Individual/; ##
  my %non;
  my @all = <$nondir*Contig*_Non_diallelic_SNPID.txt>;
  foreach (@all) {
    my $file = $_;
    my $lib = $1 if basename ($file) =~ /(Contig\d+)_Non_diallelic_SNPID.txt/;
    open (NON, "<", $file);
    while (<NON>) {
      chomp (my @line = split /\s+/,$_);
      $non{$lib}{$line[0]}++;
    }
    close NON;  
  }

  
  my @SNP = <$res*_SNP>;
  my $sample = $res."sampleID.txt";
  my $out = "Final_SNAPP";
  $out = $out . "_random_informative" if $flag eq 'ri';
  $out = $out . "_random" if $flag eq 'r';
  $out = $out . "_informative" if $flag eq 'i';
  $out = $out . "_all" if $flag eq 'a';

  
  my $file = $resdir . "Final_geno2.txt";
  open (OUT, ">", $file);
  
  foreach (@SNP) {
    my $file = $_;
    my $locus = $1 if basename($file) =~ /(Contig\S+)_SNP/;  
    open (IN, "<", $file);
    my $d = 1;  
    while (<IN>) {
      chomp (my $line = $_);
      print OUT $locus, "\t", $d, "\t", $line, "\n";
      $d++;
    }
    close IN;
  }
  close OUT;
 
 
  open (IN, "<", $file);
  my %hash;
  my @ind;  
  while (<IN>) {
    chomp (my @line = split /\s+/,$_);
    my $geno =join('',@line[2 .. $#line]);
    $geno =~ s/-1/?/g;
    my @miss = ($geno =~ m/\?/g);
    my @het = ($geno =~ m/1/g);
    my @homR =  ($geno =~ m/0/g);
    my @homA =  ($geno =~ m/2/g);
    if ($flag eq 'ri' || $flag eq "i") {  
      if (scalar(@homA) + scalar(@het) > 1 && scalar(@homR) > 0 ) {
	if (scalar(@het)/$number <= $het && scalar(@miss)/$number <= $missing && !$non{$line[0]}{$line[1]}) {
	    $hash{$line[0]}{$line[1]} = $geno;	  
	}
      }
    }
    
    if ($flag eq 'r' || $flag eq "a") {
      if ( scalar(@homR) > 0 ) {
	if (scalar(@homA) > 0 || scalar(@het) > 0) {
	  if (scalar(@het)/$number <= $het && scalar(@miss)/$number <= $missing && !$non{$line[0]}{$line[1]} ) {
	    $hash{$line[0]}{$line[1]} = $geno;
	  }
	}
      } 
    }
  }
  close IN;
 
  open (ID, "<", $sample);  
  while (<ID>) {
    chomp (my $line = $_);
    push (@ind, $line);
  }
  close ID;
  
  if ($flag eq 'r' || $flag eq 'ri') {
    my $selected = $resdir . "selected_one_SNP_per_Contig_for_SNAPP.txt";
    my $selected1 = random (\%hash, $resdir, $selected);
    my $allelefreq = $resdir . "selected_one_SNP_per_Contig_for_SNAPP.maf";
    allelfreq ($selected1, $allelefreq);
    nexus (\@ind, $selected1, $out, $number,$resdir);
  }
  if ($flag eq 'a' || $flag eq 'i') {
    my $selected1 = $resdir . "selected_ALL_SNP_chosen_for_SNAPP.txt";
    open (O, ">", $selected1);
    for my $locus (sort {$a cmp $b} keys %hash){
      for my $pos (sort {$a <=> $b} keys %{$hash{$locus}}){
	print O $locus, "\t", $pos, "\t", $hash{$locus}{$pos},"\n";
      }
    }
    close O;
    my $allelefreq = $resdir . "selected_ALL_SNP_chosen_for_SNAPP.maf";
    allelfreq ($selected1, $allelefreq);
    nexus (\@ind, $selected1, $out, $number, $resdir);
  }
  
}  

sub allelfreq {
  my ($allele, $selected) = @_;
  open (O, "<", $allele);
  open (OUT, ">", $selected);
  while (<O>) {
    chomp (my @line = split /\s+/, $_);
    chomp (my @geno = split //,$line[2]);
    print OUT $line[0], "\t", $line[1], "\t";      
    my $total;
    my $minor;
    foreach (@geno) {
      
      my $m = $_;
      if ($m eq "?") {
	print OUT $m,"\t";
      }
      else {
	$total += 2 if $m == 2;
	$total += 2 if $m == 0 ;
	$total += 2 if $m == 1 ;
	
	$minor += 2 if $m == 2;
	$minor += 1 if $m == 1;
	print OUT $m,"\t";
      }
    }
    my $freq =  $minor/$total;
    $freq = 1-$freq if $freq >= 0.5;
    print OUT sprintf("%.2f", $freq), "\n";
    
  }
  close O;
  close OUT;
}


 sub nexus {
    my ($indinfo, $selected, $final, $sample, $resdir) = @_;
    my @ind = @{$indinfo};
    
    open (IN, "<", $selected);
    open (OUT, ">", "trimmed");
    my $n = 0;
    while (<IN>) {
      
      chomp(my @line = split /\s+/, $_);
      print OUT (join "\t", @line[2..$#line]), "\n";
      $n++;   
    }
    close IN; 
    close OUT;  
    print "The number of selected sites for SNAPP is: ", $n, "\n";
    print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    
    my $out = 'transposed';
    open (OUT1, '>', $out);
    for (my $i = 1; $i <= $sample; $i++) {
      my $count = $n;
      open (IN1, '<', 'trimmed');
      foreach (<IN1>) {
	chomp();
	if ($count > 1) {
	  my @line = split //, $_;
	  print OUT1 $line[$i-1], "\t";
	  $count--;
	}
	elsif ($count == 1) {
	  my @line = split //, $_;
	  print OUT1 $line[$i-1], "\n";
	}
	else {
	  last;
	}
      }
    }
    close IN1;
    close OUT1;
    
    open (IN, "<", 'transposed');
    my $results =  $resdir . $final . ".nexus";
    open (OUT, ">", $results);
    
    print OUT "#NEXUS", "\n\n";
    print OUT "Begin data;", "\n";
    print OUT "\tDimensions ntax=$sample nchar=$n;","\n";
    print OUT "\tFormat symbols=","\"012\""," missing=?;","\n";
    print OUT "\tMatrix","\n\n"; 
    my $a = 0;
    
    while (<IN>) {
      chomp (my @line = $_);
      my @d = split /\s+/, $line[0];
      print OUT $ind[$a], "\t", @d, "\n";
      $a++;
    }
    print OUT "\t", ";","\n";
    print OUT "End;", "\n";
    system("rm trimmed transposed");
  }

