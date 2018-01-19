#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;
use Getopt::Std;
use List::Util qw[min max];
use List::Util qw(sum);
use List::MoreUtils qw(part);
use List::MoreUtils qw/ uniq /;
use Tie::Array::Packed;
use File::Temp;
use List::Util 'shuffle';

die(qq/

Usage: seqCapture align [options]

options:

-e     DIR      folder with all filtered individual fasta files from target regions
                (AAA.filtered.fasta, BBB.filtered.fasta, CCC.filtered.fasta...)
-b     FILE     A txt file that contains names of libraries that you want 
                to align together, one name per line
-a     INT      Method for alignment 1=align by DNA (non-coding); 2=align by codon [2]
-m     FLOAT    Percent missing data defined as no data [0.3]
-n     FLOAT    Percent individuals that fail to pass -m filter 
                defined as missing loci [0.3]
-l     FLOAT    l th percentile to trim off loci with extremely low depth [1]
-h     FLOAT    h th percentile to trim off loci with extremely high depth [99]
-H     FLOAT    H th percentile to trim off loci with extremely 
                high average heterozygosity [100]
-z     FLOAT    Maximum proportion of shared polymorphic sites in a locus [0.2]
-c     FLOAT    Removes all positions in the alignment with gaps in 1-c percent or 
                more of the sequences [0.5]
-t     INT      number of threads used in alignment [5]
-M              mitocondrial? [null]

\n\n/) unless (@ARGV);
  
my %opts = (e=>undef,b=>undef, t=>5, m=>0.3, n=>0.3, l=>1, h=>99, H=>100, z=>0.2, T=>1,a=>2, c=>0.5);
getopts('e:b:m:n:l:h:H:z:t:a:T:c:M', \%opts);

my $dir = redir ($opts{e});

my $refdir = $dir . "alignment/";
mkdir $refdir unless -e $refdir;
my $per = $opts{c}; 
my $thread = $opts{t}; 
my $missingData = $opts{m}; 
my $missingLoci = $opts{n}; 
my $names =  $opts{b} ;
my $trimaln = $opts{T};
my $lowper = $opts{l};  
my $highper = $opts{h}; 
my $highH = $opts{H};
my $het =  $opts{z};
my $coding = $opts{a};
my $M;
$M = 1 if $opts{M};
$M = 2 if !$opts{M};

open (IN,"<", $names);
my $sample = 0;
while (<IN>) {
  $sample ++ unless $_ =~ /^$/; 
}
close IN;

my @allh = <$dir*individual_H.txt>;
my @allcov = <$dir*loci_depth.txt>;

###filter loci based on coverage percentile
processCov (\@allcov,$highper, $lowper, $refdir);

###filter loci based on het percentile
processH (\@allh, $highH,$refdir);

sub processH {
  my ($allh, $highH, $dir) = @_;
  my @allh = @{$allh};
  
  foreach (@allh) {
    my $file = $_;
    my $lib = $1 if $file =~ /(\S+)_individual_H.txt/;
    my $name = $1 if basename ($file) =~ /(\S+)_individual_H.txt/;
    my $kept = $lib . "_individual_H_kept.txt"; 
    my $filtered = $lib . "_individual_H_filtered.txt"; 
    my $hDir;
    $hDir = $dir . "IndividualH_filtered/";
    mkdir $hDir unless -e $hDir;
    
    open (IN, "<", $file);
    my @h;
    while (<IN>) {
      chomp (my @d = split /\s+/, $_);
      push @h, $d[1];
    }
    seek IN, 0,0;
    
    my $h_percentile = $lib . "_individual_H_percentile.txt";  
    my ($low, $high) = percentile (\@h, $h_percentile, '0.5', $highH);
    $high = $high + 1000 if $highH == 100;  
    
    open (KEPT, ">", $kept);
    open (FILTER, ">", $filtered);
    my $d = 0;
    while (<IN>) {
      chomp (my @d = split /\s+/, $_);
      if ($d[1] > $high) {
	print FILTER $d[0], "\t", $d[1], "\n"; 
	$d++;
      }
      else {
	print KEPT $d[0], "\t", $d[1], "\n";
      }   
    }
    close IN;
    close KEPT;
    close FILTER;
    system ("mv $filtered $hDir");
    system ("rm $kept $h_percentile"); 
    print "In ",$name, ", ", $d, " loci were removed because they failed to pass individual heterozygosity filter!", "\n" ; 
  } ##foreach (@allh) {
}

sub  processCov {
  my ($allcov,$highper, $lowper, $dir) = @_;
  my @allcov = @{$allcov};
  
  foreach (@allcov) {
    my $file = $_;
    my $lib = $1 if $file =~ /(\S+)_loci_depth.txt/;
    my $name = $1 if  basename ($file) =~ /(\S+)_loci_depth.txt/;
    
    my $kept = $lib . "_loci_depth_kept.txt"; 
    my $filtered = $lib . "_loci_depth_filtered.txt"; 
    my $covDir;
    $covDir = $dir . "Coverage_filtered/";
    mkdir $covDir unless -e $covDir; 
    
    open (IN, "<", $file);
    my @cov;
    while (<IN>) {
      chomp (my @d = split /\s+/, $_);
      push @cov, $d[1];
    }
    seek IN, 0,0;
    
    my $gene_percentile = $lib . "_loci_depth_percentile.txt";  
    my ($low, $high) = percentile (\@cov, $gene_percentile, $lowper, $highper);
    
    $low = $low - 1000 if $lowper == 0;
    $high = $high + 1000 if $highper == 100;
    
    open (KEPT, ">", $kept);
    open (FILTER, ">", $filtered);
    my $d = 0;
    while (<IN>) {
      chomp (my @d = split /\s+/, $_);
    if ($d[1] > $high || $d[1] < $low) {
      print FILTER $d[0], "\t", $d[1], "\n"; 
      $d++;
    }
      else {
	print KEPT $d[0], "\t", $d[1], "\n";
      }   
    }
    close IN;
    close KEPT;
    close FILTER;
    system ("mv $filtered $covDir");
    
    print "In ",$name, ", ", $d, " loci were removed because they failed to pass individual coverage filter!", "\n";
    system ("rm $kept $gene_percentile");
  } ##foreach (@allh) {  
}

my $filtered = $refdir . "combined_filtered.fasta";

my $d1 = $refdir . "IndividualH_filtered/";
my $d2 = $refdir . "Coverage_filtered/";
system (" cat $d1*individual_H_filtered*  $d2*loci_depth_filtered* > $filtered ");

my %filter;
open (IN, "<", $filtered);
while (<IN>) {
  chomp (my @d = split /\s+/, $_);
  $filter{$d[0]}++;  
}
close IN;   
unlink ($filtered);

my @targetdata = <$dir*filtered.fasta>;
foreach (@targetdata) {
  my $file = $_;
  my $name = $1 if   $file =~ /(\S+)_filtered.fasta/;
  my $new = $file . ".2";
  open (IN, "<", $file);
  open (OUT, ">", $new);
  while (<IN>) {
    my $line = $_; 
    if ($line =~ m/^>(\S+)/ ) { 
      my $id = $1;
      chomp (my $seq = <IN>);	 
      unless ($filter{$id}) {
	print OUT ">", $id, "\n", $seq, "\n"; 
      }
    }
  }
  close IN;
  close OUT; 
} 

####now making MSAs
print "\nNow looking for orthologous markers...\n";
MakeAlignment ($refdir, $dir,  $names, $missingData, $missingLoci, $sample, $lowper, $highper, $highH,$het, $trimaln, $M, $coding, $thread, $per);
print "\nProgram finished successfully!\n\n";


sub MakeAlignment {
  my ($resdir, $dir, $file, $missingData, $missingLoci, $sample, $min, $max,  $H,$het, $trimaln,$M,$coding, $thread, $per) = @_;
  
  my $misDir = $resdir . "missingData_filtered/";
  mkdir $misDir unless -e $misDir;
  
  my $hetDir = $resdir . "HetSites_filtered/";
  mkdir $hetDir unless -e $hetDir;
    
  my %genos = ('Y' => '1' , 'M' => '1', 'R'=> '1', 'S' => '1', 'K' => '1',  'W' => '1', 'A'=> '2', 'C'=> '2', 'T'=> '2', 'G'=> '2', 'N' => '-1', '-' => '-1');
  my %genos2 = ('Y' => 'CT' , 'M' => 'AC', 'R'=>'AG', 'S' => 'GC', 'K' => 'GT',  'W' => 'AT', 'A'=> 'AA', 'C'=> 'CC', 'T'=> 'TT', 'G'=> 'GG', 'N' => 'NN', '-' => 'NN');
  my %genos3 = ('Y' => ['C','T'] , 'M' => ['A','C'], 'R'=>['A','G'], 'S' => ['G','C'], 'K' => ['G','T'],  'W' => ['A','T'], 'A'=> ['A','A'], 'C'=> ['C','C'], 'T'=> ['T','T'], 'G'=> ['G','G'], 'N' => ['N','N'], '-' => ['N','N']);
  
  my $start1 = time;
  my $new_master = $resdir . "combined_filtered.fasta"; 
  system (" cat $dir*filtered.fasta.2 > $new_master ");
  
  my %name;
  my @nameorder;
  my $d=1;
  open (NAME, "<", $file);
  while (<NAME>) {
    chomp (my $line = $_);
    unless (/^$/) {
      push @nameorder, $line;
      $name{$d} = $line;
      $d++;
    }
  }
  close NAME;
  $d = $d-1;
  
  my $hash = seqhash2($new_master);
  my %seq = %{$hash};
  
  unlink ($new_master);  
  
  my %alldata; 
  my $missing;
  
  foreach my $number (sort {$a cmp $b} keys %seq) {
    foreach my $lib (sort {$a cmp $b} keys %{$seq{$number}}) {   
      push @{$alldata{$number}}, $lib. "_" . $number;
    }
  }
  
  foreach my $number (sort {$a cmp $b} keys %seq) {
    my $muscle_in = $resdir  . $number;
    open (MUSIN, ">", $muscle_in );
    #    my %corr;
    foreach my $d (sort {$a <=> $b} keys %name) {
      my $yes = 0; 
      foreach my $lib (sort {$a cmp $b} keys %{$seq{$number}}) {
	if ($name{$d} eq $lib) {
	  print MUSIN ">", $name{$d}, "\n";
	  my $seq = $seq{$number}{$lib}{'seq'};
	  $seq =~ s/[Y|S|K|R|W|M|y|s|k|r|w|m]/N/g if ($M == 1);
	  $seq =~ s/[N|n]//g  if ($coding == 1);  
	  print MUSIN $seq, "\n";
	  $yes = 1;
	  last;
	}
      }
      
      if ($yes == 0) {
	print MUSIN ">", $name{$d}, "\n"  if $coding == 1;
	print MUSIN '---------------------------------', "\n" if $coding == 1;
	#       print MUSIN 'AAA', "\n" if $coding == 2;
	#	$corr{$name{$d}} = 'AAA' if $coding == 2;
	#	$corr{$name{$d}} = '---------------------------------' if $coding == 1;
      }  
    } #foreach my $d (sort {$a <=> $b} keys %name) { 
    close MUSIN;
    
    my $muscle_out = $resdir . $number . ".aln";
    
    if ($coding == 1) {
      system ("mafft --inputorder  --leavegappyregion --threadit 0  --allowshift --unalignlevel 1 --anysymbol --retree 1  --thread $thread --maxiterate 1000  --globalpair  $muscle_in  > $muscle_out "); 
      unlink ($muscle_in);   
    }
    
    if ($coding == 2) {
      my $prankin = $resdir . $number . ".prank_in";
      #system ("mv $muscle_in $prankin");           
      open (IN, "<", $muscle_in);
      open (OUT, ">", $prankin);
      while (<IN>) {
	chomp (my $line = $_);
	if ($line =~ /^>(\S+)/) {
	  my $id = $1;
	  chomp (my $seq = <IN>);
	  my $sequence;	  
	  my @a = split //, $seq;
	  my $count = int ((scalar @a) /3);
	  @a = splice (@a, 0, 3*$count);
	  while(my ($a,$b,$c) = splice(@a,0,3)) {	    
	    my $string = $a . $b . $c;
	    if ($string  =~ /N|0/) {
	      #print "yes", "\n";
	      #print $string , "\n";
	      $sequence .= "---";	      
	    }
	    elsif ($string !~ /N|0/) {
	      #print $string , "\n";
	      $sequence .= $string ;	      
	    }
	    else {
	      print $string , "\n";
	      print "really? there are other conditions?\n";
	      exit;
	    } 
	  }
	  print OUT ">", $id, "\n";
	  print OUT $sequence,"\n";
	}
      }
      close IN;
      close OUT;
      unlink ($muscle_in);
      system ("prank -d=$prankin -o=$prankin -f=" . "'". 'fasta' . "'" . " -codon -iterate=5  ");
      my $prankout = $prankin . ".best.fas";
      open (PRANKOUT, "<", $prankout);
      my %prankout;
      my $ids;
      while (<PRANKOUT>) {
	chomp (my $line = $_);
	if ($line =~ m /^>(\S+)/) {
	  $ids = $1;
	}         
	else {
	  $prankout{$ids} .= $line;
	}	
      }
      close PRANKOUT;
      unlink ($prankout);
      
      my $alnlen = 0;
      foreach my $lib ( keys %prankout) {
	$alnlen = length $prankout{$lib};
	last;
      }
      
      my $prankout2 = $prankin . ".best.fas2";
      
      open (OUT5, ">", $prankout2); 
      foreach my $name (@nameorder) {
	if ($prankout{$name}) {
	  print OUT5 ">", $name, "\n";
	  print OUT5 $prankout{$name}, "\n";
	}
	else {
	  print OUT5 ">", $name, "\n";
	  my $dash = '-' x $alnlen;
	  print OUT5 $dash, "\n";
	}	
      }  
      close OUT5;      
      unlink ($prankin);
      system ("mv $prankout2 $muscle_out");
    } ##$coding == 2;
    
    my $muscle_out1 = $resdir . $number . ".aln1";
    
    if ($trimaln == 1) {    ###use trimal to trim gaps and ambigugous regions from the alignment, parameters are not optimized.  
      #system ("trimal -in $muscle_out -out $muscle_out" . ".fa"."  -block 20  -resoverlap 0.3 -seqoverlap 30  -keepseqs -gappyout ") if $coding == 1;
      #system ("trimal -in $muscle_out -out $muscle_out" . ".fa". " -block 15  -resoverlap 0.3 -seqoverlap 30  -keepseqs -gappyout ") if $coding == 2;
      system ("trimal -in $muscle_out -out $muscle_out" . ".fa"."   -keepseqs -gt $per ") if $coding == 1;
      system ("trimal -in $muscle_out -out $muscle_out" . ".fa". "  -keepseqs -gt $per ") if $coding == 2;
      my $trim = $muscle_out . ".fa";
      
      my $counts = 0;
      open (GB, "<", $trim);
      while (<GB>) {
	chomp (my $line = $_);
	$counts ++ if $line !~ m /^$/;  
      }
      seek GB, 0,0;
      
      if ( $counts > 0) {
	open (GB, "<", $trim);
	my $ll;
	my $tid;
	my %trim1;
	while (<GB>) {
	  chomp (my $line = $_);
	  if ($line =~ m/^>(\S+)/) {
	    $tid = $1;
	  }
	  else {       
	    $trim1{$tid} .= $line;
	  }
	}
	close GB;
	$ll = length ($trim1{$tid});
	
	my $trim2 = $trim. ".1";
	open (GB2, ">", $trim2);
	
	my %gene;
	foreach my $d (sort {$a <=> $b} keys %name) { 
	  my $yes = 0; 
	  
	  foreach my $number1 (sort {$a cmp $b} keys %trim1) {    
	    if ($name{$d} eq $number1) {	
	      
	      my $seq = $trim1{$number1};
	      
	      my $a = ($seq =~ s/[A]//ig);
	      my $tcg = ($seq =~ s/[TCG]//ig);
	      
	      print GB2 ">", $number1, "\n";
	      
	      if ($a == 3 && $tcg == 0) {
		my $empty = '-' x $ll;
		print GB2 $empty, "\n";
		chomp (my @a = split /\s+/,$empty);
		$gene{$number1} = uc (join "", @a);
	      }
	      
	      else {
		chomp (my @a = split /\s+/,$trim1{$number1});
		$gene{$number1} = uc (join "", @a);
		#following codes are stupid but I don't have time to figure a smarter way to write them!
		my $seq1 = $trim1{$number1};
		$seq1 =~ s/-[A-Z]-/---/g;	      
		$seq1 =~ s/-[A-Z][A-Z]-/----/g;
		$seq1 =~ s/-[A-Z][A-Z][A-Z]-/-----/g;
		$seq1 =~ s/-[A-Z][A-Z][A-Z][A-Z]-/------/g;
		$seq1 =~ s/-[A-Z][A-Z][A-Z][A-Z][A-Z]-/-------/g;
		$seq1 =~ s/-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-/--------/g;
		$seq1 =~ s/-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-/---------/g;
		$seq1 =~ s/-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]-/----------/g;

		
		$seq1 =~ s/-[A-Z]$/--/;
		$seq1 =~ s/-[A-Z][A-Z]$/---/;
		$seq1 =~ s/-[A-Z][A-Z][A-Z]$/----/;
		$seq1 =~ s/-[A-Z][A-Z][A-Z][A-Z]$/-----/;
		$seq1 =~ s/-[A-Z][A-Z][A-Z][A-Z][A-Z]$/------/;
		$seq1 =~ s/-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]$/-------/;
		$seq1 =~ s/-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]$/--------/;
		$seq1 =~ s/-[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]$/---------/;
		
		print GB2 $seq1, "\n" 
	      }
	      
	      $yes = 1;
	      last;
	    }
	  }
	  
	  if ($yes == 0) {
	    my $empty = '-' x $ll;
	    chomp (my @a = split /\s+/,$empty);
	    $gene{$name{$d}} = uc (join "", @a);
	    print GB2 ">", $name{$d}, "\n";
	    print GB2  $empty, "\n";	  
	  } 
	  
	} #foreach my $d (sort {$a <=> $b} keys %name) { 
	close GB2;
	
	
	unlink ($trim) if $trim;
        #my $mus = $muscle_out . ".copy";
	#system ("cp $muscle_out $mus");
	system ("mv $trim2 $muscle_out") ;
	my $N;
	
	my @genearray;
	foreach my $c (sort {$a cmp $b} keys %gene) {
	  my $seq = $gene{$c};	
	  my $atgc = ($seq =~ s/[ATGCYWRMKS]//ig);
	  push @genearray, $atgc;
	} 
	my $maxLen = max (@genearray);
	my $md = 0;
	foreach my $c (sort {$a cmp $b} keys %gene) {
	  my $seq = $gene{$c};	
	  my $l = ($seq =~ s/[ATGC]//ig);	
	  my $missing = 1 - $l/$maxLen;
	  $md ++ if ($missing > $missingData);
	} 
	system ("mv $muscle_out $misDir ") if ($md/$sample > $missingLoci || $maxLen <= 10 );
	$missing++ if ($md/$sample > $missingLoci || $maxLen <= 10 );
	my $delete = $1 if $muscle_out =~ /(Contig\S+)\.aln/;
	delete $alldata{$delete} if ($md/$sample > $missingLoci);
	delete $alldata{$delete} if ($maxLen <= 100);      	
      }
      
      if ($counts == 0) {
	$missing++;
	system ("mv $muscle_out $misDir ");
	unlink ($trim) if $trim;
      }
      
    } ## if ($trimaln = 1) {
    
    
    if ($trimaln == 2) {
      my $gb = $muscle_out . ".fa";
      system ("Gblocks $muscle_out -t=c -e=.fa -p=n -b5=a -b4=10 -b3=4");
      
      open (GB, "<", $gb);
      my $yes = 0; 
      
      while (<GB>) {
        chomp (my $line = $_);
        $yes ++ if ($line !~ m /^>/ && $line !~ m /^$/);
      }
      close GB;

      if ($yes == 0) { 
	$missing++;
	my $delete = $1 if $muscle_out =~ /(Contig\S+)\.aln/;
	delete $alldata{$delete} if $alldata{$delete};
	system ("mv $muscle_out $misDir ") ;
	unlink ($muscle_out1,$muscle_out,$gb);
      }

      if  ($yes > 0) { 
	system ("mv $gb $muscle_out");
	
	
	open (IN, "<", $muscle_out);
	open (OUT, ">", $muscle_out1);
	my %gene;      
	my $id;
	my $N;
	while (<IN>) {
	  chomp (my $line = $_);
	  if ($line =~ m/^>(\S+)/) {
	    $id = $1;
	    print OUT ">", $id,"\n";
	  }
	  else {
	    chomp (my @a = split /\s+/,$line);
	    $gene{$id} .= uc (join "", @a);
	    print OUT uc (join "", @a), "\n";
	  }
	  
	}
	close IN;
	close OUT;
	
	system ("mv  $muscle_out1  $muscle_out ");
	
	my @genearray;
	foreach my $c (sort {$a cmp $b} keys %gene) {
	  my $seq = $gene{$c};	
	  my $atgc = ($seq =~ s/[ATGCYWRMKS]//ig);
	  push @genearray, $atgc;
	} 
	my $maxLen = max (@genearray);
	my $md = 0;
	foreach my $c (sort {$a cmp $b} keys %gene) {
	  my $seq = $gene{$c};	
	  my $l = ($seq =~ s/[ATGCYWRMKS]//ig);	
	  my $missing = 1 - $l/$maxLen;
	  $md ++ if ($missing > $missingData);
	} 
	system ("mv $muscle_out $misDir ") if ($md/$sample > $missingLoci);
	$missing++ if ($md/$sample > $missingLoci);
	my $delete = $1 if $muscle_out =~ /(Contig\S+)\.aln/;
	delete $alldata{$delete} if ($md/$sample > $missingLoci);
	delete $alldata{$delete} if ($maxLen <= 10);      
      }
      
    } ##($trimaln = 2) {
    
  } ## foreach my $number (sort {$a cmp $b} keys %seq) {
  
  
  
  my $time1 = int((time - $start1)/60);  
  print "\n Alignment is finished and it took $time1 minutes! \n";
  print "But... we need to filter out alignment with too much missing data and potential paralogs!!\n";
  print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  
  
  my $percentmissingLoci = (1-$missingLoci)*100;
  my $percentmissingData =  (1-$missingData)*100;
  my $missingdata= $missingData*100;

  
  print "\n", $missing, " loci are defined as loci with too much missing data by the user!\n" if $missing;
  print "\n", "In each alignement at least " , $percentmissingLoci, "% of the samples contain no more than ",$missingdata, "% missing data!\n";
  
  
  my %final;
  
  #my @cov;
  #my @h;
  
  
  foreach my $number (sort {$a cmp $b} keys %alldata) { 
    $final{$number}++
  }
  
  my @aln = <$resdir*aln>;
  my %fhet;
  
  foreach (@aln) {
    my $file = $_;   
    my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/; 
    open (IN, "<", $file);
    my $i;
    my %seq;
    while (<IN>) {
      chomp (my $line = $_);
      if ($line =~ /^>(\S+)/){
	$i = $1;
      }
      else {
	$seq{$i}{'length'} +=  length ($line);
	$seq{$i}{'seq'} .= $line;
      }
    }
    close IN;
    
    my $infile = $resdir .$locus . "_tmp1";
    open (OUT, ">", $infile);
    my $site = scalar keys %seq;
    my $seqs;
    
    my %seq2 = %seq;
    my $sampleID = $resdir . "sampleID.txt";
    open (ID, ">", $sampleID );
    my $last = 0;
    foreach my $sample (sort {$a cmp $b} keys %seq2) {
      print ID $sample, "\n";
      $last++;
      last if ($last == $site);
    }
    close ID;
    
    foreach my $sample (sort {$a cmp $b} keys %seq) { 
      $seqs = $seq{$sample}{'length'};
      my @a = split //, $seq{$sample}{'seq'};
      foreach (@a) {
	print OUT $_, "\t";
      }
      print OUT "\n";
    }      
    close OUT;
    
    my $transposed = $resdir .  $locus . "_tmp"; 
    transpose ($infile,  $seqs, $site, $transposed);
    
    open (TRANS, "<", $transposed);
    unlink ($infile, $transposed);
    
    my $indSNPs =  $resdir .  $locus . "_SNP";
    my $indGeno =  $resdir .  $locus . "_geno";
    open (SNP, ">", $indSNPs);
    open (GENO, ">", $indGeno);
    
    my $SNP_ID  =  $resdir .  $locus . "_SNPID.txt"; ##############
    open (SNPID, ">", $SNP_ID);##############
    my $dd =1;##############
    print SNPID "Position\tMajor\tMinor\n";############
    
    
    #my %genos = ('Y' => '1' , 'M' => '1', 'R'=> '1', 'S' => '1', 'K' => '1',  'W' => '1', 'A'=> '2', 'C'=> '2', 'T'=> '2', 'G'=> #'2', 'N' => '-1', '-' => '-1');
    #  my %genos2 = ('Y' => 'CT' , 'M' => 'AC', 'R'=>'AG', 'S' => 'GC', 'K' => 'GT',  'W' => 'AT', 'A'=> 'AA', 'C'=> 'CC', 'T'=> #'TT', 'G'=> 'GG', 'N' => 'NN', '-' => 'NN');
    #  my %genos3 = ('Y' => ['C','T'] , 'M' => ['A','C'], 'R'=>['A','G'], 'S' => ['G','C'], 'K' => ['G','T'],  'W' => ['A','T'], #'A'=> ['A','A'], 'C'=> ['C','C'], 'T'=> ['T','T'], 'G'=> ['G','G'], 'N' => ['N','N'], '-' => ['N','N']);
    
    
    while (<TRANS>) {
      my $minor;#############
      print SNPID $dd, "\t";#############
      $dd++;#############
      my @minor;
      chomp (my @nu = split /\s+/, $_);
      my @array;	
      my @nu2 = @nu;
      my $min_het = 0;
      
      foreach (@nu2) {
	push @array, $_ if $_ =~ /[A|T|C|G]/;
      }
      
      my ($item, $count) = most_Frequent(\@array) if (@array);
      
      if (@array) {
	foreach my $site (@nu) {	    
	  print GENO $genos2{$site}, "\t";  	    
	  if ($site eq $item ) {
	    if ($item =~ /[A|T|G|C]/) {   
	      print SNP "0", "\t";
	    }
	    else {
	      print "error! major should not be 'N' or '-' !", "\n";
	      exit;
	    }  
	  }
	  else {
	    $minor = $genos2{$site} if ($site ne 'N' && $site ne '-'); #############PHYLO
	    push @minor, $genos3{$site}[0], $genos3{$site}[1] if ($site ne 'N' && $site ne '-'); ######PHYLO
	    print SNP $genos{$site}, "\t";
	    $min_het ++ if ($genos{$site}) eq '1';
	  }
	} ##foreach my $site (@nu) {
	
	print SNP "\n";
	print GENO "\n";
	print SNPID $item, "\t";#############
	
	if ($minor) {###########
	  my @unique = uniq @minor; ################
	  foreach (@unique) {################
	    print SNPID $_,"\t" if $_ ne $item;#############
	  }###############
	  print SNPID "\n";###########	    
	}     #############    
	if (!$minor) {#############
	  print SNPID $item, "\n";#############
	} #############	  
      } ##if (@array)
      
      else {
	foreach my $site (@nu) {
	  $minor = $genos2{$site} if ($site ne 'N' && $site ne '-');  #############
	  push @minor, $genos3{$site}[0], $genos3{$site}[1] if ($site ne 'N' && $site ne '-'); ######PHYLO
	  print GENO $genos2{$site}, "\t";  
	  print SNP $genos{$site}, "\t";
	  $min_het ++ if ($genos{$site}) eq '1';
	}
	print SNP "\n";
	print GENO "\n";
	
	if ($minor) {#############
	  my @unique = uniq @minor; ################
	  foreach (@unique) {################
	    print SNPID $_,"\t"; #############
	  }###############
	  print SNPID "\n";###########
	} #############
	if (!$minor) {#############
	  print SNPID "N\tN\n" if (!@array);#############
	}#############
      } ## else
      
      
      
      if ($min_het/$site > $het) {
	delete ($final{$locus});
	system ("mv $file $hetDir");
	unlink ($indSNPs, $indGeno,  $SNP_ID); #####################
	$fhet{$locus}++;
	last;
      }
      
    } ## while (<TRANS>)
    close SNP;
    close GENO;      	
    close SNPID;
    if (-e $file) {
      my $non_dialleic = $resdir .  $locus . "_Non_diallelic_SNPID.txt"; ##################
      print $non_dialleic, "\n";
      open (ININ, "<", $SNP_ID);###########
      open (NONSNP, ">", $non_dialleic) ; ################
      
      while (<ININ>) {#############
	chomp (my $l = $_);#############
	my @line = split /\s+/,$l;###########
       
	print NONSNP $l, "\n" if (scalar (@line) > 3);#############
      }#############
      close ININ;###########
      close NONSNP;############
      unlink ( $non_dialleic) if (-z $non_dialleic);
    }
    
  }##foreach @aln
  
  my $perhet = $het * 100;
  print "\n", scalar keys %fhet, " loci have at least one site where ", $perhet, "% of the samples are heterozygous and are therefore filtered out !!\n";  
  #print "\n", "All the above loci failed to pass various filters are moved to subfolders and will not be used for the downstrean analyses!\n\n\n";
  
  my $non_di =  $resdir . "Individual_Non_diallelic/";
  mkdir $non_di unless -e $non_di;
  
  my $fake = $resdir . "empty_Non_diallelic_SNPID.txt";
  
  open (FAKE, ">", $fake);
  print FAKE "test", "\t", "1","\n";
  close FAKE;
  
  system ("mv $resdir*Non_diallelic_SNPID.txt $non_di ");
  system ("cp $resdir*sampleID.txt $non_di");
  
  my $snpid = $resdir . "Individual_SNPID/";
  mkdir $snpid unless -e $snpid;
  system ("mv  $resdir*SNPID.txt $snpid ");
  system ("cp $resdir*sampleID.txt $snpid");
  
  my $snps = $resdir . "Individual_SNPs/";
  mkdir $snps unless -e $snps;
  system ("mv  $resdir*_SNP $snps ");
  system ("cp $resdir*sampleID.txt $snps");
  
  my $geno =  $resdir . "Individual_GENOs/";
  mkdir $geno unless -e $geno;
  system ("mv  $resdir*_geno $geno");
  system ("cp $resdir*sampleID.txt $geno");
  
  my $alns = $resdir . "Individual_ALNs/";
  mkdir $alns unless -e $alns;
  system ("mv  $resdir*.aln $alns");
  system ("cp $resdir*sampleID.txt $alns");
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

sub percentile {
  my ($array, $out, $min, $max) = @_;
  
  open (OUT1, ">", $out);
  tie my @a, 'Tie::Array::Packed::Number';
  @a = @{$array};
  
  tied(@a)->sort;
  my $low;
  my $high;
  foreach my $id (0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,7.5,10,20,30,40,50,60,70,80,90,92.5,95,95.5, 96, 96.5,97,97.5,98,98.5,99,99.5,100) {  
    print OUT1 $id . " precent percentile: ",  sprintf("%.2f", $a[$#a*$id/100]), "\n";
    $low = sprintf("%.2f", $a[$#a*$id/100]) if $id == $min; 
    $high = sprintf("%.2f", $a[$#a*$id/100]) if $id == $max;     
  }
  close OUT1;
  return ($low, $high);
}

sub seqhash2 {
  my ($file) = @_;
  my %seq;
  open (IN, "<", $file); 
  my $id;
  my $locus;
  while (<IN>) {
    chomp (my $line = $_);
    if ($line =~ m /^>(\S+)_(Contig\d+)/) {   
      $id = $1;
      $locus = $2;
    }
    else {
      $seq{$locus}{$id}{'seq'} .= $line;
      $seq{$locus}{$id}{'len'} += length ($line);  
    }   
  }
  close IN;
  return (\%seq);
}

  
sub seqhash1 {
  my ($file) = @_;
  my %seq;
  open (IN, "<", $file); 
  my $id;
  while (<IN>) {
    chomp (my $line = $_);
    if ($line =~ m /^>(\S+)/) {   
      $id = $1;
    }
    else {
      $seq{$id}{'seq'} .= $line;
      $seq{$id}{'len'} += length ($line);  
    }   
  }
  close IN;
  return (\%seq);
}

sub most_Frequent {
  my ($array) = @_; 
  my @array = @{$array};
  
  my(%count);
  foreach my $value (@array) {
    $count{$value}++;
  }
  my $max_value = (sort {$count{$b} <=> $count{$a}} @array)[0];
  my $counts = $count{$max_value};

  if ($counts ==1 ) {
    my @array2;
    foreach (@array) {
      push @array2, $_;
      
    }
    my $ForMax=0;
    length ($array2[$ForMax]) >length($array2[$_]) or $ForMax = $_ for 1 .. $#array2;
    $max_value = $array2[$ForMax];
    
  }
  return ($max_value, $counts);
  
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
