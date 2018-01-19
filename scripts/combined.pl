#!/usr/bin/env perl
use warnings;
use strict;
use List::MoreUtils qw/ uniq /;

die(qq/

Usage: seqCapture combine  [options]

options:

-c     DIR      folder with all alignments of target regions
-f     DIR      folder with all alignments of flanking regions
-r     FIR      result folder

\n\n/) unless (@ARGV);

my %opts = (c=>undef,f=>undef, r=>undef);
getopts('c:f:r:', \%opts);

my $codingfolder = redir ($opts{c}) || die "must provide the path to a folder with alignments of coding regions!";
my $flankingfolder = redir ($opts{f}) || die "must provide the path to a folder with alignments of flanking regions!";
my $resdir = redir ($opts{r});
mkdir $resdir unless -e $resdir;

my %genos = ('Y' => '1' , 'M' => '1', 'R'=> '1', 'S' => '1', 'K' => '1',  'W' => '1', 'A'=> '2', 'C'=> '2', 'T'=> '2', 'G'=> '2', 'N' => '-1', '-' => '-1');
my %genos2 = ('Y' => 'CT' , 'M' => 'AC', 'R'=>'AG', 'S' => 'GC', 'K' => 'GT',  'W' => 'AT', 'A'=> 'AA', 'C'=> 'CC', 'T'=> 'TT', 'G'=> 'GG', 'N' => 'NN', '-' => 'NN');
my %genos3 = ('Y' => ['C','T'] , 'M' => ['A','C'], 'R'=>['A','G'], 'S' => ['G','C'], 'K' => ['G','T'],  'W' => ['A','T'], 'A'=> ['A','A'], 'C'=> ['C','C'], 'T'=> ['T','T'], 'G'=> ['G','G'], 'N' => ['N','N'], '-' => ['N','N']);



my @caln = <$codingfolder*aln>;
my @faln = <$flankingfolder*aln>;


my %seq1; 

foreach (@caln) {
  my $file = $_;   
  my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
  open (IN, "<", $file);
  my $i;  
  
  while (<IN>) {
    chomp (my $line = $_);
    if ($line =~ /^>(\S+)/){
      $i = $1;         
    }
    else {
      $seq1{$locus}{$i}{'length'} +=  length ($line);
      $seq1{$locus}{$i}{'seq'} .= $line;
    }
  }
      close IN;
}


my %seq2; 

foreach (@faln) {
  my $file = $_;   
  my $locus = $1 if basename($file) =~ /(Contig\S+)\.aln/;
  open (IN, "<", $file);
  my $i;      
  while (<IN>) {
    chomp (my $line = $_);
    if ($line =~ /^>(\S+)/){
      $i = $1;         
    }
    else {
      $seq2{$locus}{$i}{'length'} +=  length ($line);
      $seq2{$locus}{$i}{'seq'} .= $line;
    }
  }
  close IN;
}


my %combined;  
foreach my $marker (keys %seq1) {
  if ($seq2{$marker}) {
    foreach my $sample (keys %{$seq1{$marker}}) {
      
      $combined{$marker}{$sample}{'length'} =  $seq1{$marker}{$sample}{'length'} +  $seq2{$marker}{$sample}{'length'};
      $combined{$marker}{$sample}{'seq'} =  $seq1{$marker}{$sample}{'seq'} . $seq2{$marker}{$sample}{'seq'};
      $combined{$marker}{$sample}{'coding'} = $seq1{$marker}{$sample}{'length'};
      
    }
    delete $seq1{$marker};
    delete $seq2{$marker};
  } 
}
foreach my $marker (keys %seq1) {
  foreach my $sample (keys %{$seq1{$marker}}) {
    $combined{$marker}{$sample}{'length'} =  $seq1{$marker}{$sample}{'length'};
    $combined{$marker}{$sample}{'seq'} =  $seq1{$marker}{$sample}{'seq'} ;
    $combined{$marker}{$sample}{'coding'} = $seq1{$marker}{$sample}{'length'};
  }
}

foreach my $marker (keys %seq2) {
  foreach my $sample (keys %{$seq2{$marker}}) {
    $combined{$marker}{$sample}{'length'} =  $seq2{$marker}{$sample}{'length'};
    $combined{$marker}{$sample}{'seq'} =  $seq2{$marker}{$sample}{'seq'} ;
    $combined{$marker}{$sample}{'coding'} = 0;
  }
}





foreach my $marker (sort {$a cmp $b} keys %combined) {  
  my $site = scalar keys %{$combined{$marker}};
  
  my $sampleID = $resdir . "sampleID.txt";
  open (ID, ">", $sampleID );
  my $last = 0;
  foreach my $sample (sort {$a cmp $b}  keys %{$combined{$marker}}) {
    print ID $sample, "\n";
    $last++;
    last if ($last == $site);
  }
  close ID;
}




foreach my $locus (sort {$a cmp $b} keys %combined) {  
  my $infile = $resdir . $locus . "_tmp1";
  my $seqs;
  my $coding;
  
  my $site = scalar keys %{$combined{$locus}};
  
  my $alns = $resdir . $locus . ".aln";
  open (ALN, ">", $alns);
  
  open (OUT, ">", $infile);
  foreach my $sample (sort {$a cmp $b}  keys %{$combined{$locus}}) {
    print ALN ">", $sample, "\n";
    print ALN $combined{$locus}{$sample}{'seq'}, "\n";
    my @a = split //, $combined{$locus}{$sample}{'seq'};
    $seqs = $combined{$locus}{$sample}{'length'};
    $coding = $combined{$locus}{$sample}{'coding'};
    
    foreach (@a) {
      print OUT $_, "\t";
    }
    print OUT "\n";
  }##foreach my $sample (sort {$a cmp $b}  keys %{$combine{$marker}}) {
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
  
  my $codingend = $resdir .  $locus . "_coding_end_position.txt";
  open (CODE, ">", $codingend);
  print CODE $coding, "\n";
  
  while (<TRANS>) {
    my $minor;#############
    print SNPID $dd, "\t";#############
    $dd++;#############
    my @minor;
    chomp (my @nu = split /\s+/, $_);
    my @array;	
    my @nu2 = @nu;
    
    
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
    
    
    
  } ## while (<TRANS>)
  close SNP;
  close GENO;
  close SNPID;
  
  
  my $non_dialleic = $resdir .  $locus . "_Non_diallelic_SNPID.txt"; ##################
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
  
  
  
  
  
} ##foreach my $locus (sort {$a cmp $b} keys %combine) {



my $non_di =  $resdir . "Individual_Non_diallelic/";
mkdir $non_di unless -e $non_di;

my $fake = $resdir . "empty_Non_diallelic_SNPID.txt";

open (FAKE, ">", $fake);
print FAKE "test", "\t", "1","\n";
close FAKE;

system ("mv $resdir*Non_diallelic_SNPID.txt $non_di ");
system ("cp $resdir*sampleID.txt $non_di");

my $coding_end_id = $resdir . "Individual_coding_end_position/";
mkdir $coding_end_id  unless -e $coding_end_id;
system ("mv  $resdir*coding_end_position.txt $coding_end_id ");
system ("cp $resdir*sampleID.txt $coding_end_id");


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
