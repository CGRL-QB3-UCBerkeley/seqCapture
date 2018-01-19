#!/usr/bin/env perl

use warnings;
use strict;
use File::Basename;
use Getopt::Std;
use List::Util qw[min max];
use List::Util qw(sum);
use List::MoreUtils qw(part);
use Tie::Array::Packed;
use List::MoreUtils qw/ uniq /;
use File::Temp;
use List::Util 'shuffle';
use Cwd 'abs_path';



die (qq/

Usage: seqCapture buildcontigs [options]

Basic Options:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-t     INT     Target sequences could be one of the following: 
               1=individual exons; 
               2=cDNA (no UTR); 
               3=transcripts (including UTR); 
               4=random (such as UCEs, no need for exon identification) 
               [3]  
-a     DIR     Path to a folder with all intarget assemblies
               (AAA_targetedRegionAndFlanking.fasta);
-f     DIR     Path to a folder with all bed files
-b     DIR     A folder with all cleaned reads (AAA_1_final.fq, 
               AAA_2_final.fq, AAA_u_final.fq...)
-i     INT     Avg. Insert size [200];
-m     INT     memory limit (in MB) for the program, default 800;
               0 for unlimitted [0]   
-n     INT     number of threads [10]     
-d     INT     Minimum depth to keep a site, otherwise masked as an "N" [5]
-D     INT     Maximum depth to keep a site, otherwise masked as an "N" [100000]
-N     INT     INDEL filtering window [5]
-M     FLOAT   Discard a locus if M percent bases are Ns [0.8]
-c     INT     only keep concordant mapping for PE reads?
               1 = yes
               0 = no [1]
-r     INT     read length (bp) of original raw reads [100]
-s     INT     repeat Masking?
               1 = yes
               0 = no [0]

RepeatMasking Options: (use T or R)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-R     CHAR    Species used for repeatmasking. some examples are: human, mouse, rattus, 
               "ciona savignyi",arabidopsis, mammal, carnivore, rodentia, rat, cow, pig,
               cat, dog, chicken, fugu, danio, "ciona intestinalis", drosophila, 
               anopheles, elegans,diatoaea, artiodactyl, rice, wheat, maize, 
               "vertebrata metazoa" 
-T     CHAR    Use a custom-build repetitive library (full path) for repeat masking, 
               in this case do not use -R  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      
\n\n/) unless (@ARGV);
  
  
my %opts = (a=>undef, b=>undef, m=>0, d=>5, D=>100000, N=>5,r=>100, M=>0.8, i=>200, R=>undef, T=>undef, t=>3,f=>undef, c=>1,s=>0,n=>10);
getopts('a:m:d:D:N:s:M:R:T:b:i:t:f:r:c:n:', \%opts);

my $Rlength =  $opts{r} + $opts{r};    
my $refdir = redir ($opts{a}) ; 
my $beddir = redir ($opts{f});
my $readDir = redir($opts{b});
my $repeat = $opts{s};
my $insert =  $opts{i};
my $mem = $opts{m};
my $vcf2fqMinDepth = $opts{d};
my $vcf2fqMaxDepth = $opts{D};
my $vcf2fqIndelFilter =  $opts{N};
my $thread = $opts{n};
my $ma = $opts{M};
my $coding = $opts{t};
my $con = $opts{c};
my $cpu = $opts{n};


#sort the bed first!!!
my @bedall = <$beddir*.bed>;
foreach (@bedall) {
  my $file = $_;
  my $file1 = $file.".sorted";
  system ("sort -k 1,1 -k 2n,2 $file > $file1");
  system ("mv $file1  $file");
  
}

####identify repeatMasking library
my $rep;
my $c;

if ($opts{R}) {
  $rep = $opts{R};
  $c =1;
  my @rep = split /_/, $rep;
  if (scalar (@rep) > 1) {
    $rep = '"' . $rep[0] ." ". $rep[1] . '"';
  }
}
if ($opts{T}) {
  $rep = $opts{T};
  $c = 2;
}
die "Warning! You can only use -R or -T !\n" if ($opts{R} && $opts{T});

####
my @reads = <$readDir*_1_final.*> or die "Warning! No Input Libraries!\n";

my $indfilter = $refdir . "reconstructed_contigs/";
mkdir $indfilter unless -e $indfilter;

foreach my $read1 (@reads) {
  my $read2 = $read1; $read2 =~ s/_1_final/_2_final/;
  my $readU = $read1; $readU =~ s/_1_final/_u_final/;    
  my $lib = $1 if basename($read1) =~ m/(\S+)_1_final/;    
  my $Oritarget =  $refdir  . $lib . "_targetedRegionAndFlanking.fasta";
  my $codingbed =  $beddir  . $lib . "_coding.bed" if $coding ne "4" ;
  my $flankingbed =  $beddir  . $lib . "_flanking_ONLY.bed" if $coding ne "4" ;
  my $target_region = $beddir  . $lib . "_targeted_region_and_flanking.bed" if $coding eq "4" ;
  
  BuildContig ($lib, $refdir, $Oritarget, $read1, $read2, $readU, $vcf2fqMinDepth, $vcf2fqMaxDepth, $vcf2fqIndelFilter, $thread, $indfilter, $ma, $rep, $c, $opts{G},$opts{P}, $codingbed, $flankingbed, $insert,$Rlength, $con, $repeat,$cpu) if $coding ne "4" ;
  
  BuildContig ($lib, $refdir, $Oritarget, $read1, $read2, $readU, $vcf2fqMinDepth, $vcf2fqMaxDepth, $vcf2fqIndelFilter, $thread, $indfilter, $ma, $rep, $c, $opts{G},$opts{P}, $target_region,"0", $insert,$Rlength, $con,$repeat,$cpu) if $coding eq "4" ;
}


sub BuildContig {
  my ($name, $resdir, $template,$read1, $read2, $readCombined, $d, $D, $L, $p, $dir, $ma,$rep,$c, $GATK,$Picard, $tbed, $fbed, $insert, $Rlength, $con,$repeat,$cpu) = @_;
  my $template_fix = $template . "_final";
  open (TEM, "<",$template);
  open (TEMFINAL,">", $template_fix);
  
  my %temf;
  while (<TEM>) {
    chomp (my $line = $_);
    if ($line =~ m/^>(\S+)/) {
      chomp (my $seq = <TEM>);
      my $d = $1;
      $temf{$d}{'seq'} = $seq;
      my $length = length ($seq);
      my $seq2 = $seq;
      my $Ncount;
      $Ncount = ($seq2 =~ s/[N|n]//g);
      $temf{$d}{'N'} = $Ncount;	
      $temf{$d}{'len'} = $length;
      $temf{$d}{'ef'} =  $length - $Ncount;
      if ($Ncount) {
	$temf{$d}{'Nratio'} = $Ncount/$length;	
      }
      else {
	$temf{$d}{'Nratio'} = 0;
      }     
    }    
  }
  close TEM;
  
  foreach my $id (sort {$a cmp $b} keys %temf) {
    if ($temf{$id}{'Nratio'} <= 0.9) {
      print TEMFINAL ">", $id, "\n";    
      print TEMFINAL $temf{$id}{'seq'}, "\n";
    }
  }
  close TEMFINAL;
  
  system ("mv $template_fix $template");
  
  my $outPairedSam1 =  $resdir . $name . ".outPairedSam1";
  my $outSoloSam1 = $resdir . $name . ".outSoloSam1";
  my $paired_in_target =  $resdir . $name . ".paired_in_target.sam";
  my $solo_in_target = $resdir .$name . ".solo_in_target.sam";
  my $paired_in_target_bam = $resdir . $name .".paired_in_target.bam";
  my $solo_in_target_bam = $resdir . $name .".solo_in_target.bam";
  my $rawbam = $resdir . $name . ".raw.bam";
  my $sorted_in_target_bams = $resdir . $name . "_sorted.bam";
  
  my $indexref =  $template . ".nix";
  system ("novoindex $indexref  $template");
   
  my $std = $insert * 0.5;
  system("novoalign -R 30 -t 150  -d $indexref -f $read1  $read2  -i PE $insert, $std  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -n $Rlength -o SAM > $outPairedSam1");
  system("novoalign -R 30 -t 200 -d $indexref -f $readCombined -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -n $Rlength -o SAM > $outSoloSam1");


  if ($con == 1) {
    my $outPairedSam2 = $resdir . $name . ".outPairedSam2";
    open (IN, "<", $outPairedSam1);
    open (OUT, ">", $outPairedSam2 );
    
    while (<IN>) {
      chomp (my $line = $_);
      if ($line =~ /^@/) {
	print OUT $line, "\n";
      }
      else {
	  chomp (my @a = split /\s+/,$line);
	  if ( $a[1] == 99 || $a[1] == 147 || $a[1] == 83 || $a[1] == 163 ) {
	    print OUT $line, "\n";
	  } 
	}
    }
    close IN;
    close OUT;
    system ("mv $outPairedSam2 $outPairedSam1");
  }

  system("grep -v ZS:Z:NM $outPairedSam1 | grep -v NH:i: > $paired_in_target");
  system("grep -v ZS:Z:NM $outSoloSam1 | grep -v NH:i:  > $solo_in_target");  
  system("samtools view -bhS $paired_in_target >  $paired_in_target_bam");
  system("samtools view -bhS $solo_in_target > $solo_in_target_bam");
  system("samtools merge --reference $template $rawbam $paired_in_target_bam $solo_in_target_bam");
  system("samtools sort $rawbam --threads $cpu -o $sorted_in_target_bams");
  system("samtools index $sorted_in_target_bams"); 
  unlink ($outPairedSam1, $outSoloSam1, $paired_in_target, $solo_in_target, $paired_in_target_bam, $solo_in_target_bam, $rawbam); 
  system ("samtools faidx $template");
  
  my $refname = $1 if basename ($template) =~ m/(\S+)\.fasta/;
  my $index_ref = $resdir. $refname . '.dict'; 
  system ("picard -Xms512m -Xmx4g CreateSequenceDictionary R=$template O=$index_ref");  
  
  my $addgroupbam = $resdir . $name . ".rg.bam";
   
  system("picard -Xms512m -Xmx4g AddOrReplaceReadGroups INPUT=$sorted_in_target_bams OUTPUT=$addgroupbam RGID=$name RGLB=exonCap RGPL=illumina RGPU=lane1 RGSM=$name");
     
  my $intervals = $resdir . $name .".intervals";
  system("samtools index $addgroupbam");
  system ("gatk -Xms512m -Xmx4g -T RealignerTargetCreator -R $template -I $addgroupbam -o $intervals");
  my $final = $resdir . $name . ".bam2";
  my $sorted = $resdir . $name . "_sorted.bam";
  system ("gatk -Xms512m -Xmx4g -I $addgroupbam -R $template -T IndelRealigner -targetIntervals $intervals -o $final");
  unlink ($sorted_in_target_bams, $addgroupbam, $intervals);
  system ("rm $resdir$name*.bai $index_ref");
  system ("mv $final $sorted");
  system ("samtools index $sorted");
 
  my $fq_all = $resdir . $name .'_all.fq';
  
  print "\n\nNow reconstructing each locus!\n\n";
  
  system ("samtools mpileup -A -ugEf $template -Q 20 -q 2 $sorted | bcftools call -c - | vcfutils.pl vcf2fq -d $d -D $D -l $L >  $fq_all");
  
  
### filter out non-biallelic sites.  
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  my $loci_to_filter_target =  generateVCF ($template, $sorted , $tbed, $d) ;
  my $loci_to_filter_flanking =  generateVCF ($template, $sorted , $fbed, $d) if ($fbed ne "0");  
 
  my @filter_target = @{$loci_to_filter_target};
  my @filter_flanking = @{$loci_to_filter_flanking} if ($fbed ne "0");  

  sub generateVCF {
    my ($template, $sorted , $bed, $d) = @_;

    my @vcf = `samtools mpileup -A -t DP -I -ugEf $template $sorted -l $bed | bcftools call -c -`;
  
    my @loci_to_filter;
    foreach (@vcf) {
      chomp (my $line = $_);
      next if $line =~ /^#/;
      my $depth = $1 if  $line =~ /DP=(\d+)/;
      chomp (my @t = split /\s+/, $_);
      if ($t[3] ne 'N'){
	if ($t[4] =~ m/,/ && $depth >= $d) {
	  push @loci_to_filter, $t[0];	 
	}    
      }
    }   
    return (\@loci_to_filter);
  }

  ## note: number of sites in "loci_depth.txt" and ".unfiltered.fasta" is different because there are some loci do not have data so they don't show in "loci_depth.txt";
  
  my $fa_all = $dir . $name . '_all.fasta';

  system ("seqtk seq -q 20 -A $fq_all -n N > $fa_all");
 
  my $fa_all_fix = $dir . $name . '_fixLength_all.fasta';
  
  open (IN, "<", $fa_all);
  open (OUT, ">", $fa_all_fix);
  
  while (<IN>) {
    chomp (my $line = $_);
    if ($line =~ m/^>(\S+)/) {
      my $d = $1;
      chomp (my $seq = <IN>);
      my $ns = 'N' x ($temf{$d}{'len'} - length ($seq) );
      my $newseq = $seq . $ns;
      print OUT ">", $d , "\n";
      print OUT $newseq, "\n"; 
    }
  }
  close IN;
  close OUT;
  unlink ($fq_all, $fa_all);
  system ("mv $fa_all_fix $fa_all");
   
  process ($fa_all, $rep, $dir, $c, $name, $ma, \@filter_target,\@filter_flanking, $resdir, $sorted, $tbed, $fbed,\%temf, $repeat) if ($fbed ne "0");  
  process ($fa_all, $rep, $dir, $c, $name, $ma, \@filter_target,"0", $resdir, $sorted, $tbed, "0",\%temf,$repeat) if ($fbed eq "0"); 
  
  #system ("rm  $resdir*.fai  $resdir*sorted.bam* "); 
}


sub process {
  my ($fa, $rep,$dir,$c, $name, $ma, $filter_target, $filter_flanking, $resdir, $sorted, $tbed, $fbed,$refhash, $repeat) = @_;
  my @filter_target = @{$filter_target};
  my @filter_flanking = @{$filter_flanking} if $fbed ne "0";
  my %refhash = %{$refhash};  
  my $fa2 = $dir . $name . '_flanking_filtered.fasta' if $fbed ne "0";
  my $h =  $dir . $name . '_flanking_individual_H.txt' if $fbed ne "0";
  
  my $fa2_exon = $dir . $name . '_target_filtered.fasta';
  my $h_exon =  $dir . $name . '_target_individual_H.txt';
  
  my $exondir = $dir . "targetOnly/" ; 
  mkdir $exondir unless -e $exondir;
  
  my $alldir = $dir . "flanking/" if $fbed ne "0";
  mkdir $alldir if $fbed ne "0";

  my ($bed, $bedl, $all_length) = parsebed ($tbed);
  my ($bedf, $bedlf, $all_lengthf) = parsebed ($fbed) if $fbed ne "0";
  
  my %tbed = %{$bed};
  my %tbedl = %{$bedl};
  my %fbed = %{$bedf} if $fbed ne "0";
  my %fbedl = %{$bedlf} if $fbed ne "0";
     
  sub parsebed {
    my ($bed) = @_;
    open (BED, "<", $bed);
    my %bed;
    my %bedl;
    my $all_length;
    while (<BED>) {
      #JMSR003_indexing10_Contig1
      chomp (my @a = split /\s+/, $_);
      $bedl{$a[0]}{'ef'} += $a[2] - $a[1];
      $all_length += $a[2] - $a[1];
      push @{$bed{$a[0]}}, {'s' =>$a[1], 'e'=>$a[2]}; ;
      
    }
    close BED;
    return (\%bed, \%bedl, $all_length); 
  }

################################################################################################################################
  
  open (IN, "<", $fa);
  open (OUTF, ">", $fa2) if $fbed ne "0";
  open (OUTF2, ">", $h) if $fbed ne "0";
  open (OUTEXON, ">", $fa2_exon);
  open (OUTEXON2, ">", $h_exon);
  
  my $Ns = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN';
  
  while (<IN>) {
    chomp (my $line = $_);
    if ($line =~ /^>(\S+)/) {
      my $id = $1;
      chomp (my $seq = <IN>);
      
      if ($tbed{$id}) {
        my $Hexon = 0;
        my $full;
	$tbedl{$id}{'N'} = 0;
	
	for (my $i = 0; $i < scalar @{$tbed{$id}}; $i++) {
	  my $seqexon = substr $seq,  $tbed{$id}[$i]{'s'}, $tbed{$id}[$i]{'e'} - $tbed{$id}[$i]{'s'};
	  $seqexon =~ s/[a|c|g|t|n|y|s|k|r|w|m]/N/g ;
	  my $seq1 = $seqexon;
	  $Hexon += ($seq1 =~ s/[Y|S|K|R|W|M]//g);	 
	  if ($i == scalar @{$tbed{$id}} -1 ) {
	    $full .= $seqexon;
	  }
	  else {
	    $full .= $seqexon . $Ns;
	    $tbedl{$id}{'N'} += length ($Ns);
	  }	  
	}
	
	print OUTEXON  ">", $id, "\n"; 
	print OUTEXON  $full , "\n"; 
	print OUTEXON2 $id, "\t", sprintf("%.3f",$Hexon/$tbedl{$id}{'ef'}), "\n"; 
      }## if ($tbed{$id})
      
      if ($fbed ne "0" && $fbed{$id}) {
        my $Fexon = 0;
        my $Ffull;
	$fbedl{$id}{'N'} = 0;
	
	for (my $i = 0; $i < scalar @{$fbed{$id}}; $i++) {
	  my $seqexon = substr $seq,  $fbed{$id}[$i]{'s'}, $fbed{$id}[$i]{'e'} - $fbed{$id}[$i]{'s'};
	  $seqexon =~ s/[a|c|g|t|n|y|s|k|r|w|m]/N/g ;
	  my $seq1 = $seqexon;
	  $Fexon += ($seq1 =~ s/[Y|S|K|R|W|M]//g);	 
	  if ($i == scalar @{$fbed{$id}} -1 ) {
	    $Ffull .= $seqexon;
	  }
	  else {
	    $Ffull .= $seqexon . $Ns;
	    $fbedl{$id}{'N'} += length ($Ns);
	  }	  
	}
	print OUTF  ">", $id, "\n"; 
	print OUTF  $Ffull , "\n"; 
	print OUTF2 $id, "\t", sprintf("%.3f",$Fexon/$fbedl{$id}{'ef'}), "\n"; 
	
      } ##if ($filter_flanking ne "0" && $fbed{$id}) {    
    }##if ($line =~ /^>(\S+)/) {
  } ##while (<IN>) {
  
  close IN;
  close OUTF if $fbed ne "0";
  close OUTF2 if $fbed ne "0";
  close OUTEXON;
  close OUTEXON2;
  
  unlink ($fa);
  
  #my @filter_flanking = @{$filter_flanking} if $filter_flanking ne "0";
  if ( $fbed ne "0") {
    repeatANDdep ($fa2_exon,\%tbedl, 'exon',\%tbed, \@filter_target, $tbed, $rep, $c,$exondir, $ma, $sorted, $dir, $name, $h_exon, $repeat);
    repeatANDdep ($fa2, \%fbedl, 'flanking',\%fbed, \@filter_flanking, $fbed, $rep,$c,$alldir,$ma, $sorted, $dir, $name, $h, $repeat);
  }
  if ( $fbed eq "0") {
    repeatANDdep ($fa2_exon,\%tbedl, 'exon',\%tbed, \@filter_target, $tbed, $rep, $c,$exondir, $ma, $sorted, $dir, $name, $h_exon, $repeat);
  }
}

sub repeatANDdep {
  my ($fa2, $ref, $yes,$beds, $loci_to_filter, $bed, $rep,$c, $dir,$ma, $sorted, $oridir, $name, $h, $repeat) = @_;
 
  my %ref = %{$ref};
  my %bed = %{$beds};
  
  my @loci_to_filter = @{$loci_to_filter};
  
  my %filter_masked;
  if ($repeat == 1) {
    my $spath = $1 . 'dependencies/morePackages/'  if  (dirname(abs_path($0)) =~ /(\S+)scripts/) or die "scripts path is incorrect!\n";
    my $RepeatMasker = $spath ."RepeatMasker/RepeatMasker";
    system ("$RepeatMasker -q -no_is -noint -norna -species $rep $fa2") if ($c == 1);
    system ("$RepeatMasker -q -no_is -noint -norna  -lib $rep $fa2") if ($c == 2);
    
    my $infile = $fa2. ".masked";
    if (-z $infile ) {
      open (IN, "<", $infile);
      my $id;   
      while (<IN>) {
	chomp (my $line = $_);
	if ($line =~ /^>(\S+)/) {
	  $id = $1;
	}
	else {
	  $filter_masked{$id} .= $line;    
	}
      }
      close IN;
    }
    if (!-z $infile) {
      open (IN, "<",$fa2);
      my $id;   
      while (<IN>) {
	chomp (my $line = $_);
	if ($line =~ /^>(\S+)/) {
	  $id = $1;
	}
	else {
	  $filter_masked{$id} .= $line;    
	}
      }
      close IN;
    }
  }
  if ($repeat == 0) {
    open (IN, "<",$fa2);
    my $id;   
    while (<IN>) {
      chomp (my $line = $_);
      if ($line =~ /^>(\S+)/) {
	$id = $1;
      }
      else {
	$filter_masked{$id} .= $line;    	
      }
    }
    close IN;
  }
  
  my $fa_repeatmasked = $dir . $name . '_filtered.fasta';
  open (OUT, ">", $fa_repeatmasked);
  
  my $missing = $dir . $name . '_missing.txt';
  open (MISS, ">", $missing);   
  
  my $d;
  foreach my $id (sort {$a cmp $b} keys %filter_masked) {
    my $seq = $filter_masked{$id};
    my $length = length ($seq);
    my $seq1 = $seq;
    my $N = ($seq1 =~ s/N//ig) - $ref{$id}{'N'};
    
    if ($N/$ref{$id}{'ef'} >= $ma) {
      $d++;
    }
    unless ($N/$ref{$id}{'ef'} >= $ma){
      unless (grep {$_ eq $id} @loci_to_filter) {
	print OUT ">", $id, "\n", $seq, "\n";
	print MISS $id, "\t", $N/$ref{$id}{'ef'}, "\n";
      }
    }
  }
  close OUT;
  close MISS;
  print  $d, " Contigs were discarded from ", $fa2, " because of >= ", $ma, " missing data in their sequences!\n\n\n\n\n" if $d;
  print  "0 Contigs were discarded from ", $fa2, " because of >= ", $ma, " missing data in their sequences!\n\n\n\n\n" unless $d;
  system (" rm $fa2* ");
  
  my $hfinal = $dir . $name . "_individual_H.txt";
  system (" mv $h $hfinal");
  
  #####now calculate coverage#####
  my $file = $dir . $name . '.depth';
  my $per_gene_depth =  $dir . $name . "_loci_depth.txt";
  
  system ("samtools depth -b $bed $sorted  > $file");
      
  tie my @gene, 'Tie::Array::Packed::Number';
    
  open (IN, "<", $file);
  
  my %gene;
  while (<IN>) {
    chomp (my @l =split /\s+/, $_);
    $gene{$l[0]}{'dep'} += $l[2];
    $gene{$l[0]}{'count'} ++;
  }
  close IN;
  unlink ($file);
  
  open (OUT1, ">", $per_gene_depth);
  
  foreach my $g (sort { $a cmp $b} keys %gene) {
    print OUT1 $g, "\t",  $gene{$g}{'dep'}/$gene{$g}{'count'},"\n";
  }
  close OUT1;
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
