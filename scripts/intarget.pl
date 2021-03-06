#!/usr/bin/env perl

use warnings;
use strict;
use List::Util qw(max min);
use List::MoreUtils qw(uniq);
use Getopt::Std;
use File::Basename;
use File::Temp;
use List::Util qw(sum);
#no warnings 'recursion';


die(qq/

Usage: seqCapture intarget [options]

Basic options:

-t  FILE     Target sequence file in fasta format (.fasta or .fa)
-a  DIR      A folder with all final assemlies generated 
             by seqCapture assemble 
-m  FLOAT    How much should you cluster the targets and 
             assemblies at the get go [0.98]
-d  FLOAT    How much overlap is allowed between adjoining 
             assembled contigs mapping to the same target [0.3]
-p  INT      How similar does the assembled contig have to
             be to the target (note this is out of 100) [90]
-M  INT      Memory (in Mb) needed for cdhit [4096]
-T  INT      Number of threads used in blast and cdhit [10]
-E  FLOAT    Used in the initial BLAST step [1e-10]
-b  INT      Merging individual assemblies?
             1 = yes (for population genetic datasets)
             0 = no  (for phylogenetic datasets) [0]
-c  INT      Is the targeted loci from a mt genome (or most of it?)
             1 = yes  
             0 = no [0]    
-g  INT      For nuclear genes, retain flanking sequences or not
             1 = yes  
             0 = no [1]  
-L  INT      Min length cutoff in initial cdhit to keep 
             a assembled contig [200] 
-e  INT      Target sequences could be one of the following: 
             1=individual exons (one exon per gene); 
             2=cds (no UTR, sequence in coding frame); 
             3=transcripts (including UTR) or cDNA; 
             4=random (such as UCEs, no need for exon identification) 
             [3]          
-f  INT       +\/- flanking bp you would like to add [100]
-s  FILE      Annotated transcripts [required by e=1, 2, and 3]


\n\n/) unless (@ARGV);


my %opts = (t=>undef, a=>undef,  m=>0.98, d=>0.3, p=>90, E=>1e-10, e=>3, M=>4096, f=>100, s=>undef, b=>0, l=>0, c=>0, g=>1, T=>10,L=>200);
getopts('a:t:m:d:p:E:M:e:r:s:f:l:b:c:g:T:L:', \%opts);

my $mincut = $opts{L};
my $oriref = $opts{t};
my $cluster = $opts{m}; #how much should you cluster the targets and assemblies at the get go
my $maxOverlap =$opts{d}; #how much overlap is allowed between adjoining assembled contigs mapping to the same target
my $perMatch = $opts{p};#how similar does the assembled contig have to be to the target (note this is out of 100)
my $eval = $opts{E}; #used in the initial BLAST step
my $exon = $opts{e};
my $seqTrans1 = $opts{s} if $opts{s};
my $offset = $opts{f};
my $pop = $opts{b};
my $mem = $opts{M};
my $selfblast = $opts{l};
my $mt = $opts{c};
my $nu = $opts{g};
my $th = $opts{T};
my $dir = redir ($opts{a}); 

my @files = <$dir*.fasta> or die "can not find the merged assemblies!\n";

#make a result folder
my $resdir = $dir . "In_target/";
mkdir $resdir unless -e $resdir;

#rename the reference and make a file containing origianl names and the current names
my $seqfile = $oriref . "_rename";
my $corres = $oriref . "_rename_compared.txt";

#save the original reference sequence in a hash %ref
my %ref;
my $d = 0;
my $id;
open (IN, "<", $oriref);
while (<IN>) {
  chomp (my $line = $_);
  if ($line =~ m/^>(.*)/) {
    $id = $1;
    $d++;
  }
  else {
    $ref{$d}{'id'} = $id;
    $ref{$d}{'seq'} .= $line;
  }
}
close IN;

#rename the reference and make a file containing origianl names and the current names
open (OUT1, ">", $seqfile);
open (OUT2, ">", $corres);

foreach my $d (sort {$a <=> $b} keys %ref) {
  print OUT1 ">Contig" , $d, "\n";
  print OUT1 $ref{$d}{'seq'}, "\n";
  print OUT2 ">Contig" , $d, "\t", $ref{$d}{'id'},"\n";    
}
close OUT1;
close OUT2;

#final clustering the original reference
system ("cd-hit-est -i $seqfile -M $mem -o tmp  -c $cluster -T $th");
system ("mv tmp $seqfile");
system ("rm tmp*");

## open $corres and save names in a hash
my %protname;
my %trans;

if ($exon ne "4") {
  ##read target file.
  ##header musy be ENS\S+\d+)_\S+
  open (NAME, "<", $corres);
  
  while (<NAME>) {
    chomp (my $l = $_);
    chomp (my @line = split m/\s+/, $l);
    my $id = $1 if $line[0] =~ m/^>(\S+)/;
    my $prot = $1 if $l =~ m/\S+_(ENS\S+\d+)_\S+/;
    $protname{$prot} = $id if $prot;  
  } 
  close NAME;
  
  ##open transcript file
  open (TRANS, "<", $seqTrans1);
  
  my %trans2;
  my $match2;
  
  while(<TRANS>) {
    chomp(my $line = $_);
    next if $line =~ /^$/;
    if ($line =~ m/^>(\S+)/) {
      my $start2 = "";
      $start2 = $1 if $line =~ m/gs(\d+)_/;
      my $end2 = "";
      $end2 = $1 if $line =~ m/_ge(\d+)/;			
      $match2 = $1 if $line =~ m/(ENS\S+\d+)/;
      $trans2{$match2}{'start'} = $start2;
      $trans2{$match2}{'end'} = $end2;
    }
    else {
      $trans2{$match2}{'seq'} .= $line;
    }    
  }
  close TRANS;
  
  foreach my $match (keys %trans2) {      
    my $seq = $trans2{$match}{'seq'};
    
    if ($trans2{$match}{'start'} && $trans2{$match}{'end'} ) {
      my $length = $trans2{$match}{'end'} - $trans2{$match}{'start'} + 1;
      
      unless ($length % 3) {
	$length = 3 * int($length / 3);
      }    
      
      $seq = substr $seq, $trans2{$match}{'start'}  - 1, $length;
      my $aa = translate($seq);
      
      if ($protname{$match}) {
	$trans{$protname{$match}}{'protseq'} = $aa;
	$trans{$protname{$match}}{'transeq'} = $trans2{$match}{'seq'};
	$trans{$protname{$match}}{'protid'} = $match;
	$trans{$protname{$match}}{'u5'} = $trans2{$match}{'start'};
	$trans{$protname{$match}}{'u3'} = $trans2{$match}{'end'};    
	#print $protname{$match}, "\t", $match, "\n"; 
      }	  
    }
    
    else {
      my $embossin = $resdir . "embossIN.fa";
      my $embossout = "embossOUT.fa";
      open (OUT, ">", $embossin);
      print OUT ">" , $match, "\n";
      print OUT $seq, "\n";	
      close OUT;
      
      system ("getorf -sequence $embossin  -outseq $embossout -minsize 30 -maxsize 10000000 -find 2 -reverse N -auto -warning N -error N  -fatal N -die N");
      unless (-z $embossout) { 
	open (EMBOSS, "<", $embossout );
	unlink $embossout;
	unlink $embossin;
	my %emboss;
	my $id;
	my $number;
	while (<EMBOSS>) {
	  chomp (my $line = $_);
	  if ($line =~ m/^>(\S+)_(\d+)\s\[(\d+)\s\-\s(\d+)\]/) {
	    $id = $1;
	    $number = $2;
	    my $length = $4 - $3 + 1;
	    
	    $emboss{$id}{$number} = {'len' => $length};
	  }
	  else {
	    $emboss{$id}{$number}{'seq'} .= $line;
	  }
	}
	close EMBOSS;
	
	foreach my $number (sort {$emboss{$id}{$b}{'len'} <=> $emboss{$id}{$a}{'len'}} keys %{$emboss{$id}}) {
	  my $aa = translate ($emboss{$id}{$number}{'seq'});
	  if ($protname{$match}) {
	    $trans{$protname{$match}}{'protseq'} = $aa;
	    $trans{$protname{$match}}{'transeq'} = $emboss{$id}{$number}{'seq'};
	    $trans{$protname{$match}}{'protid'} = $match;
	    $trans{$protname{$match}}{'u5'} = 1;
	    $trans{$protname{$match}}{'u3'} = length ($emboss{$id}{$number}{'seq'});	      
	  }
	  last;
	}
      } ##unless (-z $embossout) { 
      if (-z $embossout) {
	print "can not find orf in $protname{$match} and $match ?!", "\n";
	exit;
      }
    } #else
  }     
}

#process each de novo assemblies  
my $sf = readSeq ($seqfile);
system("makeblastdb -in $seqfile -dbtype nucl");

my %refseq = %{$sf};
foreach (@files) {
  my $assembly = $_;
  my $lib = $1 if basename ($assembly) =~ m/(\S+)\.fa/;
  my $finalSeq = $resdir . $lib . "_intargetPremasked.fa";
  
  ##final clustering the de novo assemblies
  my $original =  $assembly . ".original";
  open (AS, "<", $assembly);
  open (OR, ">", $original);
  my $ddi = 1;
  while (<AS>) {
    chomp (my $l = $_);
    if ($l =~ /^>/) {
      my $id = "contig" . $ddi;
      print OR ">", $id, "\n";
    }
    else {
      print OR $l, "\n";
    }
    $ddi ++;
  }
  close OR;
  close AS;
  my $tmp = $original. ".tmp";
  system("cd-hit-est -i $original -M $mem -o $tmp -l $mincut -c $cluster -T $th");
  system("mv $tmp $assembly");
  system("rm $tmp*");
  
  my $assem = readSeq ($assembly);
  my %assem = %{$assem};
  
  #blast denovo assemblies and original target 
  my $assem_cut = $assembly. ".cut";
  blastcut (\%assem, $seqfile, $assembly, $assem_cut, $eval,  $resdir, $lib, $mem, $cluster, $th, $mincut);
  
  ###########new stuff here: first round blast###########
  
  my $outcut = $resdir . $lib. "_blast_cut.out";
  system ("blastn -query $assem_cut -db $seqfile -out $outcut -evalue $eval -outfmt 6 -num_threads $th");
  open(IN, "<$outcut");
  
  my %merge1;
  while (<IN>) {
    chomp(my $line = $_);
    my @d = split(/\s+/,$line);
    unless (grep {$_ eq $d[0]} @{$merge1{$d[1]}}) {	 
      push(@{$merge1{$d[1]}},$d[0]) if $d[2] >= $perMatch && $d[3] >= 50;
    }
    delete $merge1{$d[1]} if scalar (@{$merge1{$d[1]}}) == 0;
  }
  close(IN);
  system ("rm $outcut");
  #################done here####################################
  
  $assem = readSeq ($assem_cut);
  %assem = %{$assem};
  
  my $assembly1 = $assembly . ".in_target";
  open (OUT, ">", $assembly1);
  
  my $tracker = 1; 
  foreach my $gene (sort {$a cmp $b} keys %merge1) {    
    my $capin = $resdir . $lib . "_" . $gene . ".in";
    open (CAP, ">", $capin);
    for (my $i = 0; $i < scalar (@{$merge1{$gene}}); $i++) {
      print CAP ">", $merge1{$gene}[$i], "\n";
      print CAP $assem{$merge1{$gene}[$i]}, "\n";
    }
    close CAP;
    system ("cap3 $capin -o 16 -y 9999999999");
    
    my $assembled = $capin . ".cap.contigs";
    my $singlets = $capin . ".cap.singlets";
    
    my $id;
    
    if (! -z $singlets) {
      open(SIN, "<$singlets");
      while(<SIN>) {
	chomp(my $line = $_);
	if ($line =~ m/^>/) {
	  print OUT ">contig" . $tracker, "\n";
	  $tracker++;
	}
	else {
	  print OUT $line, "\n";
	}
      }
      close SIN;
    }
    
    
    if (! -z $assembled) {
      open(ASS, "<$assembled");
      while(<ASS>) {
	chomp(my $line = $_);
	if ($line =~ m/^>/) {
	  print OUT ">contig" . $tracker, "\n";
	  $tracker++;
	  }
	else {
	  print OUT $line, "\n";
	}
      }
      close ASS;
    }     
    system("rm $capin" . "*");	
  } 
  close OUT;
  
  my $out1 = $resdir . $lib. "_blastNEW.out1";
  system("blastn -query $assembly1 -db $seqfile -out $out1 -evalue $eval -outfmt 6 -max_target_seqs 1 -num_threads $th");
  
  if ($mt == 1 || $nu == 0) {
    my $assem2 = readSeq ($assembly1);
    my %assem2 = %{$assem2};
    my $assembly2 = $assembly1 . ".mt";
    open (OUT, ">", $assembly2);
    my $dd = 0;
    open(IN, "<", $out1);
    while (<IN>) {
      chomp(my $line = $_);
      $dd++;
      my @a = split(/\s+/,$line);
      print OUT ">contig", $dd, "\n";
      print OUT substr ($assem2{$a[0]},$a[6]-1, $a[7] - $a[6] +1), "\n"; 
    }
    close IN;
    close OUT;
    system ("mv $assembly2 $assembly1");
    system("blastn -query $assembly1 -db $seqfile -out $out1 -evalue $eval -outfmt 6 -max_target_seqs 1 -num_threads $th");
    mains($out1, $assembly1, $seqfile, $finalSeq, $resdir, $lib, \%trans, $exon,$eval,$offset,$pop,$perMatch,$maxOverlap, $selfblast);
  }
  
  else {
    mains($out1, $assembly1, $seqfile, $finalSeq, $resdir, $lib, \%trans, $exon,$eval,$offset,$pop,$perMatch,$maxOverlap, $selfblast);
  }    
}


if ($pop == 1) {   
  my $resdir = $dir . "In_target/";
  my @premask = <$resdir*_noChem.fasta>;
  my $precombined = $resdir . "combined_noChem.fasta";
  system ("cat $resdir*_noChem.fasta > $precombined");
  $pop = 0;
  $dir = $resdir;
  
  ###add stuff here
  my $ori =  $precombined . ".original";
  my $tmp =  $precombined . ".original.tmp";
  system ("cp $precombined $ori");
  system("cd-hit-est -i $precombined -M $mem -o $tmp -l 100 -c $cluster -T $th");
  system("mv $tmp $precombined");
  system("rm $tmp*");
  
  my %seqcontig;
  open (IN, "<", $precombined);
  my $ids;
  while (<IN>) {
    chomp (my $line = $_);
    if ($line =~ /^>(\S+)_\d+$/) {
      $ids = $1;
    }
    else {
      push @{$seqcontig{$ids}}, $line;
    }      
  }
  close IN;
  
  my $assembly3 = $resdir .  "combined_in_target.fasta";
  open (OUT1, ">", $assembly3);
  my $tracker = 1;
  
  foreach my $id (sort {$a cmp $b} keys %seqcontig) {
    my $seq = $resdir . $id . "_combined.seq.fasta";
    open (SEQS, ">", $seq);
    
    for (my $i = 0; $i < scalar (@{$seqcontig{$id}}); $i++) {
      print SEQS ">", $id, "_", $i+1, "\n";
      print SEQS $seqcontig{$id}[$i], "\n";
    }
    close SEQS;
    
    system ("cap3 $seq -o 16 -y 9999999999");
    
    my $assembled = $seq . ".cap.contigs";
    my $singlets = $seq . ".cap.singlets";
    
    if (! -z $singlets) {
      open(SIN, "<$singlets");
      while(<SIN>) {
	chomp(my $line = $_);
	if ($line =~ m/^>/) {
	  #	    print OUT ">contig" . $tracker, "\n";
	  print OUT1 ">contig" . $tracker, "\n";
	  $tracker++;
	}
	else {
	  #	    print OUT $line, "\n";
	  print OUT1 $line, "\n";
	}
      }
      close SIN;
    }
    
    if (! -z $assembled) {
      open(ASS, "<$assembled");
      while(<ASS>) {
	chomp(my $line = $_);
	if ($line =~ m/^>/) {
	  #	    print OUT ">contig" . $tracker, "\n";
	  print OUT1 ">contig" . $tracker, "\n";
	  $tracker++;
	}
	else {
	  #	    print OUT $line, "\n";
	  print OUT1 $line, "\n";
	}
      }
      close ASS;
    }
    
    system("rm $seq" . "*");
    
  } # foreach my $id (sort {$a cmp $b} keys %seqcontig) {
  
  close OUT1;
  system("makeblastdb -in $seqfile -dbtype nucl");
  my $out1 = $resdir . "combined_". "blastNEW.out1";
  system("blastn -query $assembly3 -db $seqfile -out $out1 -evalue $eval -outfmt 6 -max_target_seqs 1 -num_threads $th");      
  system("rm $seqfile.n* "); 
  
  my $lib = "combined";
  my $resdir1 = $resdir . "Final/";
  mkdir $resdir1 unless -e $resdir1;
  my $finalSeq = $resdir1 . "combined" . ".premask.fasta";
  mains($out1, $assembly3, $seqfile, $finalSeq, $resdir1, $lib, \%trans,  $exon,$eval,$offset,$pop,$perMatch,$maxOverlap,$selfblast);
}



sub blastcut {
  my ($assem, $seqfile, $assembly, $assem_cut, $eval, $resdir, $lib, $mem ,$cluster, $th,$mincut) = @_;
  
  my %assem = %{$assem};
  open (OUT, ">", $assem_cut);
  my $out = $resdir . $lib. "_blast.out";
  system ("blastn -query $assembly -db $seqfile -out $out -evalue $eval -outfmt 6 -num_threads $th");

  open(IN, "<$out");
  my %dup;
  while (<IN>) {
    chomp(my $line = $_);
    my @d = split(/\s+/,$line);
    push @{$dup{$d[0]}},$d[1]; 
  }
  seek IN, 0, 0;
  
  my %parse;
  foreach my $id (keys %dup) {
    my $length = scalar (uniq @{$dup{$id}}); 
    if ($length > 1) {
      $parse{$id}++;
    }
  }
  my $ned = 1;
  while (<IN>) {
    chomp(my $line = $_);
    my @d = split(/\s+/,$line);
    if ($parse{$d[0]}) {
      my $start = $d[6];
      my $end = $d[7];
      my $newseq = substr ($assem{$d[0]}, $start-1, $end -$start + 1);
      print OUT ">contig", $ned, "\n";
      print OUT $newseq,"\n";
      $ned++;
    }
    else {
      print OUT ">contig", $ned, "\n";
      print OUT $assem{$d[0]},"\n";
      $ned++;
    }
  }
  close IN;
  close OUT;
  my $tmp =   $assem_cut . "_blast_tmp.out";
  system("cd-hit-est -i $assem_cut -M $mem -o $tmp -l $mincut -c $cluster -T $th");
  system("mv $tmp $assem_cut");
  system("rm $tmp*");
  
  system ("rm $out"); ###remove the orginal blastout
  ########nee blast with trimmed contigs: second round#########
}


sub mains {
  my ($blastcon, $assembly1, $seqfile, $finalSeq, $resdir, $lib, $trans, $exon,$eval,$offset,$pop,$perMatch,$maxOverlap,$selfblast) = @_;
  my %trans = %{$trans};
  my %tmp;
  
  open (SEQA, ">", $finalSeq);
  
  open(IN, "<", $blastcon);
  while (<IN>) {
    chomp(my $line = $_);
    my @d = split(/\t/,$line);
    push(@{$tmp{$d[0]}},\@d) if $d[2] >= $perMatch && $d[3] >= 30; 
  }
  close(IN);
  
  my $contigs = readSeq($assembly1);
  my %contigs = %{$contigs};
  
  system ("rm $assembly1");
  system ("rm $blastcon");
  
  #parse the blast matches
  ### this part of code does not make sense any more. I keep it just for later use, just incase.
  my %matches;
  foreach my $id (sort {$a cmp $b} keys %tmp) {	
    #my $mArray = removeOverlap1($tmp{$id});
    my @mArray = @{$tmp{$id}};
    for (my $i = 0; $i < scalar(@mArray); $i++) {
      push(@{$matches{$mArray[$i][1]}}, \@{$mArray[$i]});
    }
  }	
  undef %tmp;	
   
  my %keep;
  foreach my $id (sort {$a cmp $b} keys %matches) {
    my ($mArray, $new, $overlap) = removeOverlap(\@{$matches{$id}}, \%contigs, \%keep, $maxOverlap);   
    
    $matches{$id} = $mArray;
    %keep = (%keep, %{$new});  
    
    if (@{$overlap}) {
      my @overlap = @{$overlap};
      
      @overlap = sort { $b->[0] <=> $a->[0]} @overlap; 
      
      for (my $i = 1; $i < (scalar @overlap); $i++) {
	my $ids = @{$overlap[$i]}[3];
	my $start;
	my $end;
	
	if (@{$overlap[$i]}[2] > @{$overlap[$i]}[1]) {
	  $start = @{$overlap[$i]}[1];
	  $end = @{$overlap[$i]}[2];	
	} 
	else {
	  $start =  @{$overlap[$i]}[2];
	  $end = @{$overlap[$i]}[1];	
	}
	my $seq = $contigs{$ids};
	
	substr ($seq, $start-1, $end-$start+1) = "N" x length (substr ($seq, $start-1, $end-$start+1)) ;
	#print $seq, "\n" if ($matches{$id}[0][1] eq "contig40");
	$contigs{$ids} = $seq;
      }
    }
  }
  
  #save sequences of original targets and de novo assemblies in hashes.
  my $seq = readSeq($seqfile); 
  my %seq = %{$seq};
  
  my %print;	
  
  foreach my $id (sort {$a cmp $b} keys %seq) {   
    if ($matches{$id}) {
      my %length;
      for (my $i = 0; $i < scalar(@{$matches{$id}}); $i++) {
	
	#next if $keep{$matches{$id}[$i][0]}; 
	my $start = $matches{$id}[$i][8];
	my $end = $matches{$id}[$i][9];
	for (my $n = min($start,$end); $n <= max($start,$end); $n++) {
	  $length{$n}++;
	}			
	$print{$matches{$id}[$i][0]}{$id}++;
      }
      my $overlap = sprintf("%.3f",scalar(keys %length) / length($seq{$id}));	
      #print ERR $id, "\t", $overlap, "\n";	
    }
    else {
      #print ERR $id, "\t", "NA\n";
    }
  }	
  
  ##print all the sequences in premasked.fa
  my %ids;	
  foreach my $c (sort {$a cmp $b} keys %print) {
    my $newid = join("_", keys %{$print{$c}}) . "_1";
    if ($ids{$newid}) {
      $ids{$newid}++;	
      if ($newid =~ m/(\S+)_(\d*)/) {
	my $core = $1;
	my $no = $ids{$newid};
	$newid  = $core . '_' . $no;
      }
    }
    else {
      $ids{$newid}++;
    }	
    print SEQA ">", $newid, "\n", $contigs{$c}, "\n"; 
  }
  #close ERR; 
  close SEQA;
  
  
  #now comes the interesting part...main function for looking for target sequences      
  Process ($finalSeq, $resdir, $lib, \%trans,$offset, $seqfile, $exon, $eval, $selfblast, $pop); 
  
  if ($pop == 1) {
    my $hash = seqhash($finalSeq);#no chimeric sequences included
    my %Nochem = %{$hash}; ## discard all chimeric sequences
    
    ## print non-chimeric sequences in $final1   
    my $final1 = $resdir . $lib . "_noChem.fasta";
    
    open (OUT , ">", $final1);
    foreach my $id (sort {$a cmp $b} keys %Nochem) {
      print OUT ">", $id, "\n";   
      print OUT $Nochem{$id}{'seq'}, "\n";
    }
    close OUT;
  }
}
#sub


sub annotateContigs {
  my ($contig, $refprot, $resdir, $lib, $contigid, $trans, $gene, $seq, $flag) = @_;
  #  $final3, $ref, $resdir, $name, $contig, \%trans, $trans{$seqid}{'protid'},\%seq
  my %bed;
  my %trans = %{$trans} unless $flag eq "5";
  my %seq = %{$seq};
  my @call;
  @call = `exonerate -m protein2genome -t $contig -q $refprot --showtargetgff TRUE --showalignment NO --showvulgar 0 -n 1 | grep 'cds' | sort -k 4,4 -n` unless $flag eq "5";

  @call = `exonerate -m est2genome -t $contig -q $refprot --showtargetgff TRUE --showalignment NO --showvulgar 0 -n 1 | grep 'exon'  | grep 'insertions' | sort -k 4,4 -n | head -1` if $flag eq "5";
  if (@call) {	 
	
    my (@gs, @ge);
    my $reverse = 0;
    my $forward = 0;
    foreach my $line (@call) {
   
      my @d = split(/\s+/,$line);
      my ($gs,$ge);
      if ($d[6] =~ m/-/) {
	my $length = length($seq{$contigid});
	$gs = $length - $d[4] + 1;
	$ge = $length - $d[3] + 1;	    
	push(@gs,$gs); push(@ge,$ge);
	$reverse ++;												
      }
      else {
	$gs = $d[3]; $ge = $d[4];
	push(@gs,$gs); push(@ge,$ge);
	$forward ++;
      }
    } ##foreach my $line (@call)	
    
    if ($reverse && ! $forward) {
      my $seqq = $seq{$contigid};
      $seqq = reverse($seqq);
      $seqq =~ tr/atgcATGC/tacgTACG/;
      $seq{$contigid} = $seqq; 
      @gs = reverse(@gs);
      @ge = reverse(@ge);
    } ##if ($reverse && !$foward)

    if ($reverse && $forward) {
      delete $seq{$contigid};
      print "$contigid is chimeric!!!\n";
    }
    
    unless ($reverse && $forward) {	
      for (my $x = 0; $x < scalar(@gs); $x++) {
	#print $gs[$x], "\t", $ge[$x], "\n";
 	my ($newstart, $newend) = refineorf ($gs[$x], $ge[$x], $seq{$contigid}, $resdir, $lib);
	
	push @{$bed{$contigid}},  [$newstart-1,$newend];
	
      }
    }    
  } ##if (@call)

  
  #a blast match but no exonerate match, check if it is that because it is a UTR ?	
  else {
    if ($flag eq "1" ||  $flag eq "2" || $flag eq "5") {
      delete $seq{$contigid};
      print "$contigid has no annotations!!!\n";
    }

    if ($flag eq "3") {
    my $seqid = $1 if $contigid =~ m/(\S+)_\d+/;

    my $seqTrans =  $resdir. $lib . "_" . $contigid .   "_trans.fa";
    open (COUT, ">", $seqTrans);
    print COUT ">", $trans{$seqid}{'protid'}, "\n", $trans{$seqid}{'transeq'} , "\n";
    close COUT;
        
    my $call = system("blat $seqTrans $contig blatTmp -out=blast8");
   
    open(IN, "<blatTmp");
    chomp(my $match = <IN>) unless (-z "blatTmp"); ###only take the first match
    
    if ($match) {
      my @d = split(/\s+/,$match);
      if ($d[10] <= 1e-5) {
	my $max;
	my $min;
	my $start;
	my $end;
	if ($d[8] < $d[9]) {
	  $max = $d[9];
	  $min = $d[8];
	  $start = $d[6];
	  $end = $d[7];
	} ##if ($d[8] < $d[9])
	else {
	  my $seqq = $seq{$contigid};
	  $seqq = reverse($seqq);
	  $seqq =~ tr/atgcATGC/tacgTACG/;
	  $seq{$contigid} = $seqq;  
	  $max = $d[8];
	  $min = $d[9];
	  my $len = length ($seqq);
	  
	  $end = $len - $d[6] + 1;
	  $start = $len - $d[7] + 1;
	} ##else

	
	my $overlap = 0;
	for (my $i = $min; $i <= $max; $i++) {
	  if ($i > $trans{$seqid}{'u5'} && $i < $trans{$seqid}{'u3'}) {
	    $overlap++;
	  }
	}	
	
	push @{$bed{$contigid}},  [0, 0] ;
	push @{$bed{$contigid}},  [$start-1, $end] ;

      }	## if ($d[10] <= 1e-5) {

      else {
	delete $seq{$contigid} ;
	print "$contigid has no annotations!!!\n";
      } ##else

    
    } ##if $match
    
    else {
      delete $seq{$contigid};
      print "$contigid has no annotations!!!\n";

    }
    
    system("rm blatTmp");
    unlink($seqTrans);
   } ## if ($flag eq "3") {
  } ###else 

  return(\%bed, \%seq);
}

  
 
sub Process {
   
  my ($new_master, $resdir, $name, $trans, $offset, $seqfile, $flag, $eval, $selfblast, $pop) = @_;

  ## $new_master is $finalseq-->premasked
  ## $seqfile is renamed reference
  ## $protein reference in $prot
  
  my %trans = %{$trans} unless $flag eq "4";
  
  ## discard chimeric sequences
  my $hash = seqhash($new_master);#no chimeric sequences included
  my %Nochem = %{$hash}; ## discard all chimeric sequences
  
  ## print non-chimeric sequences in $final1   
  my $final1 = $resdir . $name . "_noChem.fasta";
  
  open (OUT , ">", $final1);
  foreach my $id (sort {$a cmp $b} keys %Nochem) {
    print OUT ">", $id, "\n";   
    print OUT $Nochem{$id}{'seq'}, "\n";
  }
  close OUT;
  unlink ($new_master);
  
  ## save sequence of references in %ref. and print how many sequences in this file.
  my $dd;
  my %ref;
  my $idcontig;
  open (IN, "<", $seqfile);
  while (<IN>) {
    chomp (my $line = $_);  
    if ($line =~ m/^>(\S+)/) {
      $dd++;
      $idcontig = $1; 
    }
    else {
      $ref{$idcontig} .= $line;
    }    
  }
  close IN;
  #print "There are $dd sequences in the original target file!", "\n";
  
  ## blast premasked to itself to remove redundent part. 
  my $final2 = $final1 . "_copy";
  system ("cp $final1 $final2" );
  my $blastout = $resdir . $name .'.blast.out';
  my $call1 = system("makeblastdb -in $final1 -dbtype nucl > log");
  my $call2 = system("blastn -db $final1 -query $final2 -evalue 1e-20 -outfmt 6 -out $blastout");
  system("rm $final1.n* $final2 log");
  
  ## save sequences in premasked.fasta in %seq 
  my $hash1 = readSeq($final1);
  unlink ($final1);
  my %seq = %{$hash1};
  
  ## save blastout results in %tmp  
  my %tmp;
  open(IN, "<$blastout");
  while (<IN>) {
    chomp(my $line = $_);
    my @d = split(/\s+/,$line);
    if ($d[0] ne $d[1]) {
      push (@{$tmp{$d[0]}},\@d);
    }
    else {
      if ((max ($d[6], $d[7]) != max ($d[8], $d[9])) or (min ($d[6], $d[7]) != min ($d[8], $d[9])) ) {
	push (@{$tmp{$d[0]}},\@d); 
      }
    }
  }
  close(IN);
  system ("rm $blastout");
  
  
  #parse blastout results  
  foreach my $id (sort {$a cmp $b} keys %tmp) {
    if (scalar(@{$tmp{$id}}) > 1 ) {
      for (my $i = 0; $i < scalar(@{$tmp{$id}}); $i++) {  # from the second match. The first match is itself 
	my $start1;
	my $end1;
	my $start2; 
	my $end2; 
	next if $tmp{$id}[$i][0] eq $tmp{$id}[$i][1];
	if ($tmp{$id}[$i][6] < $tmp{$id}[$i][7]) {	    
	  $start1 = $tmp{$id}[$i][6];
	  $end1 = $tmp{$id}[$i][7];
	}
	if ($tmp{$id}[$i][6] > $tmp{$id}[$i][7]) {	    
	  $start1 = $tmp{$id}[$i][7];
	  $end1 = $tmp{$id}[$i][6];
	}
	
	if ($tmp{$id}[$i][8] < $tmp{$id}[$i][9]) {	    
	  $start2 = $tmp{$id}[$i][8];
	  $end2 = $tmp{$id}[$i][9];
	}
	if ($tmp{$id}[$i][8] > $tmp{$id}[$i][9]) {	    
	  $start2 = $tmp{$id}[$i][9];
	  $end2 = $tmp{$id}[$i][8];
	}
	
	my $seq1 = $seq{$tmp{$id}[$i][0]}  ;
	my $seq2 = $seq{$tmp{$id}[$i][1]}  ;
	
	substr ($seq1, $start1-1, $end1-$start1+1) = "n" x length (substr ($seq1, $start1-1, $end1-$start1+1)) ;
	substr ($seq2, $start2-1, $end2-$start2+1) = "n" x length (substr ($seq2, $start2-1, $end2-$start2+1)) ;
	
	$seq{$tmp{$id}[$i][0]} = $seq1 ;  
	$seq{$tmp{$id}[$i][1]} = $seq2 ;
	
	#print $seq1 , "\n" if  $tmp{$id}[$i][1] eq "Contig288_4_1";
	#print $seq2 , "\n" if  $tmp{$id}[$i][1] eq "Contig288_4_1";
      }
    }
  } ##  foreach my $id (sort {$a cmp $b} keys %tmp) {
  
  if ($selfblast == 0) {
    %seq = %{$hash1};
  }
  
  
  my %bed;
  my $dd1;
  
  
  foreach my $contig (sort {$a cmp $b} keys %seq) {
    
    my $yes = 0;
    my $seqid = $1 if $contig =~ m/(\S+)_\d+/;
    my $seq1 =  $seq{$contig}; 
    
    my $atgc = ($seq1 =~ s/[ATGCYWRMKS]//ig);
    
    delete $seq{$contig} if $atgc < 30; ### short sequence is no longer kept
    
    unless ( $flag  eq "4" ||  $flag  eq "5") {
      if ($trans{$seqid} && $seq{$contig} ) {
	##### if there is a match in the de novo assemblies, print that sequence
	my $final3 = $resdir . $name . "_" . $contig . "_target.fasta";    
	open (TAR, ">", $final3);
	print TAR ">",$contig, "\n";
	
	$dd1++;
	
	my $seq2 =  $seq{$contig};
	$seq2 = removeN ($seq2);
	$seq{$contig} = $seq2;
	
	print TAR $seq{$contig}, "\n";
	#print $seq{$contig}, "\n";
	### print the corresponding protein reference sequence in the target sequence file
	my $ref = $resdir . $name . "_" . $contig . "_ref.fasta";
	open (REF, ">", $ref);
	print REF ">", $trans{$seqid}{'protid'}, "\n";
	print REF $trans{$seqid}{'protseq'},"\n";
	
	close TAR;
	close REF;
	my $protid = $trans{$seqid}{'protid'}; 
	my ($newbed, $newseq) = annotateContigs ($final3, $ref, $resdir, $name, $contig, \%trans,$protid  ,\%seq, $flag);
	unlink ($final3, $ref);
	%seq = %{$newseq};
	%bed = (%bed, %{$newbed});      
      }
    } ##unless ( $flag  eq "4" ||  $flag  eq "5")
    
    
    if ( $flag  eq "5" ) {
      if ( $seq{$contig} ) {
	
	##### if there is a match in the de novo assemblies, print that sequence
	my $final3 = $resdir . $name . "_" . $contig . "_target.fasta";    
	open (TAR, ">", $final3);
	print TAR ">",$contig, "\n";
	
	my $seq2 =  $seq{$contig};
	#print $seq2 , "\n" if $contig eq "Contig288_4_1";
	$seq2 = removeN ($seq2);
	#print $seq2 , "\n" if $contig eq "Contig288_4_1";
	#exit if $contig eq "Contig288_4_1";
	$seq{$contig} = $seq2;
	
	print TAR $seq{$contig}, "\n";
	
	my $ref = $resdir . $name . "_" . $contig . "_ref.fasta";
	my $id2 = $1 if ($contig =~ m/(\S+_\d+)_\d+$/);
	open (REF, ">", $ref);
	print REF ">", $id2, "\n";
	print REF $ref{$id2},"\n";
	
	close TAR;
	close REF;
	
	my ($newbed, $newseq) = annotateContigs ($final3, $ref, $resdir, $name, $contig, "0", "0" ,\%seq, "5");
	unlink ($final3, $ref);
	%seq = %{$newseq};
	%bed = (%bed, %{$newbed});      
      }
    } ## unless ( $flag  eq "4" ) {
    
    
    if ($flag eq "4") {
      if ($seq{$contig} ) {
	my $final3 = $resdir . $name . "_" . $contig . "_target.fasta";    
	open (TAR, ">", $final3);
	print TAR ">",$contig, "\n";
	
	my $seq2 =  $seq{$contig};	
	$seq2 = removeN ($seq2);
	$seq{$contig} = $seq2;
	
	print TAR $seq{$contig}, "\n";
	#print $seq{$contig}, "\n";
	
	
	### print the corresponding reference sequence in the target sequence file
	my $ref2 = $resdir . $name . "_" . $contig . "_ref.fasta";
	my $id2 = $1 if ($contig =~ m/(\S+)_\d+/);
	open (REF, ">", $ref2);
	print REF ">", $id2, "\n";
	print REF $ref{$id2},"\n";
	
	close TAR;
	close REF;
	
	my $BEDblastout = $resdir . $name . "_" . $contig .'.blastBED.out';
	my $BEDblastout1 = $resdir . $name . "_" . $contig .'.blastBED.out.sorted';
	
	my $call5 = system("makeblastdb -in $final3 -dbtype nucl > log");
	my $call6 = system("blastn -db $final3 -query $ref2 -evalue $eval -outfmt 6 -out $BEDblastout ");
	system("rm $final3.n*  log");
	system ("sort -k 1,1 -k 7n,12 $BEDblastout > $BEDblastout1");
	
	open (IN, "<", $BEDblastout1);
	unlink ($final3, $ref2);
	unlink ($BEDblastout, $BEDblastout1);
	
	my $count =0;
	while (<IN>) {
	  $count++;		
	}	
	seek IN, 0, 0;
	##this is mostly becuase a lot of contig_XX are masked with ns.
	
	print $contig, " no blast match ???!!!! OK, $contig is removed from the marker list. " , "\n" if $count == 0;
	
	if ($count == 0) {
	  delete $seq{$contig};
	} ## if ($count == 0) {
	
	
	if ($count > 0) {
	  my $bedfile = $resdir . $name . "_" . $contig . ".tmpbed";
	  open (BED, ">", $bedfile);
	  my $reverse = 0;
	  my $forward = 0;
	  my @gs;
	  my @ge;
	  while (<IN>) {
	    my @d = split(/\s+/,$_);
	    my ($gs,$ge);
	    my $id = $d[1];   
	    if ($d[8] > $d[9] ) {
	      my $length = length($seq{$contig});
	      $gs = $length - $d[8] + 1;
	      $ge = $length - $d[9] + 1;    
	      push(@gs,$gs); push(@ge,$ge);
	      $reverse ++;												
	    }
	    else {
	      $gs = $d[8]; $ge = $d[9];
	      push(@gs,$gs); push(@ge,$ge);
	      $forward ++;
	    }
	  } ## while <IN>	
	  close IN;
	  if ($reverse && ! $forward) {
	    
	    my $seqq = $seq{$contig};
	    $seqq = reverse($seqq);
	    $seqq =~ tr/atgcATGC/tacgTACG/;
	    $seq{$contig} = $seqq; 
	    #@gs = reverse(@gs);
	    #@ge = reverse(@ge);
	  } ##if ($reverse && !$foward)
	  
	  if ($reverse && $forward) {
	    delete $seq{$contig};
	    #print "$contig is chimeric and removed from dataset!!!\n";
	  }
	  
	  unless ($reverse && $forward) {	
	    for (my $x = 0; $x < scalar(@gs); $x++) {              
	      print BED $contig, "\t", $gs[$x]-1, "\t", $ge[$x], "\n";
	    }
	  }
	  close BED;
	  if (!-z $bedfile) {
	    MakeBed ($bedfile, 0) ;
	    open (IN2, "<", $bedfile);
	    
	    while (<IN2>) {
	      chomp (my @line = split /\s+/, $_);
	      push @{$bed{$line[0]}},  [$line[1] ,$line[2]];
	      
	    }
	    close IN2;
	  }
	  unlink ($bedfile);
	} ##if ($count > 0) {  
      }
    } ##if ($flag eq "4") {
    
  } ## foreach my $contig (sort {$a cmp $b} keys %seq)
  
  
  my ($exonflankingbed, $flankingbed ,$seq1, $bed2) = modbed (\%bed, $resdir, $name, $offset, \%seq); 
  
  my %exonflankingbed = %{$exonflankingbed} ;
  
  my %flankingbed= %{$flankingbed};
  
  %seq = %{$seq1};
  
  %bed = %{$bed2};
  
  #### remove utr from %bed if there is utr present;
  
  unless ($flag eq "4") {
    foreach my $contig (sort {$a cmp $b} keys %bed) {
      if ($bed{$contig}[0][1] == 0) {
	my $start = $exonflankingbed{$contig}[0][0];
	my $end = $exonflankingbed{$contig}[0][1];
	delete $flankingbed{$contig};
	$flankingbed{$contig}[0][0] = $start;
	$flankingbed{$contig}[0][1] = $end;
	delete $bed{$contig};
      }
    }
  }
  
  
  my $Ns = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'; ##39 Ns
  my %position;
  my %everything;
  
  
  foreach my $contig (sort {$a cmp $b} keys %ref) {
    
    my $final3 = $resdir . $name . "_" . $contig . "_target.fasta";    
    my $ref = $resdir . $name . "_" . $contig . "_ref.fasta";
    
    open (TAR, ">", $final3);
    
    open (REF, ">", $ref);
    print REF ">", $contig, "\n";
    print REF $ref{$contig},"\n";
    
    my $yes = 0;
    my $d;
    foreach my $id (sort {$a cmp $b} keys %seq) { 
      
      my $seqid = $1 if $id =~ m/(\S+)_\d+$/;
      
      if ($contig eq $seqid) { 
	
	print TAR ">",$id, "\n";
	print TAR $seq{$id}, "\n";
	$yes++;
      }
    }
    close REF;
    close TAR;
    
    unlink ($final3, $ref) if $yes == 0;
    
    
    if ($yes > 0) { ### if there is a match
      
      my $BEDblastout = $resdir . $name . "_" . $contig .'.blastBED.out';
      my $BEDblastout1 = $resdir . $name . "_" . $contig .'.blastBED.out.sorted';
      
      my $call5 = system("makeblastdb -in $final3 -dbtype nucl > log");
      my $call6 = system("blastn -db $final3 -query $ref -evalue $eval -outfmt 6 -out $BEDblastout ");
      system("rm $final3.n*  log");
      system ("sort -k 1,1 -k 7n,12 $BEDblastout > $BEDblastout1");
      
      open (IN, "<", $BEDblastout1);
      unlink ($final3, $ref);
      unlink ($BEDblastout, $BEDblastout1);
      
      my $count =0;
      while (<IN>) {
	$count++;		
      }	
      seek IN, 0, 0;
      
      if ($count == 0) {
	foreach my $id (sort {$a cmp $b} keys %seq) { 
	  my $seqid = $1 if $id =~ m/(\S+)_\d+$/;
	  if ($contig eq $seqid) { 
	    delete $seq{$id};	    
	  }
	}
      } ## if ($count == 0) {
      print $contig, " has no blast match???!!!! OK, $contig is removed from the marker list. " , "\n" if $count == 0;
      
      if ($count > 0) {
	
        my $ds = 0;
	my $addlength;
	my $previous = 0;
	my @match;
	
	while (<IN>) {
	  my @a = split /\s+/,$_;
	  my $refContig = $a[0];
	  my $assContig = $a[1];
	  
  	  my $answer = match ($assContig,\@match);
	  
	  push @match, $assContig unless ($answer eq 'yes') ;
	  $everything{$refContig}{'seq'} .=  $seq{$assContig} . $Ns unless ($answer eq 'yes' ); ## when $assContig is a new contig in the array
	  
          $previous = length ($seq{$assContig} . $Ns) unless ($answer eq 'yes'); ## when $assContig is a new contig in the array	  
          
	  # print $assContig , "\t", $previous, "\n" unless ($answer eq 'yes');
	  
	  $addlength = length ($everything{$refContig}{'seq'}) - $previous unless ($answer eq 'yes'); ## when $assContig is a new contig in the array
          
	  
	  unless ($answer eq 'yes') {
            if ($bed{$assContig}){
	      for (my $i = 0; $i < scalar @{$bed{$assContig}}; $i++) {
		#print $assContig, "\t", $bed{$assContig}[$i][0], "\t", $bed{$assContig}[$i][1], "\n";
		$bed{$assContig}[$i][0] = $bed{$assContig}[$i][0] + $addlength  ;
		$bed{$assContig}[$i][1] = $bed{$assContig}[$i][1] + $addlength ;
		#print $assContig, "\t", $bed{$assContig}[$i][0], "\t", $bed{$assContig}[$i][1], "\n";
	      }
	    }
            if ($exonflankingbed{$assContig}) {
	      for (my $i = 0; $i < scalar @{$exonflankingbed{$assContig}}; $i++) {
		#print $assContig, "\t", $exonflankingbed{$assContig}[$i][0], "\t", $exonflankingbed{$assContig}[$i][1], "\n";
		$exonflankingbed{$assContig}[$i][0] = $exonflankingbed{$assContig}[$i][0] + $addlength ;
		$exonflankingbed{$assContig}[$i][1] = $exonflankingbed{$assContig}[$i][1] + $addlength ;	
		#print $assContig, "\t", $exonflankingbed{$assContig}[$i][0], "\t", $exonflankingbed{$assContig}[$i][1], "\n";      
	      }
	    }
	    
	    if ($flankingbed{$assContig}) {
	      for (my $i = 0; $i < scalar @{$flankingbed{$assContig}}; $i++) {
                #print $assContig, "\t", $flankingbed{$assContig}[$i][0], "\t", $flankingbed{$assContig}[$i][1], "\n";
		$flankingbed{$assContig}[$i][0] = $flankingbed{$assContig}[$i][0] + $addlength ;
		$flankingbed{$assContig}[$i][1] = $flankingbed{$assContig}[$i][1] + $addlength ;
                #print $assContig, "\t", $flankingbed{$assContig}[$i][0], "\t", $flankingbed{$assContig}[$i][1], "\n";	      
	      }
	    }
	    $ds++;
	  }
	  
	} ## while (<IN>)
	
	close IN;	
	$everything{$contig}{'seq'} = substr($everything{$contig}{'seq'}, 0, length($everything{$contig}{'seq'})-length ($Ns));
	my $final4 = $resdir . $name . "_" . $contig . "_target.fasta1";
	my $ref4 = $resdir . $name . "_" . $contig ."_ref.fasta1";
	
	open (TAR, ">", $final4);
	open (REF, ">", $ref4);
	
	print REF ">", $contig, "\n";
	print REF $ref{$contig},"\n";
	
	print TAR ">",$contig . "_ass", "\n";
	print TAR $everything{$contig}{'seq'}, "\n";
	
	close REF;
	close TAR;
	
	my $BEDblastout3 = $resdir . $contig .'.blastBED.out1';
	my $BEDblastout4 = $resdir .  $contig .'.blastBED.out.sorted1';
	system("makeblastdb -in $final4 -dbtype nucl > log");
	system("blastn -db $final4 -query $ref4 -evalue $eval -outfmt 6 -out $BEDblastout3 ");
	system("rm $final4.n*  log");
	system ("sort -k 1,1 -k 7n,12 $BEDblastout3 > $BEDblastout4");
	
	open (IN2, "<", $BEDblastout4);
	unlink ($BEDblastout3, $BEDblastout4);
	unlink ( $final4, $ref4); 
	
	while (<IN2>) {
	  chomp (my @a = split /\s+/, $_);
	  my $refContig = $a[0];
	  my $start; 
	  my $end;
	  if ($a[8]  < $a[9]) { 
	    $start = $a[8];
	    $end = $a[9];
	  }
	  if ($a[8]  > $a[9]) { 
	    $start = $a[9];
	    $end = $a[8];
	  }
	  print $contig, " has inverted coordinates !\n" if ($a[8]  > $a[9]);
          
	  push @{$position{$refContig}},[$start , $end];
	}
	close IN2;
	
      } ## if ($count > 0) {
    } ##if ($yes > 0)
  } ##foreach my $contig (sort {$a cmp $b} keys %ref) {	
  
  ###print all different output files. 
  my $exonbed;
  $exonbed  = $resdir . $name ."_coding.bed" unless $flag eq "4";
  $exonbed = $resdir . $name ."_targeted_region.bed" if $flag eq "4";
  
  
  open (EXON, ">", $exonbed);
  foreach my $contig (sort {$a cmp $b} keys %bed) {
    for (my $i = 0; $i < scalar @{$bed{$contig}}; $i++) {
      my $seqig = $1 if $contig =~ /(\S+)_\d+$/;
      print EXON  $name, "_",$seqig, "\t", $bed{$contig}[$i][0], "\t", $bed{$contig}[$i][1], "\n" if $seq{$contig};
    }
  }
  
  close EXON;
  
  my $exonAndflankingbed;
  $exonAndflankingbed = $resdir . $name . "_coding_and_flanking.bed" unless $flag eq "4";
  $exonAndflankingbed = $resdir . $name . "_targeted_region_and_flanking.bed" if $flag eq "4";
  
  open (EXONANDFLANK, ">", $exonAndflankingbed);
  
  foreach my $contig (sort {$a cmp $b} keys %exonflankingbed) {
    for (my $i = 0; $i < scalar @{$exonflankingbed{$contig}}; $i++) {
      my $seqig = $1 if $contig =~ /(\S+)_\d+$/;
      print EXONANDFLANK $name, "_", $seqig, "\t", $exonflankingbed{$contig}[$i][0], "\t", $exonflankingbed{$contig}[$i][1], "\n" if $seq{$contig};
    }
  }
  
  close EXONANDFLANK;
  
  my $flankingbed1;
  $flankingbed1 = $resdir . $name . "_flanking_ONLY.bed" unless $flag eq "4";
  $flankingbed1 = $resdir . $name . "_flanking_ONLY.bed" if $flag eq "4";
  
  
  open (FLANKONLY, ">", $flankingbed1);
  
  foreach my $contig (sort {$a cmp $b} keys %flankingbed) {
    for (my $i = 0; $i < scalar @{$flankingbed{$contig}}; $i++) {
      my $seqig = $1 if $contig =~ /(\S+)_\d+$/;
      print FLANKONLY  $name, "_", $seqig, "\t", $flankingbed{$contig}[$i][0], "\t", $flankingbed{$contig}[$i][1], "\n" if $seq{$contig};
    }
  }
  
  close FLANKONLY;
  
  my $AllBed  =  $resdir . $name .'_targetedRegionforExonCapEval.bed';
  open (BED, ">", $AllBed);
  
  foreach my $item (sort {$a cmp $b} keys %position) {
    for (my $i = 0 ; $i < scalar @{$position{$item}}; $i++ ) {
      print BED $name, "_", $item, "\t";
      print BED $position{$item}[$i][0]-1, "\t", $position{$item}[$i][1],"\n";
    }
  }
  close BED;
  
  MakeBed ($AllBed, 0);
  
  my $contigbed =  $resdir . $name .'_allcontig.bed';
  open (ALLBED2, ">", $contigbed  );
  
  my $finalfasta  =  $resdir . $name .'_targetedRegionAndFlanking.fasta';
  open (FASTA, ">", $finalfasta);
  
  #$everything{$contig}{'seq'}
  foreach my $item (sort {$a cmp $b} keys %everything) {   
    print FASTA ">", $name, "_", $item, "\n";
    print FASTA $everything{$item}{'seq'}, "\n"; 
    print ALLBED2   $name, "_", $item , "\t", "0", "\t", length ($everything{$item}{'seq'}), "\n";
  }
  close FASTA ;
  close ALLBED2 ;
  
  ### perform a final selfblast
  my $finalfasta2 = $finalfasta . "_copy";
  system ("cp $finalfasta $finalfasta2" );
  my $blastout1 = $resdir . $name .'.final.blast.out';
  system("makeblastdb -in $finalfasta2 -dbtype nucl > log");
  system("blastn -db $finalfasta2 -query $finalfasta -outfmt 6 -evalue 1e-10 -out $blastout1");
  system("rm $finalfasta2.n* $finalfasta2 log");
  
  my %rmsites;
  open(IN, "<$blastout1");
  while (<IN>) {
    chomp(my $line = $_);
    my @d = split(/\s+/,$line);
    if ($d[0] ne $d[1]) {
      #if ($d[2] > 85 || $d[2] < 100 ) {
      push @{$rmsites{$d[0]}}, {'start' => $d[6], 'end' => $d[7]};
      push @{$rmsites{$d[1]}}, {'start' => $d[8], 'end' => $d[9]} if $d[9] > $d[8] ;
      push @{$rmsites{$d[1]}}, {'start' => $d[9], 'end' => $d[8]} if $d[8] > $d[9] ; 
      #}
    }
    else {
      if ((max ($d[6], $d[7]) != max ($d[8], $d[9])) || (min ($d[6], $d[7]) != min ($d[8], $d[9]))) {
	push @{$rmsites{$d[0]}}, {'start' => $d[6], 'end' => $d[7]};
	push @{$rmsites{$d[1]}}, {'start' => $d[8], 'end' => $d[9]} if $d[9] > $d[8] ;
	push @{$rmsites{$d[1]}}, {'start' => $d[9], 'end' => $d[8]} if $d[8] > $d[9] ; 
      }
    }
  }
  close(IN);
  system ("rm $blastout1");
  
  
  my $sites_to_remove = $resdir . $name . "_sites_to_remove.txt";
  open (RM, ">", $sites_to_remove);
  print RM "chr\tsite\n";
  my %s;
  foreach my $id (sort {$a cmp $b} keys %rmsites) {
    for (my $i = 0; $i < scalar (@{$rmsites{$id}}); $i++) {
      for (my $j = $rmsites{$id}[$i]{'start'}; $j <= $rmsites{$id}[$i]{'end'}; $j++) {
	if (!$s{$id}{$j}) {
	  print RM $id, "\t", $j, "\n";
	  $s{$id}{$j}++;
	}
      }
    }
  }
  close RM;
  
  my $sites_to_remove1 = $sites_to_remove . ".sorted";
  system ("sort -k1,1 -k2,2n $sites_to_remove > $sites_to_remove1");
  system ("mv $sites_to_remove1 $sites_to_remove");  
  
  
  my $finalfadir = $resdir . "fasta/";
  my $finalbeddir =  $resdir . "bed/";
  mkdir $finalfadir unless -e $finalfadir;
  mkdir $finalbeddir unless -e $finalbeddir ;
  system ("mv $resdir*.fasta $finalfadir");
  system ("mv $resdir*.bed $resdir*.txt  $finalbeddir");
} #Process

sub removeOverlap1 {
  my ($array) = @_;  
  for (my $i = 0; $i < scalar(@$array); $i++) {
    my $start1 = $array->[$i]->[6];
    my $end1 = $array->[$i]->[7];
    my %bp;
    for (my $n = min($start1,$end1); $n <= max($start1,$end1); $n++) {
      $bp{$n}++;
    }
    for (my $j = $i+1; $j < scalar(@$array); $j++) { 	
      my $start2 = $array->[$j]->[6];
      my $end2 = $array->[$j]->[7];
      my $overlap = 0;
      for (my $n = min($start2,$end2); $n <= max($start2,$end2); $n++) {
	$overlap++ if $bp{$n};
      }
      $overlap = $overlap / min(abs($start1 - $end1),abs($start2 - $end2));	
      if ($overlap > 0.2) {
        #undef @{$array}; 	
      } ##if ($overlap > 0.2) {
    }##for (my $j = $i+1; $j < scalar(@$array); $j++) { 	
  } ##for (my $i = 0; $i < scalar(@$array); $i++) {	
  return($array);	
}

sub readFile {
  my ($file,$hash,$base) = @_;
  if (-f $file) {
    open(IN, "<$file");	
    my $id; my $tracker = 1;
    while(<IN>) {
      chomp(my $line = $_);
      if ($line =~ m/>(\S+)/) {
	$id = $base . $tracker;
	$tracker++;
      }
      else {
	$hash->{$id} .= $line;
      }
    }
    close(IN);	
  }	
  return($hash);
}	

sub readSeq {
  my ($seqfile) = @_;
  my %seq; my $id;
  open(IN, "<$seqfile");
  while(<IN>) {
    chomp(my $line = $_);
    if ($line =~ m/>(\S+)/) {
      $id = $1;
    }
    else {
      $line =~ s/\r//g;
      $seq{$id} .= $line;
    }
  }
  close(IN);
  return(\%seq);
}


sub seqhash {
  my ($file) = @_;
  my %seq;
  open (IN, "<", $file); 
  my $id;
  my $d = 0;
  while (<IN>) {
    chomp (my $line = $_);
    if ($line =~ m /^>(Contig\d+_\d+)$/) {   
      $id = $1;
      chomp (my $seq = <IN>);
      $seq{$id}{'seq'} = $seq;
    }
    if ($line =~ m /^>Contig\d+_Contig\d+/) {
      print $line, "\n";
      $d++;    
    }  
    else {
    next;
    }
  }
  print "\n\nIn ", $file, ", " ,  $d, " chimeric sequences are discarded!", "\n\n";
  close IN;
  return (\%seq);
}


sub MakeBed {
  my ($bed, $offset) = @_;

  open (BED, "<", $bed);
  my @Bed;
  my $count =0;
  while (<BED>) {
    chomp (my @name = split /\s+/, $_);
    #add bed to matrix
    
      $Bed[$count] -> [0] = $name[0]; #contig name
      $Bed[$count] -> [1] = $name[1]; #start coordinate
      $Bed[$count] -> [2] = $name[2]; #end coordinate
      $count++;
    
  }
  close BED;
  
  @Bed = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]} @Bed; 
  my ($tmp_offset, $tmp);

  #create temp files to save space!  
  $tmp_offset = File::Temp->new(SUFFIX=>'.offset');
  $tmp = File::Temp->new(SUFFIX=>'.original');

  foreach (@Bed) {
    print $tmp @{$_}[0], "\t", @{$_}[1], "\t", @{$_}[2], "\n";   
    if (@{$_}[1] >= $offset) {
      print $tmp_offset @{$_}[0], "\t", @{$_}[1]-$offset, "\t", @{$_}[2]+$offset, "\n";
    }
    else {
      print $tmp_offset @{$_}[0], "\t", 0, "\t",  @{$_}[2]+$offset, "\n";
    }
  }
  
  seek $tmp, 0, 0; #rewind $region list to beginning for main loop
  seek $tmp_offset, 0, 0;
  
  my $out1 = $bed . "_offset_final";
  mining ($tmp_offset, $out1);
  my $out2 = $bed . "_final";
  mining ($tmp, $out2);
  unlink ($tmp_offset, $tmp);
  system ("mv $out2 $bed");
  system ("rm $out1");
}

sub mining {
  my ($file, $out) = @_;
  open (IN, "<", $file);
  my $howmany = 0;
  while (<IN>) {
    $howmany++;
  }
  close IN;
  
  if ($howmany > 1) {
    open (IN, "<", $file);
    open (OUT, ">", $out);
    chomp (my $first = <IN>);
    my @d = split /\s+/, $first;
    
    my $gene = $d[0];
    my $start = $d[1];
    my $end = $d[2];
    
    while (<IN>) { ##do not forget to process the end of file!!
      
      chomp (my @line = split /\s+/, $_);
      if ($line[0] eq $gene) {
	
	if (!eof){ ## if not the end of the file
	  
	  if ($line[1] > $end) {
	    print OUT $gene, "\t", $start, "\t", $end, "\n";
	    $start = $line[1];
	    $end = $line[2];
	  } 	
	  if ($line[1] <= $end && $line[1] >= $start) {
	    if ($line[2] >= $end) {
	      $end = $line[2];	  
	    }	 	  
	  }
	  if ($line[1] < $start) {
	    $start = $line[1];
	    if ($line[2] >= $end) {
	      $end = $line[2];	  
	    }	  
	  }
	} ##if (!eof){
	
	else { #end of the file; need to print both lines. Same as below
	  if ($line[1] > $end) {
	    print OUT $gene, "\t", $start, "\t", $end, "\n";
	    print OUT $gene, "\t", $line[1], "\t", $line[2], "\n";
	  } 	
	  if ($line[1] <= $end && $line[1] >= $start) {
	    if ($line[2] >= $end) {
	      print OUT $gene, "\t", $start, "\t", $line[2], "\n";	      
	    } 
	    else {
	      print OUT $gene, "\t", $start, "\t", $end, "\n";
	    }
	  }
	  if ($line[1] < $start) {
	    if ($line[2] >= $end) {
	      print OUT $gene, "\t", $line[1], "\t", $line[2], "\n";
	    }
	    else {
	      print OUT $gene, "\t", $line[1], "\t", $end, "\n";
	    } 
	  }
	} ## else	
      } ##if ($line[0] eq $gene) {
      
      
      if ($line[0] ne $gene) {
	if (!eof) {
	  print OUT $gene, "\t", $start, "\t", $end, "\n";	
	  $gene = $line[0];
	  $start = $line[1];
	  $end = $line[2];   
	}
	else {
	  print OUT $gene, "\t", $start, "\t", $end, "\n";
	  $gene = $line[0];
	  $start = $line[1];
	  $end = $line[2]; 
	  print OUT $gene, "\t", $start, "\t", $end, "\n";
	}
      
       } ##if ($line[0] ne $gene) {
    } ## while (<IN>) {
    close IN;
    close OUT;

    } ## if $howmany >1
    elsif ($howmany  == 1)  {
      open (IN, "<", $file);
      open (OUT, ">", $out);
      chomp (my $first = <IN>);
      my @d = split /\s+/, $first;
      
      my $gene = $d[0];
      my $start = $d[1];
      my $end = $d[2];
      
      print OUT $gene, "\t", $start, "\t", $end ,"\n";
      close IN;
      close OUT;
      
    } ##elsif
    
    else {
      print $file, "\n";
      #exit; 
    }
  } 
  
	
sub translate {
  my $string = shift;
  $string = uc($string);
  my @codons = $string =~ m/(\S\S\S)/g;
  my %codons = (	'ATG'=>'M','ACG'=>'T','CTG'=>'L','CCG'=>'P','GTG'=>'V','GCG'=>'A','TTG'=>'L','TCG'=>'S',
			'ATA'=>'I','ACA'=>'T','CTA'=>'L','CCA'=>'P','GTA'=>'V','GCA'=>'A','TTA'=>'L','TCA'=>'S',
			'ATC'=>'I','ACC'=>'T','CTC'=>'L','CCC'=>'P','GTC'=>'V','GCC'=>'A','TTC'=>'F','TCC'=>'S',
			'ATT'=>'I','ACT'=>'T','CTT'=>'L','CCT'=>'P','GTT'=>'V','GCT'=>'A','TTT'=>'F','TCT'=>'S',
			'AGG'=>'R','AAG'=>'K','CGG'=>'R','CAG'=>'Q','GGG'=>'G','GAG'=>'E','TGG'=>'W','TAG'=>'*',
			'AGA'=>'R','AAA'=>'K','CGA'=>'R','CAA'=>'Q','GGA'=>'G','GAA'=>'E','TGA'=>'*','TAA'=>'*',
			'AGC'=>'S','AAC'=>'N','CGC'=>'R','CAC'=>'H','GGC'=>'G','GAC'=>'D','TGC'=>'C','TAC'=>'Y',
			'AGT'=>'S','AAT'=>'N','CGT'=>'R','CAT'=>'H','GGT'=>'G','GAT'=>'D','TGT'=>'C','TAT'=>'Y');
  my $translate;
  foreach(@codons) {
    if ($codons{$_}) {
      $translate = $translate . $codons{$_};
    }
    else {
      #print "ERROR: ILLEGAL PASS TO CODON TRANSLATION: $_ is not a codon!\n";
      $translate = $translate . 'X';
    }
  }
  return($translate);
}


sub blast {
  my ($exome,$prot, $eval2, $dir) = @_; 
  my $out = $dir . "blastResult.out";
  my $call1 = system("makeblastdb -in $prot -dbtype prot") unless (-f $prot . ".phr");
  my $call2 = system("blastx -db $prot -query $exome -outfmt 6 -out $out -evalue $eval2") unless (-f $out);

  my %blast;
  open(IN, "<$out");
  while(<IN>) {
    chomp(my $line = $_);
    my @d = split(/\t/,$line);
    $blast{$d[0]} = $d[1] unless $blast{$d[0]};
  }
  close(IN);	
  unlink ($out);
  return(\%blast);
}

## refine the bed file, remove redudancies    
sub modbed {
  my ($bed, $resdir, $name, $offset, $seq) = @_;
  
  my %newbed = %{$bed};
  my %seq = %{$seq};
  my %exonflankingbed;
  my %flankingbed;
  
  my $rebed = $resdir . $name . "_bed.txt";
  open (BED, ">", $rebed); 
    
  foreach my $contig (sort {$a cmp $b} keys %newbed) { 
    for (my $i = 0; $i < scalar (@{$newbed{$contig}}); $i++) {
      print BED $contig, "\t", $newbed{$contig}[$i][0], "\t", $newbed{$contig}[$i][1], "\n"      
    }
  }
  close BED;
  MakeBed ($rebed, 0);

  my %bed;
  
  open (BED, "<", $rebed);
  while (<BED>) {
    chomp (my @l = split /\s+/,$_);
    push @{$bed{$l[0]}},  [ $l[1],$l[2]];

  }
  close BED;
  unlink ($rebed);
  
  foreach my $contig (sort {$a cmp $b} keys %bed) { 
    my $tmpallbed = $resdir . $name . "_" . $contig . "_all_tmp.bed";
    my $tmpflankingbed = $resdir . $name . "_" . $contig . "_flanking_tmp.bed";
    
    open (BEDALL, ">",$tmpallbed);
    open (BEDF, ">", $tmpflankingbed);
    
    my $contiglength = length ($seq{$contig});
    
    for (my $i = 0; $i < scalar (@{$bed{$contig}}); $i++) {
      
      my $start;
      my $end;
      my $nextS;
      my $formerE;
      
      $start = $bed{$contig}[$i][0];
      $end = $bed{$contig}[$i][1];
           
     
      $nextS = $bed{$contig}[$i+1][0] if ($i < scalar (@{$bed{$contig}}) -1);
      $nextS = $contiglength if ($i == scalar (@{$bed{$contig}}) -1);
      $formerE = 0 if ($i == 0);
      $formerE = $bed{$contig}[$i-1][1] if  ($i > 0);
      
      my $news = $start - $offset;
      $news = $formerE if $news < $formerE;
      
      my $newe = $end + $offset;
      $newe = $nextS if $newe > $nextS;	  
      print "BED:", $contig, "\t", $start, "\t", $end, "\n" if ($newe < $news);
      #print "ALLBED:",  $contig, "\t", $news, "\t", $newe, "\n" if ($newe < $news);
      print $contig, "\t", " is inverted bed file and removed!", "\n" if ($newe < $news);
      #last if ($newe < $news);
      delete $seq{$contig} if ($newe < $news);
      
      print BEDALL $contig, "\t", $news, "\t", $newe, "\n" if ($bed{$contig}[$i][1] != 0 && $newe > $news );
      
      
      #exit if ($newe < $news);
      unless ($bed{$contig}[$i][1] == 0) {
	my $newf1s = 0;
	my $newf1e = 0;
	my $newf2s = 0;
	my $newf2e = 0;
      
	if ($news < $start) {
	  $newf1s = $news;
	  $newf1e = $start;	    
	}
	if ($newe > $end) {
	  $newf2s = $end;
	  $newf2e = $newe;	    
	}
	print BEDF  $contig, "\t", $newf1s, "\t", $newf1e, "\n" if (($newf1e - $newf1s) > 20 && $newe > $news) ;
	#print $contig, "\t", $newf1s, "\t", $newf1e, "\n" if (($newf1e - $newf1s) > 20) ;
	print BEDF  $contig, "\t", $newf2s, "\t", $newf2e , "\n" if (($newf2e - $newf2s) > 20 && $newe > $news);
	#print  $contig, "\t", $newf2s, "\t", $newf2e , "\n" if (($newf2e - $newf2s) > 20);
      }
      
    } ## or (my $i = 0; $i < scalar (@{$bed{$contig}}); $i++) {	
    close BEDALL;
    close BEDF;
    

    #my $tmpallbed1 = $tmpallbed . ".1";
    #system ("cp $tmpallbed $tmpallbed1");
    
    MakeBed ($tmpallbed, 0) unless -z $tmpallbed;
    MakeBed ($tmpflankingbed, 0) unless -z $tmpflankingbed;
    
    unless (-z $tmpallbed ) {
      open (NEWBEDALL, "<", $tmpallbed);
      while (<NEWBEDALL>) {
	chomp (my @l = split /\s+/, $_);
	push @{$exonflankingbed{$l[0]}}, [$l[1] , $l[2]] if ($l[2] - $l[1] > 20 );
      }
      close NEWBEDALL;
    }
    unless (-z $tmpflankingbed ) {
      open (NEWBEDFL, "<", $tmpflankingbed);
      while (<NEWBEDFL>) {
	chomp (my @l = split /\s+/, $_);
	push @{$flankingbed{$l[0]}}, [$l[1], $l[2]] if ($l[2] - $l[1] > 20 );
      }
      close NEWBEDFL;
    }
    unlink ($tmpallbed, $tmpflankingbed);
    
  } #foreach my $contig (sort {$a cmp $b} keys %bed)
  
  return (\%exonflankingbed, \%flankingbed, \%seq,\%bed);
  
  
}

sub removeN { 
  my ($line) = @_;
  if ($line =~ /^[n|N]*[A|C|G|T|Y|S|K|R|W|M]+/ || $line =~ /[A|C|G|T|Y|S|K|R|W|M]+[n|N]*$/) {
    $line =~ s/^[n|N]*//g if $line =~ /^[n|N]*/;
    $line = reverse $line;
    $line =~ s/^[n|N]*//g if $line =~ /^[n|N]*/;
    $line = reverse $line;
  }
  return $line; 
}  

sub match {
  my ($item, $array) = @_;
  my $answer = 'no';
  
  my @array = @{$array};
  foreach (@array) {
    my $line = $_;
    $answer = 'yes' if $line eq $item;
  }
  return ($answer);
}

###parsing the blast to find best matches
sub removeOverlap {
  my ($array, $contigs, $a, $maxOverlap) = @_;
  my %contigs = %{$contigs};
  my %a = %{$a};
  my @overlap;
  my $c = 0;
  
  for (my $i = 0; $i < scalar(@$array); $i++) {
    #length of the assembled contig
    my $leni = length ($contigs {$array->[$i]->[0]}) ;
    
    my %bp;
    for (my $n = min($array->[$i]->[8],$array->[$i]->[9]); $n <= max($array->[$i]->[8],$array->[$i]->[9]); $n++) { 
      $bp{$n}++;
    }
    for (my $j = $i+1; $j < scalar(@$array); $j++) { 
      my $lenj = length ($contigs {$array->[$j]->[0]}) ;
      my $overlap = 0;
      
      for (my $m = min($array->[$j]->[8],$array->[$j]->[9]); $m <= max($array->[$j]->[8],$array->[$j]->[9]); $m++) {
	$overlap++ if $bp{$m};
      }         
      $overlap = $overlap / min(abs($array->[$i]->[8] - $array->[$i]->[9]), abs($array->[$j]->[8] - $array->[$j]->[9]));
      
      if ($overlap > $maxOverlap) {
	if  ($array->[$i] && $array->[$j]) {	 
	  if ( ($array->[$i]->[0]) eq ($array->[$j]->[0]) ) {
	    #print $array->[$i]->[0] , "\t", $array->[$j]->[0], "\n";
	    $overlap[$c]-> [0] = abs($array->[$i]->[9] - $array->[$i]->[8]);
	    $overlap[$c]-> [1] = $array->[$i]->[6];
	    $overlap[$c]-> [2] = $array->[$i]->[7];
	    $overlap[$c]-> [3] = $array->[$i]->[0];
	    $c++;
	    $overlap[$c]-> [0] = abs($array->[$j]->[9] - $array->[$j]->[8]);
	    $overlap[$c]-> [1] = $array->[$j]->[6];
	    $overlap[$c]-> [2] = $array->[$j]->[7];
	    $overlap[$c]-> [3] = $array->[$j]->[0];
	    $c++;
	    
	    if ( abs($array->[$i]->[8] - $array->[$i]->[9]) > abs($array->[$j]->[8] - $array->[$j]->[9]) ) {
	      splice(@$array,$j,1);	      
	    }
	    elsif ( abs($array->[$i]->[8] - $array->[$i]->[9]) <= abs($array->[$j]->[8] - $array->[$j]->[9]) ) {
	      splice(@$array,$i,1);	     
	    }
	    else {
	      next;
	    }
	  }#if ( ($array->[$i]->[0]) eq ($array->[$j]->[0]) ) {


	  ### here I use 20bp as a cutoff, too stringent? 
	  else  {
	    if ( abs($array->[$i]->[8] - $array->[$i]->[9]) - abs($array->[$j]->[8] - $array->[$j]->[9]) > 20 ) { 
	      splice(@$array,$j,1);				
	    }
	    elsif ( abs($array->[$j]->[8] - $array->[$j]->[9]) - abs($array->[$i]->[8] - $array->[$i]->[9]) > 20 ) {
	      splice(@$array,$i,1);			
	    }
	    
	    elsif ( abs (abs($array->[$j]->[8] - $array->[$j]->[9]) - abs($array->[$i]->[8] - $array->[$i]->[9])) <= 20 ) {

	      if ( $lenj > $leni ) {
		my $name = $array->[$i]->[0];
		for (my $m = 0; $m < scalar(@$array); $m ++) {
		  splice(@$array,$m,1) if $array->[$m]->[0] eq $name;		  
		}
	      }
	      
	      elsif ( $lenj <= $leni ) {
		my $name = $array->[$j]->[0];
		for (my $m = 0; $m < scalar(@$array); $m ++) {		  
		  splice(@$array,$m,1) if $array->[$m]->[0] eq $name;		  
		}
	      }	
	    }
	  } #else
	} #if  ($array->[$i] && $array->[$j])
	removeOverlap($array, \%contigs, \%a,$maxOverlap);
	last;
      }#if ($overlap > $maxOverlap) {
    }#for (my $j = $i+1; $j < scalar(@$array); $j++) { 
  }
  return($array, \%a, \@overlap);	
}



######################################################
sub refineorf {
  my ($start, $end, $contig, $dir, $lib) = @_;
  my $newstart = $start;
  my $newend = $end;
  my $embossin = $dir . $lib. ".subin";
  my $embossout = $dir . $lib. ".subout";
  
  my $subseq = substr ($contig, $start-1, $end-$start+1);

  open (SUBIN, ">", $embossin);
  print SUBIN ">subseq\n";
  print SUBIN  $subseq, "\n";
  close SUBIN;
  
  system ("getorf -sequence  $embossin  -outseq $embossout -minsize 10 -maxsize 10000000 -find 2 -reverse N -auto -warning N -error N  -fatal N -die N");
  my %emboss;
  unless (-z $embossout) { 
    open (EMBOSS, "<", $embossout); 
    while (<EMBOSS>) {
      chomp (my $line = $_);
      if ($line =~ m/^>(\S+)_(\d+)\s\[(\d+)\s\-\s(\d+)\]/) {
	my $number = $2;
	my $length = $4 - $3 + 1;
	$emboss{$number}{'s'} = $3;
	$emboss{$number}{'e'} = $4;
	$emboss{$number}{'len'} = $length;
      }
    }
    close EMBOSS;
    unlink ($embossin, $embossout);    
  }
  else {
    print "emboss output is empty?!\n";
    #exit;
  }
  foreach my $number (sort {$emboss{$b}{'len'} <=> $emboss{$a}{'len'}} keys %emboss) {
    #print $emboss{$number}{'s'}, "\t", $emboss{$number}{'e'}, "\n";
    $newstart = $start + $emboss{$number}{'s'} -1;     
    $newend = $start + $emboss{$number}{'e'} - 1;
    last;
  }
  return ($newstart, $newend);
  
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

