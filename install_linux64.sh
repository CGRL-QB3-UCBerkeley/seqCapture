#! /bin/bash
printf "\n"
echo "now installing dependencies for the sequence capture analysis workflow!"
cd dependencies
printf "\n"

if ! (which conda 1> /dev/null 2>&1)
then
   echo "miniconda2 is not installed! please run bash Miniconda2-latest-Linux-x86_64.sh to install!"
   printf "\n"
   exit 1
fi

CHECKPATH=FLASH2/
if [ ! -d "$CHECKPATH" ]
then
   conda create -n seqCapture
else
   echo "it seems that most dependencies have been installed. now checking status of all these packages..."
   printf "\n"
fi

conda install -n seqCapture --file seqCapture_v0.0.1_linux64_install.txt
source activate seqCapture
printf "\n"


CHECKPATH=FLASH2/ 
if [ ! -d "$CHECKPATH" ]; then
   echo "now installing additional packages needed by the pipeline!" 
   tar -xf morePackages.tar.gz
   printf "\n"
fi

###########flash2################
if ! (which  FLASH2/flash2 1> /dev/null 2>&1)
then 
   echo "now installing flash2!"
   cd FLASH2
   make clean
   make
   cd ..
   printf "\n"
fi
 
############rcorrector#############
if  ! (which Rcorrector/rcorrector 1>/dev/null 2>&1 )
then
    echo "now installing Rcorrector!"
    cd Rcorrector
    make
    cd ..
    printf "\n"
fi

###########supe_dedupper############
if ! ( which Super-Deduper/super_deduper 1>/dev/null 2>&1 )
then
    echo "now installing super_deduper!"
    cd Super-Deduper
    make clean
    make
    cd ..
    printf "\n"
fi

##########RepeatMasker configuration########
cd RepeatMasker
PERL=`which perl`
BINDIR=${PERL%/perl}

#RepeatMaskerConfig.pm:    $RMBLAST_DIR   = "media12345";
if ! grep -q $BINDIR RepeatMaskerConfig.pm
then
   echo "now configuring RepeatMasker!"
   perl changename.pl
   printf "\n"
else
   echo "RepeatMasker is installed, looking good!"
   printf "\n"
fi

cd ..
#########check if gatk is registered###########

GATK=${PERL%/bin/perl}'/opt/gatk-3.8/GenomeAnalysisTK.jar'

if [ -e $GATK ]
then
    echo "GATK is registered, looking good!"
else
    gatk-register GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar
fi
printf "\n"

#########add seqCapture path to ~/.bashrc ##########
cd ..
PWD=`pwd .`
if ! grep -q $PWD ~/.bashrc
then
   echo "export \"PATH=$PWD:\$PATH\"" >> ~/.bashrc
else
   echo "seqCapture is already in your PATH, looking good!"
fi #  !BATCH
printf "\n"

###############add env PERL5LIB PATH ##############3
PERLLIBPATH=${PERL%/bin/perl}'/lib/perl5/5.22.0/lib/perl5' 
PERLINSTALLBASE=${PERL%/bin/perl}'/lib/perl5/5.22.0/'
if ! grep -q $PERLLIBPATH ~/.bashrc
then 
    echo "export \"PERL5LIB=$PERLLIBPATH:\$PERL5LIB\"" >>  ~/.bashrc
else
    echo "$PERLLIBPATH is already in your PERL5LIB PATH, looking good!"
fi #  !BATCH
printf "\n"

#################start installing perl modules using cpanm ################
if  perl -MList::MoreUtils -e 'print "$List::MoreUtils\n"' >/dev/null 2>&1
then 
     echo "List::MoreUtils is already installed, looking good!" 
else 
     echo "installing perl module List::MoreUtils..." 1>&2
     cpanm -q List::MoreUtils  -l $PERLINSTALLBASE
fi
printf "\n"
if  perl -MTie::Array::Packed -e 'print "$Tie::Array::Packed\n"' >/dev/null 2>&1
then 
     echo "Tie::Array::Packed is already installed, looking good!" 
else 
     echo "installing perl module Tie::Array::Packed..." 1>&2
     cpanm -q Tie::Array::Packed  -l $PERLINSTALLBASE
fi
printf "\n"
if  perl -MText::Soundex -e 'print "$Text::Soundex\n"' >/dev/null 2>&1
then 
     echo "Text::Soundex is already installed, looking good!" 
else 
     echo "installing perl module Text::Soundex..." 1>&2
     cpanm  -q Text::Soundex  -l $PERLINSTALLBASE
fi
printf "\n"
if  perl -MSort::Packed -e 'print "$Sort::Packed\n"' >/dev/null 2>&1
then 
     echo "Sort::Packed is already installed, looking good!" 
else 
     echo "installing perl module Sort::Packed..." 1>&2
     cpanm  -q Sort::Packed -l $PERLINSTALLBASE
fi
printf "\n"
if  perl -MBio::PopGen::Statistics -e 'print "$Bio::PopGen::Statistics\n"' >/dev/null 2>&1
then 
     echo "Bio::PopGen::Statistics is already installed, looking good!" 
else 
     echo "installing perl module Bio::PopGen::Statistics..." 1>&2
     cpanm  -q Bio::PopGen::Statistics -l $PERLINSTALLBASE
fi
printf "\n"
if  perl -MStatistics::Basic -e 'print "$Statistics::Basic\n"' >/dev/null 2>&1
then 
     echo "Statistics::Basic is already installed, looking good!" 
else 
     echo "installing perl module Statistics::Basic..." 1>&2
     cpanm  -q Statistics::Basic -l $PERLINSTALLBASE
fi
printf "\n"
if  perl -MMath::Cephes -e 'print "$Math::Cephes\n"' >/dev/null 2>&1
then 
     echo "Math::Cephes is already installed, looking good!" 
else 
     echo "installing perl module Math::Cephes..." 1>&2
     cpanm  -q Math::Cephes -l $PERLINSTALLBASE
fi
printf "\n"
if  perl -MStatistics::Distributions -e 'print "$Statistics::Distributions\n"' >/dev/null 2>&1
then 
     echo "Statistics::Distributions is already installed, looking good!" 
else 
     echo "installing perl module Statistics::Distributions..." 1>&2
     cpanm  -q Statistics::Distributions -l $PERLINSTALLBASE
fi
printf "\n"
if  perl -MText::NSP::Measures::2D::Fisher2::twotailed -e 'print "$Text::NSP::Measures::2D::Fisher2::twotailed\n"' >/dev/null 2>&1
then 
     echo "Text::NSP::Measures::2D::Fisher2::twotailed is already installed, looking good!" 
else 
     echo "installing perl module Text::NSP::Measures::2D::Fisher2::twotailed..." 1>&2
     cpanm  -q Text::NSP::Measures::2D::Fisher2::twotailed -l $PERLINSTALLBASE
fi
printf "\n"
if  perl -MSort::Naturally -e 'print "$Sort::Naturally\n"' >/dev/null 2>&1
then 
     echo "Sort::Naturally is already installed, looking good!" 
else 
     echo "installing perl module Sort::Naturally..." 1>&2
     cpanm  -q Sort::Naturally -l $PERLINSTALLBASE
fi
printf "\n"

echo "do not forget to run source ~/.bashrc right after installation or start a new terminal session"
printf "\n"

#############deactivate#######################
source deactivate seqCapture


