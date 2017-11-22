#! /bin/bash

echo "uninstall seqCapture!"
cd dependencies
rm -r FLASH2 Rcorrector	Super-Deduper GenomeAnalysisTK-3.8-0 RepeatMasker	
cd ..
source activate seqCapture
DIR=`pwd .`
PERL=`which perl`
PERLLIBPATH=${PERL%/bin/perl}'/lib/perl5/5.22.0/lib/perl5'

if [ -e ~/.bashrc ] && ( grep -q $DIR ~/.bashrc )
then
   grep -v $DIR ~/.bashrc | grep -v $PERLLIBPATH  >  ~/.bashrc1
   mv  ~/.bashrc1 ~/.bashrc 
elif [ -e ~/.bash_profile ] && ( grep -q "$DIR" ~/.bash_profile )
then
   grep -v $DIR  ~/.bash_profile | grep -v $PERLLIBPATH  > ~/.bash_profile1
   mv ~/.bash_profile1 ~/.bash_profile   
fi

source deactivate

conda-env remove -n seqCapture

if [ -e ~/.bashrc ]
then 
   echo "Uninstallation finished; do not forget to run source ~/.bashrc"
else
   echo "Uninstallation finished; do not forget to run source ~/.bash_profile" 
fi
printf "\n"
