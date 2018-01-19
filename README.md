# seqCapture
Bioinformatic workflow for analyzing sequence capture NGS dataset

Detailed documentaion please read here XXXXX

Install:

1. To obtain the lastest version of seqCapture:

        $ git clone https://github.com/CGRL-QB3-UCBerkeley/seqCapture.git

2. Download your own GATK from https://software.broadinstitute.org/gatk/download/. You will need to be registered and logged into the GATK forum to access the download. Place your GenomeAnalysisTK.jar to the diretory "seqCapture/dependencies/morePackages":
    
        $ cp /path/to/GenomeAnalysisTK.jar  seqCapture/dependencies/morePackages

3. Install seqCapture:

        $ cd seqCapture/
    
    for MacOS users:
    
        $ ./install_osx64.sh
    
    for linux users
    
        $ ./install_linux64.sh
    
Uninstall seqCapture (which is sad):

    $ cd seqCapture/
    
    $ ./uninstall.sh

    

