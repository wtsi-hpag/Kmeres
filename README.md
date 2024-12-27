# Kmeres v1.0
Kmeres - Identification of Centromeres and Telomeres using Kmer Profiles  
 
### Download and Compile:

    $ git clone  https://github.com/wtsi-hpag/Kmeres.git 
    $ cd kmeres 
    $ ./make
		
### Run the pipelines

#### Run kmeres:
           $ /nfs/users/nfs_z/zn1/src/kmeres/kmeres hash -kmer 31 -sample Human CHM13v2.fa meres.fa  > cover.out \
           
	       Parameters:
             kmer:         Size of kmer word  [ default = 31 ]
             sample:       Sample name 
             CHM13v2.fa:   Input genome assembly file 
             meres.fa:     Output centromere sequence file 

    Please contact Zemin Ning ( zn1@sanger.ac.uk ) for any further information. 
