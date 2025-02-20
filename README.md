# Kmeres v1.0
Kmeres - Identification of Centromeres and Telomeres using Kmer Profiles  
 
### Download and Compile:

    $ git clone  https://github.com/wtsi-hpag/Kmeres.git 
    $ cd Kmeres 
    $ make
		
### Run the pipelines

#### Run kmeres:
           $ /nfs/users/nfs_z/zn1/src/kmeres/kmeres hash -kmer 31 -sample Human CHM13v2.fa meres.fa  > cover.out \
           
	       Parameters:
             kmer:         Size of kmer word  [ default = 31 ]
             sample:       Sample name 
             CHM13v2.fa:   Input genome assembly file 
             meres.fa:     Output centromere sequence file 

### Test dataset:

A test dataset on the T2T gemome assembly of Arabidopsis Thaliana has been put at 

https://ftp.sanger.ac.uk/pub/users/zn1/assembly/Thaliana/

The detected centromeres are in meres.fa

egrep ">" meres.fa

>Chr1_032540122_014841792_017131968_002290176
>Chr2_022217084_003826368_006048000_002221632
>Chr3_025743512_013595904_015736896_002140992
>Chr4_021578073_004201344_006979392_002778048
>Chr5_029480885_011829888_014220864_002390976

Chr1:                                    \
chromosome length 32540122               \
centromere starts from 14841792          \
centromere ends at 17131968              \
centromere length 2290176                \ 


The centromeres predected by Kmeres are very close to those detected by Chip-seq, see

Naish M. et al. 2021 "The genetic and epigenetic landscape of the Arabidopsis centromeres". Science. 2021;374:eabi7489. doi: 10.1126/science.abi7489. \ 

	The users may download the assembly and try it first with a local computer.   \ 
	Then conpare with the output files at the FTP site.     \
	The coverage profile see: Thaliana.png       \

### FTP site for test dataset:

	ftp://ftp.sanger.ac.uk/pub/users/zn1/assembly/Thaliana

 
    Please contact Zemin Ning ( zn1@sanger.ac.uk ) for any further information. 
