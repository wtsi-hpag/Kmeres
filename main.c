/**************************************************************************
 *
 *-------------------------------------------------------------------------
 * Project Notes : Assembles anchored reads onto genetic data
 *
 *
 *-------------------------------------------------------------------------
 #######################################################################
 # This software has been created by Genome Research Limited (GRL).    # 
 # GRL hereby grants permission to use, copy, modify and distribute    # 
 # this software and its documentation for non-commercial purposes     # 
 # without fee at the user's own risk on the basis set out below.      #
 # GRL neither undertakes nor accepts any duty whether contractual or  # 
 # otherwise in connection with the software, its use or the use of    # 
 # any derivative, and makes no representations or warranties, express #
 # or implied, concerning the software, its suitability, fitness for   #
 # a particular purpose or non-infringement.                           #
 # In no event shall the authors of the software or GRL be responsible # 
 # or liable for any loss or damage whatsoever arising in any way      # 
 # directly or indirectly out of the use of this software or its       # 
 # derivatives, even if advised of the possibility of such damage.     #
 # Our software can be freely distributed under the conditions set out # 
 # above, and must contain this copyright notice.                      #
 #######################################################################
 *
 *  Author : James C. Mullikin
 *
 *	Copyright (C) 1998-2001 by Genome Research Limited, All rights reserved.
 *
 **************************************************************************/


#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"
#include "ssaha.h"

static int KLEN=17;

int main(int argc, char **argv)
{
	int i;
    int SEG_LEN=KLEN;
//	__delayed_free = 0;//so we can free memory right away

	printf("SsahaHist Version 1.02\n");
	if(argc < 2)
	{
		printf("Usage: %s command [options] \n",argv[0]);
		printf("\n");
		printf("Commands:\n");
		printf("          hash     generate hash files from sequencing reads\n");
		printf("          Options: phusion2 hash <-kmer 31>\n");
		printf("\n");
		printf("          sort     sort hash files files\n");
		printf("          Options: phusion2 sort <-kmer 31> reads_0*.fastq (all the read files)\n");
		printf("\n");
		printf("          edge     generate relation matrix file (sparse)\n");
		printf("          Options: phusion2 edge <-kmer 31> <-depth 20> <-edge 50> reads_0*.fastq\n");
		printf("\n");
		printf("          matrix   generate relation matrix file (compact)\n");
		printf("          Options: phusion2 matrix <-kmer 31> <-edge 50> reads_0*.fastq\n");
		printf("\n");
		printf("          clust    cluster reads into groups:\n");
		printf("          Options: phusion2 clust <-kmer 31> <-depth 20> <-match 4> <-set 60000> > clust-d20.out \n");
		printf("\n");
		printf("          remap    remap the reads which were not used initially\n");
		printf("          Options: phusion2 remap <-kmer 31> <-depth 20> <-match 2> <-set 60> clust-d20.dat clust_remap.dat \n");
		printf("\n");
		printf("          merge    merge clusters using remapped reads:\n");
		printf("          Options: phusion2 merge <-kmer 31> <-depth 20> <-match 2> <-set 60> clust_remap.dat clust_merge.dat \n");
		printf("\n");
		printf("          fmate    output read files in patches:\n");
		printf("          Options: phusion2 fmate <-kmer 31> <-edge 100> <-split 1> clust_merge.dat reads_0*.fastq\n");
		printf("\n");
		printf("Note: for all default values, the use of any parameters is optional.\n");
		exit(1);
	}

        for(i=1;i<argc;i++)
        {
           if(!strcmp(argv[i],"-kmer"))
           {
             sscanf(argv[++i],"%d",&SEG_LEN);
           }
        } 
	ssaha_init(argv, argc, SEG_LEN);
        return EXIT_SUCCESS;
}
