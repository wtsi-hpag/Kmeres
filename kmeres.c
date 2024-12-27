/***********************************************************************\
 *                                                                     * 
 *                        PROJECT        Kmere                         *
 *                                                                     * 
 *---------------------------------------------------------------------*
 *                                                                     * 
 *        Kmere - Identification of Centromere and Telomere            * 
 *                       using Kmer Profiles                           *
 *                                                                     *
 *                                By                                   *
 *                                                                     *
 *                            Zemin Ning                               *
 *                                                                     *
 *          Copyright (C) 2024-2025 by Genome Research Limited         *
 *                         All rights reserved                         *
 *                                                                     *
 *---------------------------------------------------------------------*
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
 *---------------------------------------------------------------------*/


/****************************************************************************/

 
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

static int *list,*readlength;
static B64_long *hist;
static signed char *head;
static B64_long n_Entry;
static char sample[30]; 
static char Cov_Name[50]; 
static char Pix_Name[50]; 
static char Out_Name[50]; 
static char Plot_Name[50]; 
static char *SCGname=NULL;//SingleCopyGenome
static B64_long baseSize;
static double set_rate = 10.0;
static int n_step = 0;
static int n_pixel = 32768;
static int n_pixlen = 94000;
static int n_denoise = 0;
static int num_reads=40000000;

static void hpsortl(B64_long n, unsigned B64_long *ra)
{
	B64_long l,j,ir,i;
	unsigned B64_long rra;

	if(n < 2)
          return;
        l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l];
		} else {
			rra=ra[ir];
			ra[ir]=ra[1];
			if (--ir == 1) {
				ra[1]=rra;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir)	{
			if (j < ir && ra[j] < ra[j+1]) ++j;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		ra[i]=rra;
	}
}


/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void fastaSort(char **argv, int argc, int SEG_LEN)
/* =============================================  */
{
     B64_long i,j,k,n_Sbase=SEG_LEN;
     B64_long ps=0,ns=0,IntSeg=0,IntBase=0,lii;
     fasta *seq;
     int nSeq;
     B64_long totalBases;
     B64_long mhistc = 0;
     B64_long mhistcc = 0;
     int nsorts = 1024;
     int nshift = (n_Sbase<<1)-10;//2^10=nsorts
     int qthresh=23;
//     int tb = 0;
     int args,ac;
     int m_st,m_ed,step_len = 2,pmod;
     int nseq = 0,seqc,n_patch;
     int sshift;
     unsigned B64_long nmask;
     int id_read[64],num_sect,n_reads; //id_read will never exceed mhist
     int rd_st,rd_ed,n_input;
     char line[500] = {0},outName[60]={0},outFast[60]={0},syscmd[200]={0};
     char *ptr,base[20],zero[20]={0},line2[500]={0};
     FILE *fpMate,*fpOutname,*fpOutfast;
     int rd,read_pair[200];
     B64_long kmask = (1L<<(n_Sbase<<1))-1;
     double Qerr[100];
     B64_long gtBases=0,nclip=0;
     Qerr[0] = 1.0;
     void Hash_Process(int brr,int len, char **argv, int argc, int args);

     for(i=1;i<100;i++) Qerr[i] = pow((double)10.,(double)-.1*i);

/*   sort all the names of genes or name entries   */
     printf("Input data starts ...\n");
     args=2;
     memset(sample,'\0',30);
     strcpy(sample,"Seq-sample");

     for(i=1;i<argc;i++)
     {
       if(!strcmp(argv[i],"-kmer"))
       {
	 i++;
         args=args+2;
       }
       else if(!strcmp(argv[i],"-sample"))
       {
         memset(sample,'\0',30);
         sscanf(argv[++i],"%s",sample);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-denoise"))
       {
         sscanf(argv[++i],"%d",&n_denoise);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-step"))
       {
         sscanf(argv[++i],"%d",&n_step);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-pixel"))
       {
         sscanf(argv[++i],"%d",&n_pixel);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-rate"))
       {
         sscanf(argv[++i],"%lf",&set_rate);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-SCG"))
       {
	 SCGname = argv[++i];
         args=args+2;
       } 
     }
     for(ac=0;ac<argc;ac++)
	printf("%s ",argv[ac]);
     printf("\n");
     printf("kmer size:         %ld\n",n_Sbase);

     if((hist = (B64_long *)calloc(400,sizeof(B64_long))) == NULL)
     {
        printf("ERROR ssaha_init: calloc - hist\n");
        exit(1);
     }

     pmod=1;

     nseq = 0;
     if(strcmp(argv[1],"hash")==0)
       Hash_Process(nseq,SEG_LEN,argv,argc,args);
     else
     {
       printf("Invalid input command!\n");
       exit(1);
     }
     printf("All jobs finished\n");
     fflush(stdout);
     system("ps aux |date");
     return;
}


int ssaha_init( char **argv, int argc, int SEG_LEN)
{

    fastaSort(argv,argc,SEG_LEN);
    return(1);   
}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Hash_Process(int nRead,int SEG_LEN, char **argv, int argc, int args)
/* =============================================  */
{
     B64_long i,j,k,w,iseq,n_Sbase=SEG_LEN;
     B64_long ps=0,ns=0,IntSeg=0,IntBase=0,lii;
     FILE *namef,*namef2,*namef3;
     char *b;
     fasta *seqp;
     fasta *seq;
     char nameout_hash[100],nameout_idex[100],nameout_loci[100],nameout_size[100],nameout_list[100];
     int nSeq;
     B64_long totalBases,nsize[10],n_pix,n_pbase,max_cover;
     B64_long mhistc = 0;
     B64_long mhistcc = 0;
     int de_noise[202],de_index[202];
     int nsorts = 1024;
     int nshift; 
     int cover_avechr; 
     int n_meres,*chr_index,*chr_idpix,*chr_locus,*chr_drate;
     void Plot_Process(char **argv, int argc, int args);
     void Out_Meres(fasta *seq,int nSeq,int n_meres,int *chr_index,int *chr_idpix,int *chr_locus,char **argv,int argc,int args);
     void ArraySort_Int2(int n, int *arr, int *brr);
     B64_long offset,sizesm,dsize_4;
     int stopflag,num_copy,m,n;
     int n_ungroup,split_flag;
     int ac,max_nhit,*ctg_list,*ctg_cover;
     int seqc,num_refseq;
     unsigned B64_long seqcharc = 0;
     int step_flag,NLINK,max_contig,sshift;
     unsigned B64_long nmask;
     B64_long *patch_s_array,*patch_head,*ctg_head,*ray;
     int *patch_index,*patch_list,*patch_ofset,*patch_cdex,*patch_iddex,*dex; 
     int *n_list,*rd_used,*rd_rept,*rd_name,*rd_link;
     int id_read[64],num_sect,n_reads; //id_read will never exceed mhist
     int *cover_ind,*cover_max,*list_pix,*head_pix;
     int rd_st,rd_ed,n_input,n2_contig;
     int **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     unsigned int **uimatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     B64_long **limatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     void ArraySort_Mix(int n, B64_long *arr, int *brr);
     char syscmd[200]={0};
     B64_long kmask = (1L<<(n_Sbase<<1))-1;
     double Qerr[100],*cover_pix,*cover_chr;
     B64_long gtBases=0,nclip=0;
     Qerr[0] = 1.0;
     FILE *fp,*fpHASH,*fpIDEX,*fpLOCI,*fpLIST,*fpSIZE;
     fasta *segg;
     B64_long Size_q_pdata;
     int num_seqque;
     char *pdata;

     if((fp=fopen(argv[args],"rb"))==NULL) printf("Cannot open file\n");
       fseek(fp, 0, SEEK_END);
     Size_q_pdata = ftell(fp) + 1;
     fclose(fp);
     if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL)
       printf("calloc pdata\n");
     num_seqque = extractFastq(argv[args],pdata,Size_q_pdata);
     if((segg=(fasta*)calloc((num_seqque+1),sizeof(fasta)))==NULL)
       printf("calloc segg\n");
     if((seq=decodeFastq(argv[args],&num_seqque,&totalBases,pdata,Size_q_pdata,segg))==NULL)
       printf("no query data found.\n");
     nRead = num_seqque;
     nSeq = num_seqque;
     fastaUC(seq,nRead);
     printf("Number of Chromosomes: %d \n",nSeq);
     printf("Number of Bases: %ld \n",totalBases);

     if((ctg_head = (B64_long *)calloc(nSeq+100,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_head\n");
       exit(1);
     }
     if((ctg_list = (int *)calloc(nSeq+100,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_list\n");
       exit(1);
     }
     if((ctg_cover = (int *)calloc(totalBases+100,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_list\n");
       exit(1);
     }
     mhistc = 0;
     gtBases = 0;
     nsorts = num_reads;
     for(i=0;nsorts > 10;i++)
     {
        nsorts = nsorts>>1;
     } 
     nsorts = 1L<<i;
      
     if(n_Sbase <= 31)
       nshift = (n_Sbase<<1)-i;
     else
       nshift = (31<<1)-i;
     nshift = 40;
			 
     if((patch_list= (int *)calloc(nsorts,sizeof(int))) == NULL)
     {
       printf("phusion: calloc - patch_list\n");
       exit(1);
     }
     if((patch_head= (B64_long *)calloc(nsorts,sizeof(B64_long))) == NULL)
     {
       printf("phusion: calloc - patch_head\n");
       exit(1);
     }
     num_refseq = 0;
     for(iseq=0;iseq<nSeq;iseq++)
     {
           double e=0.0,eh[32];
 	   B64_long IntBaseRC;
           B64_long IntSegRC=0;
	   int pos=0;
           char *q;

           IntSeg=0;
           seqp=seq+iseq;
           b = (char*)seqp->data;
           if(seqp->finished) {
             q = NULL;
           } else {
             q = NULL; //(char*)seqp->qual;
           }

           ns=seqp->length;
	   ctg_list[iseq] = ns;
           printf("chrlen: %ld %ld %d\n",iseq+1,ns,nshift);
           seqcharc += strlen(seqp->name);
           k=n_Sbase;
	   if(k > ns) k=ns;

           for(;k>0;k--,b++,pos++)
           {
              if     (*b == 'A') IntBase=0;
              else if(*b == 'C') IntBase=1;
              else if(*b == 'G') IntBase=2;
              else if(*b == 'T') IntBase=3;
              else
              {
	        e += eh[pos] = 1.;
	        continue;
              }
              e += eh[pos] = (q != NULL) ? Qerr[q[pos]] : 0;
	      IntBaseRC=(IntBase^3)<<((n_Sbase-k)<<1);
              IntBase=IntBase<<((k-1)<<1);
              IntSeg+=IntBase;
	      IntSegRC+=IntBaseRC;
           }
           for(j=ns-n_Sbase;--j>=0;b++,pos++)
           {
	      int p=pos%n_Sbase;
	      if(e < .99) 
              {
	        if(IntSeg > IntSegRC) 
                {
                  int itt = IntSegRC>>nshift;
                  patch_list[itt]++;
	        } 
                else 
                {
                  int itt = IntSeg>>nshift;
                  patch_list[itt]++;
	        }
	        gtBases++;
	      }
	      e -= eh[p];
	      IntSeg = (IntSeg<<2)&kmask;
	      IntSegRC = (IntSegRC>>2)&kmask;
              if     (*b == 'A') IntBase=0;
              else if(*b == 'C') IntBase=1;
              else if(*b == 'G') IntBase=2;
              else if(*b == 'T') IntBase=3;
              else
              {
	        e += eh[p] = 1.;
	        continue;
              }
	      e += eh[p] = (q != NULL) ? Qerr[q[pos]] : 0; 
	      IntBaseRC=(IntBase^3)<<((n_Sbase-1)<<1);
              IntSeg+=IntBase;
	      IntSegRC+=IntBaseRC;
           }

	   if(e < .99) 
           {
	     if(IntSeg > IntSegRC) 
             {
               int itt = IntSegRC>>nshift;
               patch_list[itt]++;
	     } 
             else 
             {
               int itt = IntSeg>>nshift;
               patch_list[itt]++;
	     }
	     gtBases++;
           }
     }

     printf("data reading finished ...\n");
     fflush(stdout);
     system("ps aux | date");

     ctg_head[0] = 0;
     for(i=1;i<nSeq;i++)
        ctg_head[i] = ctg_head[i-1]+ctg_list[i-1];

     printf("Input data finished one set, nsorts=%d, totalgoodBP=%ld\n",nsorts,gtBases);
     fflush(stdout);
     system("ps aux | grep HASH; date");
     max_nhit=0;
     patch_head[0]=0;
     for(i=0;i<nsorts;i++) 
     {
        if(patch_list[i]>max_nhit)
          max_nhit = patch_list[i];
        if(i>0)
          patch_head[i] = patch_head[i-1] + patch_list[i-1];
//        patch_list[i] = 0;
     }
     memset(patch_list,0,4*nsorts);  
     mhistcc=patch_head[nsorts-1]+patch_list[nsorts-1]+20000;

     printf("setting array memory: %ld\n",mhistcc);
     fflush(stdout);
     system("ps aux | grep HASH; date");
     if((patch_index= (int *)calloc(mhistcc,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - patch_index\n");
       exit(1);
     }
     if((patch_ofset= (int *)calloc(mhistcc,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - patch_index\n");
       exit(1);
     }

     if((patch_s_array= (B64_long *)calloc(mhistcc+50000,sizeof(B64_long))) == NULL)
     {
       printf("ssaha: calloc - patch_s_array\n");
       exit(1);
     }
     if((patch_cdex= (int *)calloc(mhistcc+50000,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - patch_index\n");
       exit(1);
     }
     fflush(stdout);
     system("ps aux | grep HASH; date");

     seqc = 0;
       printf("Here: %ld\n",mhistcc);

     mhistcc=0;

     for(iseq=0;iseq<nSeq;iseq++,seqc++)
     {
        char *q;
        double e=0,eh[32];
        B64_long IntBaseRC;
        B64_long IntSegRC=0;
        B64_long IntSeg=0;
        int ic,pos=0;

        seqp=seq+iseq;
        b = (char*)seqp->data; 
        if(seqp->finished) {
          q = NULL;
        } else {
          q = NULL; //(char*)seqp->qual;
        } 
        ns=(seqp->length);

        for(k=n_Sbase;k>0;k--,b++,pos++)
        {
           if     (*b == 'A') IntBase=0;
           else if(*b == 'C') IntBase=1;
           else if(*b == 'G') IntBase=2;
           else if(*b == 'T') IntBase=3;
           else
           {
	     e += eh[pos] = 1.;
	     continue;
           }
           e += eh[pos] = (q != NULL) ? Qerr[q[pos]] : 0;
	   IntBaseRC=(IntBase^3)<<((n_Sbase-k)<<1);
           IntBase=IntBase<<((k-1)<<1);
           IntSeg+=IntBase;
	   IntSegRC+=IntBaseRC;
        }
        for(j=ns-n_Sbase;--j>=0;b++,pos++)
        {
	   int p=pos%n_Sbase;
	   if(e < .99) 
           {
             B64_long patch_pos;
    	     if(IntSeg > IntSegRC) 
             {
               int itt = IntSegRC>>nshift;
               patch_pos = patch_head[itt]+patch_list[itt];
               patch_index[patch_pos]=seqc;
               patch_ofset[patch_pos]=pos-n_Sbase+1;
               patch_s_array[patch_pos] = IntSegRC;
               patch_list[itt]++;
	       mhistcc++;
	     } 
             else 
             {
               int itt = IntSeg>>nshift;
               patch_pos = patch_head[itt]+patch_list[itt];
               patch_index[patch_pos]=seqc;
               patch_ofset[patch_pos]=pos-n_Sbase+1;
               patch_s_array[patch_pos] = IntSeg;
               patch_list[itt]++;
	       mhistcc++;
	     }
//       printf("Herenn: %ld %ld\n",j,mhistcc);
	     gtBases++;
	   }
	   e -= eh[p];
	   IntSeg = (IntSeg<<2)&kmask;
	   IntSegRC = (IntSegRC>>2)&kmask;
           if     (*b == 'A') IntBase=0;
           else if(*b == 'C') IntBase=1;
           else if(*b == 'G') IntBase=2;
           else if(*b == 'T') IntBase=3;
           else
           {
	     e += eh[p] = 1.;
             continue;
           }
	   e += eh[p] = 0;
	   IntBaseRC=(IntBase^3)<<((n_Sbase-1)<<1);
           IntSeg+=IntBase;
	   IntSegRC+=IntBaseRC;
        }
        if(e < .99) 
        {
          B64_long patch_pos;
	  if(IntSeg > IntSegRC) 
          {
            int itt = IntSegRC>>nshift;
            patch_pos = patch_head[itt]+patch_list[itt];
            patch_index[patch_pos]=seqc;
            patch_ofset[patch_pos]=pos-n_Sbase+1;
            patch_s_array[patch_pos] = IntSegRC;
            patch_list[itt]++;
	    mhistcc++;
	  } 
          else 
          {
            int itt = IntSeg>>nshift;
            patch_pos = patch_head[itt]+patch_list[itt];
            patch_index[patch_pos]=seqc;
            patch_ofset[patch_pos]=pos-n_Sbase+1;
            patch_s_array[patch_pos] = IntSeg;
            patch_list[itt]++;
	    mhistcc++;
	  }
	  gtBases++;
        }
     }

       printf("Here: %ld\n",mhistcc);

     ray = patch_s_array;
     for(i=0;i<nsorts;i++) 
     {
           if(patch_list[i] > 1)
           {
	     B64_long patch_offpos;
             if(i > 0)
               patch_offpos = patch_head[i]-patch_head[0];
             else
               patch_offpos = 0;
             for(k=0;k<patch_list[i];k++)
		patch_cdex[k] = k;  
             ArraySort_Mix(patch_list[i],ray+patch_offpos,patch_cdex);
             for(k=0;k<patch_list[i];k++)
             {
                B64_long patch_pos = patch_head[i]-patch_head[0];
                stopflag = 0;
	     	j=k+1;
	      	while((j<patch_list[i])&&(stopflag==0))
	      	{
		    if(patch_s_array[patch_pos+j] == patch_s_array[patch_pos+k])
		    {
			j++;
	     	    }
		    else
		     	stopflag=1;
	       	}
		num_copy = j-k;
		for(m=0;m<num_copy;m++)
		{
                   B64_long idd = patch_pos+patch_cdex[k+m];
	  	   int ictg = patch_index[idd];
                   B64_long idt = ctg_head[ictg]+patch_ofset[idd];

		   ctg_cover[idt]=num_copy;
		  
//                   if(patch_list[i] >= 9)
//                     printf("kmer: %ld %d %d %ld %d %d || %ld %d %d %ld %ld\n",k,patch_list[i],ictg,patch_s_array[patch_pos+k+m],patch_ofset[idd],patch_index[idd],idt,ctg_cover[idt],num_copy,idt,ctg_head[ictg]);
		}
		k = j-1;
             }
             
           }
	   else if(patch_list[i] == 1)
	   {
                B64_long patch_pos = patch_head[i]-patch_head[0];
                B64_long idd = patch_pos;
	        int ictg = patch_index[idd];
                B64_long idt = ctg_head[ictg]+patch_ofset[idd];
		ctg_cover[idt]=1;
//                     printf("kmer: %ld %d %d %ld %d %d || %ld %d %d %ld %ld\n",i,patch_list[i],patch_s_array[patch_pos+k+m],patch_ofset[idd],patch_index[idd],idt,ctg_cover[idt],num_copy,idt,ctg_head[ictg]);
	   }
     }

     memset(Plot_Name,'\0',50);
     sprintf(Plot_Name,"plot-%s.sh",sample);
     memset(Out_Name,'\0',50);
     sprintf(Out_Name,"chr2meres-%s.dat",sample);
     memset(Pix_Name,'\0',50);
     sprintf(Pix_Name,"chr2pix-%s.dat",sample);
     memset(Cov_Name,'\0',50);
     sprintf(Cov_Name,"cover-pix-%s.dat",sample);

     if((namef = fopen(Pix_Name,"w")) == NULL)
     {
       printf("ERROR main:: args \n");
       exit(1);
     }
     if((cover_chr= (double *)calloc(nSeq,sizeof(double))) == NULL)
     {
       printf("phusion: calloc - patch_head\n");
       exit(1);
     }
     if((cover_pix= (double *)calloc(n_pixel,sizeof(double))) == NULL)
     {
       printf("phusion: calloc - patch_head\n");
       exit(1);
     }
     if((cover_ind = (int *)calloc(n_pixel,sizeof(int))) == NULL)
     {
       printf("phusion: calloc - chr_index\n");
       exit(1);
     }
     if((cover_max = (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("phusion: calloc - chr_index\n");
       exit(1);
     }
     if((list_pix = (int *)calloc(n_pixel,sizeof(int))) == NULL)
     {
       printf("phusion: calloc - chr_index\n");
       exit(1);
     }
     if((head_pix = (int *)calloc(n_pixel,sizeof(int))) == NULL)
     {
       printf("phusion: calloc - chr_index\n");
       exit(1);
     }
     if((chr_index = (int *)calloc(n_pixel,sizeof(int))) == NULL)
     {
       printf("phusion: calloc - chr_index\n");
       exit(1);
     }
     if((chr_idpix = (int *)calloc(n_pixel,sizeof(int))) == NULL)
     {
       printf("phusion: calloc - chr_index\n");
       exit(1);
     }
     if((chr_locus = (int *)calloc(n_pixel,sizeof(int))) == NULL)
     {
       printf("phusion: calloc - chr_index\n");
       exit(1);
     }
     if((chr_drate = (int *)calloc(n_pixel,sizeof(int))) == NULL)
     {
       printf("phusion: calloc - chr_index\n");
       exit(1);
     }

     for(j=0;j<nSeq;j++)
     { 
	double ave_cover; 
	int idd = ctg_list[j];
	long int sum_cover; 
	 
	sum_cover = 0;
        for(i=0;i<idd;i++)
        {
	      sum_cover = sum_cover+ctg_cover[ctg_head[j]+i];
	}
        ave_cover = sum_cover;
	cover_chr[j] = ave_cover/idd;
     }

     n_pix = n_pixel-10;
     n_pbase = totalBases;
     n_pbase = n_pbase/n_pix;
     n_pixlen = n_pbase;

     n_pix = 0;
     max_cover = 0;
     printf("Genome2pix: %ld %d\n",n_pix,1);
     fprintf(namef,"%ld %d\n",n_pix,1);

     if(n_denoise == 10)
     {
       for(j=0;j<nSeq;j++)
       {
	  int winsize = 200; 
	  int idd = ctg_list[j];
	  long int sum_locs;
	 
          for(i=0;i<(ctg_list[j]-201);i++)
          {
                double cover_ave = 0.0;
                double rate = 0.0;
	        int n_lows,n_high,n_max;

                memset(de_noise,'\0',808);
                memset(de_index,'\0',808);
                sum_locs = 0;
		n_max = 0;
                for(m=0;m<200;m++)
                {
                   de_noise[m] = ctg_cover[ctg_head[j]+i+m];
		   sum_locs = de_noise[m]+sum_locs;
                   de_index[m] = m;
		   if(de_noise[m] > n_max)
		     n_max = de_noise[m];
                }
                ArraySort_Int2(200,de_noise,de_index);
                cover_ave = sum_locs;
                cover_ave = cover_ave/200;
		rate = n_max;
		rate = n_max/cover_chr[j];

		if((j==0)&&(i>691500)&&(i<691700))
     printf("hit: %ld %ld %d || %ld %lf %lf || %d %lf\n",j,i,n_max,sum_locs,rate,cover_ave,ctg_cover[ctg_head[j]+i],cover_chr[j]);
                if(rate >= 10.0)
		{
		  n_lows = 0;
		  n_high = 0;
                  for(m=0;m<100;m++)
                  {
		     if(de_index[m] > 100)
                       n_lows++;
		     if(de_index[m] < 100)
                       n_high++;
                  }
     printf("www: %ld %ld %d || %ld %lf %lf || %d %d\n",j,i,n_max,sum_locs,rate,cover_ave,n_lows,n_high);
		  if((n_lows > 20)&&(n_high > 20))
		  {
                    for(m=0;m<200;m++)
		       ctg_cover[ctg_head[j]+i+m] = 1;
		    i = i+200;
		  }
	        }
//     printf("yyy: %ld %ld %d\n",j,i,ctg_cover[ctg_head[j]+i]);
          }
       }
     }  

     for(j=0;j<nSeq;j++)
     { 
	long int ave_cover; 
	int idd = ctg_list[j];
	long int sum_cover; 
	
	idd = idd/n_pbase;
        list_pix[j] = idd;	
        cover_max[j] = 0;
        for(i=1;i<idd;i++)
        {
	   int idt = i*n_pbase;
	   int idi = (i-1)*n_pbase;
           sum_cover = 0;
	   for(k=idi;k<idt;k++)
	   {
	      sum_cover = sum_cover+ctg_cover[ctg_head[j]+k];
	   }
	   ave_cover = sum_cover;
	   ave_cover = ave_cover/n_pbase;
	   cover_pix[n_pix+i] = ave_cover;
	   cover_ind[n_pix+i] = j+1;
	   if(ave_cover > cover_max[j])
	     cover_max[j] = ave_cover;
	  if(ave_cover > max_cover)
             max_cover = ave_cover;
           printf("Chr: %ld %ld %d %ld %ld %ld || %lf %d\n",j+1,i,idt,n_pix+i,ave_cover,ctg_head[j]+idi,cover_chr[j],cover_max[j]);
        }
        printf("Genome2pix: %ld %d\n",n_pix+idd,1);
        fprintf(namef,"%ld %d\n",n_pix+idd,1);
	n_pix = n_pix + idd;
     }
     fclose(namef);

     cover_avechr = 0;
     for(j=0;j<nSeq;j++)
	cover_avechr = cover_avechr+cover_max[j];
     cover_avechr = cover_avechr/nSeq;

     for(j=0;j<nSeq;j++)
     {
	if(cover_max[j] < 1000)
	  cover_max[j] = cover_avechr;
     }

     head_pix[0] = 0;
     for(j=1;j<nSeq;j++)
	     head_pix[j] = head_pix[j-1]+list_pix[j-1];

     if(n_denoise > 0)
     {
       for(j=0;j<nSeq;j++)
       {
 	  for(i=0;i<(list_pix[j]-10);i++)
	  {  
	     int n_max = 0;
	     int idt = head_pix[j]+i;
             for(m=0;m<20;m++)
             {  
	        if(cover_pix[idt+m] > 2.0*cover_chr[j])
	          n_max++;
             }
             if(n_max < 12)
	       cover_pix[idt] = cover_chr[j]/10;
	     else if(cover_pix[idt] < cover_chr[j])
	       cover_pix[idt] =  5.0*cover_chr[j];
	  } 
       }
     }

     if((namef2 = fopen(Cov_Name,"w")) == NULL)
     {
       printf("ERROR main:: args \n");
       exit(1);
     }
     if((namef3 = fopen(Out_Name,"w")) == NULL)
     {
       printf("ERROR main:: args \n");
       exit(1);
     }
     n_pix = 0;
     n_meres = 0;
     for(j=0;j<nSeq;j++)
     { 
	int idd = ctg_list[j];
	idd = idd/n_pbase;

        for(i=1;i<idd;i++)
        {
	   int idt = i*n_pbase;
	   int idi = (i-1)*n_pbase;
	   double dimpix,dimpix2;

	   dimpix = cover_pix[n_pix+i];
	   dimpix = dimpix*100;
	   dimpix2 = dimpix;
	   dimpix2 = dimpix2/cover_max[j];
//	   dimpix = dimpix/cover_max[j];
	   dimpix = dimpix/max_cover;
           printf("Chr-Denoise: %ld %ld %d %lf || %lf %f %lf || %d %lf %ld\n",j+1,n_pix+i,idt,cover_pix[n_pix+i],dimpix,set_rate,cover_chr[j],cover_max[j],dimpix2,max_cover);
           fprintf(namef2,"%ld %lf\n",n_pix+i,dimpix);
	   if((dimpix > set_rate)||(dimpix2 > 80))
	   {
             fprintf(namef3,"%ld %ld %d %ld %lf\n",j+1,i,idt,n_pix+i,dimpix);
	     chr_index[n_meres] = j+1;
	     chr_idpix[n_meres] = i;
	     chr_locus[n_meres] = idt;
	     n_meres++;
	   }
        }
	n_pix = n_pix + idd;
     }
     fclose(namef2);
     fclose(namef3);

     if(n_step > 0)
     { 
       for(i=0;i<nSeq;i++)
       {
          for(k=0;k<ctg_list[i];k++)
	  {
	     if(k%n_step==0)
	       printf("Cover: %ld %ld %ld %d %d %ld\n",i+1,k,ctg_head[i]+k,ctg_cover[ctg_head[i]+k],ctg_list[i],totalBases);
	  }
       }
     } 
     free(patch_list);
     free(patch_head);

     free(patch_s_array);
     free(patch_index);

            printf("here1: %d %d\n",nSeq,n_meres);
     Out_Meres(seq,nSeq,n_meres,chr_index,chr_idpix,chr_locus,argv,argc,args);
     Plot_Process(argv,argc,args);

}

/*   Subroutine to output sequences of centromeres and telomeres   */
/* =============================================  */
void Out_Meres(fasta *seq,int nSeq,int n_meres,int *chr_index,int *chr_idpix,int *chr_locus,char **argv, int argc, int args)
/* =============================================  */
{
     int i,j,k,m,n,numMeres,num_hits;
     int stopflag;
     char *st;
     FILE *namef;
     fasta *seqp;


     if((namef = fopen(argv[args+1],"w")) == NULL)
     {
       printf("ERROR main:: args \n");
       exit(1);
     }
     numMeres = 0;
     for(i=0;i<n_meres;i++)
     {
        stopflag=0;
        j=i+1;
        while((j<n_meres)&&(stopflag==0))
        {
          if((chr_index[i] == chr_index[j])&&((chr_idpix[j]-chr_idpix[j-1]) <= 8))
          {
            j++;
          }
          else
            stopflag=1;
        }
        num_hits = j-i;
        printf("Meres: %d %d || %d %d || %d %d || %d %d\n",chr_index[j-1],i,chr_idpix[j-1],chr_idpix[i],numMeres,chr_locus[j-1]-chr_locus[i],chr_locus[i],chr_locus[j-1]);
	numMeres++;
        if((num_hits >= 2)) 
        {
          int i_st,i_ed,nline,seq_len;
          int set_len  = 100;   
	  seqp = seq+chr_index[i]-1;

	  if(num_hits == 1)
	  {
            i_st = chr_locus[i]-n_pixlen;
	    i_ed = chr_locus[i];
	  }
	  else
	  {
            i_st = chr_locus[i];
	    i_ed = chr_locus[j-1];  
	  }
	  if(i_st <= 0)
	    i_st = 1;
	  if(i_ed > seqp->length)
	    i_ed = seqp->length;

	  seq_len = i_ed-i_st;
          nline = seq_len/60;
          st = seqp->data;
          if(seq_len>=set_len)
          {
            fprintf(namef,">%s_%09d_%09d_%09d_%09d\n",seqp->name,seqp->length,i_st,i_ed,seq_len);
            for(m=0;m<nline;m++)
            {
               for(k=0;k<60;k++,st++)
                  fprintf(namef,"%c",*st);
               fprintf(namef,"\n");
            }
            for(k=0;k<(seq_len-(nline*60));k++,st++)
               fprintf(namef,"%c",*st);
            if(seq_len%60!=0)
              fprintf(namef,"\n");
          } 
        }
	i = j-1;
     }	
     fclose(namef);
}

/*   Subroutine to plot coverage profiles to reveal centromeres   */
/* =============================================  */
void Plot_Process(char **argv, int argc, int args)
/* =============================================  */
{
     int i,j,k,m,n,n_chr,n_oando;
     FILE *namef;
     char KKK1[100],KKK2[100],KKK3[100],KKK4[100],KKK5[100],KKK6[100],KKK7[100],KKK8[100];
     char *st,*ed,Sam_name[30],syscmd[2000],Chr_name[30];

     memset(Sam_name,'\0',30);
     strcpy(Sam_name,sample);
     strcpy(KKK1,"set terminal svg");     
     strcpy(KKK2,"set style line 1 lt 1 lw 3 pt 3 linecolor rgb \\\"black\\\"");     
     strcpy(KKK3,"set style line 2 lt 2 lw 5 pt 1 linecolor rgb \\\"red\\\"");     
     strcpy(KKK4,"set xlabel \\\"Pixel / Chromosome coordinates\\\"");     
     strcpy(KKK5,"set ylabel \\\"Dimensionless Kmer Coverage\\\"");
     sprintf(KKK7,"%s%d%s","plot [ 0 to ",n_pixel," ] [ 0 to 100 ]");
     strcpy(KKK6,"with lines ls 1,");
     strcpy(KKK8,"with points ls 2");
     k = 0;

     if((namef = fopen(Plot_Name,"w")) == NULL)
     {
       printf("ERROR main:: args \n");
       exit(1);
     }

     fprintf(namef,"#!/bin/bash\n");
     fprintf(namef,"\n"); 
     fprintf(namef,"function plotcmd\n");
     fprintf(namef,"\{\n");
     fprintf(namef,"printf \"%s\\n\"\n",KKK1);
     fprintf(namef,"printf \"%s\\n\"\n",KKK2);
     fprintf(namef,"printf \"%s\\n\"\n",KKK3);
     fprintf(namef,"printf \"%s\\n\"\n",KKK4);
     fprintf(namef,"printf \"%s\\n\"\n",KKK5);
//     fprintf(namef,"printf \"%s \\\"%s\\\" title \\\"%s %s \\\" %s\"\n",KKK7,"cover-pix.dat",Sam_name,"Kmer Coverage",KKK6);
       fprintf(namef,"printf \"%s \\\"%s\\\" title \\\"%s %s \\\" %s \\\"%s\\\" title \\\"%s %s \\\" %s\"\n",KKK7,Cov_Name,Sam_name,"Kmer Coverage",KKK6,Pix_Name,Sam_name,"Chromosome",KKK8);
     fprintf(namef,"}\n");
     fprintf(namef,"plotcmd | gnuplot > data.svg\n");
     fprintf(namef,"inkscape -z --export-text-to-path --export-type=pdf -o data.pdf data.svg\n");
     fprintf(namef,"gs -r600 -dNOPAUSE -dBATCH -sDEVICE=png256 -sOutputFile=%s.png data.pdf\n",Sam_name);
     fprintf(namef,"\n"); 
     fprintf(namef,"\n"); 
     fclose(namef);

     memset(syscmd,'\0',2000);
     sprintf(syscmd,"bash %s > try.out",Plot_Name);
     if(system(syscmd) == -1)
     {
        printf("System command error:\n");
     }

     exit(1);


}

/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int  **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
unsigned int     **uimatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=((nrh-nrl+1)/100 + 1)*100,ncol=nch-ncl+1;
        unsigned int  **m;
	B64_long nri,nrn=nrow/100;

        /* allocate pointers to rows        */
        if((m=(unsigned int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }

        /* allocate rows and set pointers to them        */
	/* allocate in 100 batches to use freed memory */
	nrl = 0;
	for(nri=0;nri<100;nri++,nrl+=nrn) {
           if((m[nrl]=(unsigned int *)calloc(nrn*ncol,sizeof(int)))==NULL)
           {
              printf("error imatrix: calloc error No. 2 \n");
              return(NULL);
           }

           for(i=1;i<nrn;i++)
              m[i+nrl]=m[i+nrl-1]+ncol;
	}
       /* return pointer to array of pointers to rows   */
        return m;
}

/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
B64_long     **limatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        B64_long  **m;

        /* allocate pointers to rows        */
        if((m=(B64_long **)calloc(nrow,sizeof(B64_long*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(B64_long *)calloc(nrow*ncol,sizeof(B64_long)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C   */   

/* =============================== */
void ArraySort_Mix(int n, B64_long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* creat char matrix with subscript ange cm[nrl...nrh][ncl...nch]  */
unsigned char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        unsigned char **cm;

        /* allocate pointers to rows        */
        if((cm=(unsigned char **)calloc(nrow,sizeof(unsigned char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        cm+=0;
        cm-=nrl;

        /* allocate rows and set pointers to them        */
        if((cm[nrl]=(unsigned char *)calloc(nrow*ncol,sizeof(unsigned char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        cm[nrl]+=0;
        cm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           cm[i]=cm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return cm;
}

/* =============================== */
void ArraySort_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

