/***********************************************************************\
 *                                                                     * 
 *                      PROJECT        SSAHA                            *
 *                                                                     * 
 *---------------------------------------------------------------------*
 *                                                                     * 
 *          Sequence Search and Alignment by Hashing Algorithm         *
 *                                                                     *
 *                                By                                   *
 *                                                                     *
 *                    Zemin Ning & James C. Mullikin                   *
 *                                                                     *
 *          Copyright (C) 1999-2000 by Genome Research Limited         *
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


static int *list;
static B64_long *hist;
static signed char *head;
static B64_long n_Entry;
static int K_LEN=17;
static int mhist=10;
static int outmod=0;
static int NEDGE=500;
static int nmatch=30;
static int N_SET=6000;


static void hpsortil(B64_long n, int *ra, B64_long *rb)
{
	B64_long l,j,ir,i;
	B64_long rra;
	B64_long rrb;

	if(n < 2)
        return;
        l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l];
			rrb=rb[l];
		} else {
			rra=ra[ir];
			rrb=rb[ir];
			ra[ir]=ra[1];
			rb[ir]=rb[1];
			if (--ir == 1) {
				ra[1]=rra;
				rb[1]=rrb;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir)	{
			if (j < ir && ra[j] < ra[j+1]) ++j;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				rb[i]=rb[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		ra[i]=rra;
		rb[i]=rrb;
	}
}

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
void fastaSort(char **argv, int argc, int i_Fasta, int SEG_LEN)
/* =============================================  */
{
     B64_long i,j,k,n_Sbase=SEG_LEN;
     B64_long ps=0,ns=0,IntSeg=0,IntBase=0,lii;
     char *b,*ch,*b0;
     fasta *seqp,*expp[500];
     fasta *seq,*seqr;
     int nSeq;
     B64_long totalBases;
     B64_long mhistc = 0;
     B64_long mhistcc = 0;
     unsigned B64_long *sarr, *sarrp;
     int qthresh=23;
     B64_long CG = 6;
     int CGneg=0;
//     int tb = 0;
     int args,ac;
     int nseq = 0,seqc;
     unsigned B64_long seqcharc = 0;
     char **names,**mates;
     int sshift,rshift=6;
     unsigned B64_long nmask;
     int **r_matrix,*n_list,*rd_used,*rd_rept,*rd_name;
     int id_read[64],num_sect,n_contig,n_reads,max_edge; //id_read will never exceed mhist
     int *num_front,*ind_front,*out_list,*ctg_list;
     int *num_pair,*mat_pair,*ins_pair,*dev_pair;
     int n_patch,n_outreads,new_patch,n_ctg;
     int n_grouped,rd_st,rd_ed,n_input;
     int **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     char line[500] = {0},outName[60]={0},outFast[60]={0},syscmd[200]={0};
     char *ptr,base[20],zero[20]={0},line2[500]={0};
     FILE *fpMate,*fpOutname,*fpOutfast;
     int rd,read_pair[200];

/*   sort all the names of genes or name entries   */
     printf("Input data starts ...\n");
     fflush(stdout);
     system("date");
     n_input=0;
     args=1;
     for(i=1;i<argc;i++)
     {
       if(!strcmp(argv[i],"-kmer"))
       {
         sscanf(argv[++i],"%d",&K_LEN);
         args=args+2;
         printf("kmer size: %d\n",K_LEN);
       } 
       else if(!strcmp(argv[i],"-depth"))
       {
         sscanf(argv[++i],"%d",&mhist);
         args=args+2;
         printf("cutoff number: %d\n",mhist);
       } 
       else if(!strcmp(argv[i],"-match"))
       {
         sscanf(argv[++i],"%d",&nmatch);
         args=args+2;
         printf("match number: %d\n",nmatch);
       } 
       else if(!strcmp(argv[i],"-output"))
       {
         sscanf(argv[++i],"%d",&outmod);
         args=args+2;
         printf("output mode: %d\n",outmod);
       } 
       else if(!strcmp(argv[i],"-matrix"))
       {
         sscanf(argv[++i],"%d",&NEDGE);
         args=args+2;
         printf("matrix size: %d\n",NEDGE);
       } 
       else if(!strcmp(argv[i],"-set"))
       {
         sscanf(argv[++i],"%d",&N_SET);
         args=args+2;
         printf("outreads size: %d\n",N_SET);
       } 
     }
     max_edge=NEDGE;
     n_Sbase=K_LEN; 
     sshift = (64-(n_Sbase<<1));
     nmask = (1<<sshift)-1;

     for(ac=(args+1);ac<argc;ac++)
     {
	seq = decodeFastq ( argv[ac], &nSeq, &totalBases, qthresh);
	if(seq == NULL)
	{
		printf("ERROR pileup: no data found\n");
		exit(1);
	}
	fastaUC(seq,nSeq);
	nseq += nSeq;
        seqp = seq;
        b0= (char*)seq->name;
        for(i=0;i<nSeq;i++,seqp++)
        {
           b = (char*)seqp->data;
           ns=(seqp->length)-n_Sbase;
   	   seqcharc += strlen(seqp->name);
           for(j=0;j<ns;j++)
           {
              B64_long IntSegRC=0;
              IntSeg=0;
              ch=b;
              for(k=n_Sbase;k>0;k--)
              {
	         B64_long IntBaseRC;
                 if     (*ch == 'A') IntBase=0;
                 else if(*ch == 'C') IntBase=1;
                 else if(*ch == 'G') IntBase=2;
                 else if(*ch == 'T') IntBase=3;
                 else
                 {
  		   break;
                 }
	         IntBaseRC=(IntBase^3)<<((n_Sbase-k)<<1);
                 IntBase=IntBase<<((k-1)<<1);
                 IntSeg=IntSeg+IntBase;
		 IntSegRC=IntSegRC+IntBaseRC;
                 ch++;
              }
	      if(0) 
              {
	        B64_long mask=3;
		B64_long jj;
	        char *base = "ACGT";
	        for(jj=n_Sbase;--jj>=0;)
	   	   printf("%c",base[(IntSeg>>(jj<<1))&mask]);
		printf(" ");
		for(jj=n_Sbase;--jj>=0;)
		   printf("%c",base[(IntSegRC>>(jj<<1))&mask]);
	        printf("\n");
              }
              b++;
	      if(IntSeg > IntSegRC) IntSeg = IntSegRC;
	      if(k != 0) continue;
	      if(CGneg)
              {
	        B64_long is=IntSeg;
	        for(k=SEG_LEN-1;--k>=0;)
	        {
	       	   if((is&15)==CG) break;
	 	   is >>= 2;
	        }
		if(k!=-1) 
                {
 		  if(head[IntSeg] != -127) head[IntSeg]--;
		  continue;
		}
 	      }
              if(head[IntSeg] != 127) head[IntSeg]++;
           }
        }
//	 printf("ssaha: %ld %ld\n",seq->name,seq);
        free(seq->name);
	if(seq) free(seq);
        n_input++;
     }
     printf("Input data finished one set\n");
     fflush(stdout);
     system("date");

     for(i=n_Entry;--i>=0;)
     {
        B64_long mask=3;
        B64_long is=i;
	char *base = "ACGT";
//      for(j=n_Sbase;--j>=0;)
//      printf("%c",base[(i>>(j<<1))&mask]);
//      printf(" %8d\n",head[i]);
        if(CGneg)
        {
          for(k=SEG_LEN-1;--k>=0;)
          {
   	     if((is&15)==CG) break;
	       is >>= 2;
	  }
	  if(k!=-1) 
          {
	    head[i]--;
	  }
	}
        if(head[i] < 128 && head[i] >= -128)
	  hist[head[i]+128]++;
        if(head[i] < 0) head[i] = -(head[i]+1);
        if(head[i] > mhist || head[i] <= 1 )continue;
        mhistc += head[i];
     }
     for(i=0;i<256;i++)
	printf("%5d %12ld\n",i-128,hist[i]);

     fflush(stdout);
     if(0) return;
     sarrp = sarr = (unsigned B64_long *) malloc(mhistc*sizeof(unsigned B64_long));
     names = (char **) malloc(nseq*sizeof(char *));
     names[0] = (char *) calloc(seqcharc+nseq,1);
//     tb = 0;
     seqc = 0;
     for(ac=(args+1);ac<argc;ac++)
     {
	seq = decodeFastq ( argv[ac], &nSeq, &totalBases, qthresh);
	if(seq == NULL)
	{
		printf("ERROR pileup: no data found\n");
		exit(1);
	}
	fastaUC(seq,nSeq);
        seqp = seq;
        for(i=0;i<nSeq;i++,seqp++,seqc++)
        {
//	   B64_long apr=0;
	   strcpy(names[seqc],(char*)seqp->name);
	   names[seqc+1] = names[seqc]+strlen((char*)seqp->name)+1;
           b = (char*)seqp->data;
           ns=(seqp->length)-n_Sbase;
           for(j=0;j<ns;j++)
           {
              B64_long IntSegRC=0;
              IntSeg=0;
              ch=b;
              for(k=n_Sbase;k>0;k--)
              {
	         B64_long IntBaseRC;
                 if     (*ch == 'A') IntBase=0;
                 else if(*ch == 'C') IntBase=1;
                 else if(*ch == 'G') IntBase=2;
                 else if(*ch == 'T') IntBase=3;
                 else
                 {
	           break;
                 }
		 IntBaseRC=(IntBase^3)<<((n_Sbase-k)<<1);
                 IntBase=IntBase<<((k-1)<<1);
                 IntSeg=IntSeg+IntBase;
	         IntSegRC=IntSegRC+IntBaseRC;
                 ch++;
              }
              b++;
	      if(IntSeg > IntSegRC) IntSeg = IntSegRC;
	      if(k != 0) continue;
	      if(head[IntSeg] > mhist || head[IntSeg] <= 1 ) continue;
//	      apr++;
//	      tb++;
	      mhistcc++;
	      *sarrp++ = (IntSeg<<sshift) + seqc;
           }
//         printf(" %-30s %4d\n",names[seqc],apr);
        }
//	 printf("ssaha: %ld %ld\n",seq->name,seq);
        free(seq->name);
        if(seq) free(seq);
//      freeRead(&seq);
     }

     free(head);

     printf("Sorting starts ...\n");
     fflush(stdout);
     system("date");
     hpsortl(mhistcc,sarr-1);
     printf("Sorting finished\n");
     printf("Start read groups %ld %ld",mhistc,mhistcc);
     fflush(stdout);
     system("date");
     lii = -1;
/*   allocate memery for arrays to calculate the link number   */

     if((n_list= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - n_list\n");
       exit(1);
     }
     r_matrix=imatrix(0,nseq,0,max_edge);

     num_sect=0;
     for(i=0;i<mhistcc;i++)
     {
        int j = (int)(sarr[i]&nmask);
        B64_long ii=(B64_long)sarr[i]>>sshift;
        int match,idk,idi,kk,ki,km,a_index,a_links;

        if(ii != lii)
        {
          for(kk=0;kk<num_sect;kk++)
          {
             for(ki=0;ki<num_sect;ki++)
             {
                if(ki!=kk)
                {
                  idk=id_read[kk];
                  idi=id_read[ki];
                  if(n_list[idk]<max_edge)
                  {
                    match=0;
                    for(km=0;km<n_list[idk];km++)
                    {
                       a_index=r_matrix[idk][km]>>rshift;
                       a_links=r_matrix[idk][km]&077;
//                       a_links=r_matrix[idk][km]-(a_index<<rshift);
                       if(a_index==idi)
                       {
                         if(a_links<60)
                           r_matrix[idk][km]++;
                         match=1;
                         break;
                       }
                    }
                    if(match==0)
                    {
                      r_matrix[idk][n_list[idk]]=(idi<<rshift)+1;
                      n_list[idk]++;
                    }
                  }
                }
             }
          }
          num_sect=0;
//            putchar ('\n');
        }
        else
        {
//            putchar(' ');
        }
        lii = ii;
//          printf("%s",names[j]);
        id_read[num_sect]=j;
        num_sect++;
     }
     putchar('\n');
//     free(sarrp);dont free this this is a temporary pointer to index sarr
     free(sarr);
     if((ind_front= (int *)calloc(nseq,sizeof(int))) == NULL)//worst case all reads are linked together
     {
       printf("ssaha: calloc - ind_front\n");
       exit(1);
     }
     if((num_front= (int *)calloc(nseq,sizeof(int))) == NULL)//worst case all reads are linked together
     {
       printf("ssaha: calloc - num_front\n");
       exit(1);
     }
     if((out_list= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - out_list\n");
       exit(1);
     }
     if((ctg_list= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - ctg_list\n");
       exit(1);
     }
     if((rd_name= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - rd_name\n");
       exit(1);
     }
     if((rd_used= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - rd_used\n");
       exit(1);
     }
     if((rd_rept= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - rd_rept\n");
       exit(1);
     }

/*   get the readnames for each contigs   */
     n_contig=0;
     n_reads=0;
     n_grouped=0;
     
     for(i=0;i<nseq;i++)
     {
        int id_pt,n_front,n_midd,n_end,n_over,NLINK=nmatch;
        int a_index,a_links;

        n_reads=0;
        n_over=0;
        n_front=n_list[i];
        if((rd_used[i]==0)&&(n_front>0))
        {
          for(j=0;j<n_front;j++)
          {
             a_index=r_matrix[i][j]>>rshift;
             a_links=r_matrix[i][j]&077;
             ind_front[j]=a_index;
             num_front[j]=a_links;
             if(a_links>=NLINK)
             {
               n_over++;
             }
          }
          if(n_over==0) 
            n_front=0;
          else
          {
            rd_used[i]=1;
            rd_name[n_reads]=i;
            n_reads++;
            rd_rept[i]=1;
            for(j=0;j<n_front;j++)
            {
               a_index=r_matrix[i][j]>>rshift;
               a_links=r_matrix[i][j]&077;
               if((a_links>=NLINK)&&(rd_used[a_index]==0))
               {
                 rd_name[n_reads]=a_index;
                 rd_rept[a_index]=1;
                 n_reads++;
               }
            }
          }
          while(n_front>0)
          {
            n_midd=0;
            for(j=0;j<n_front;j++)
            {
               id_pt=ind_front[j];
               if((num_front[j]>=NLINK)&&(rd_used[id_pt]==0))
               {
                 rd_used[id_pt]=1;
                 for(k=0;k<n_list[id_pt];k++)
                 {
                    a_index=r_matrix[id_pt][k]>>rshift;
                    a_links=r_matrix[id_pt][k]&077;
                    if((a_links>=NLINK)&&(rd_rept[a_index]==0)&&(rd_used[a_index]==0))
                    {
                      rd_name[n_reads]=a_index;
                      n_reads++;
                      rd_rept[a_index]=1;
                      num_front[n_front+n_midd]=a_links;
                      ind_front[n_front+n_midd]=a_index;
                      n_midd++;
                    }
                 }
               }
            }
            n_end=0;
            for(j=0;j<(n_front+n_midd);j++)
            {
               if((rd_used[ind_front[j]]==0)&&(num_front[j]>=NLINK))
               {
                 num_front[n_end]=num_front[j];
                 ind_front[n_end]=ind_front[j];
                 n_end++; 
               }
            }
            n_front=n_end;
          }
          if(n_reads>=2)
          {
            if(outmod==0)
              printf("readnames for contig %d %s %d \n",n_contig,names[i],n_reads);
            for(j=0;j<n_reads;j++)
            {
               if(outmod==0)
                 printf("%s\n",names[rd_name[j]]);
//                 printf("%s %d %d\n",names[rd_name[j]],rd_name[j],num_pair[rd_name[j]]);
               out_list[n_grouped]=rd_name[j];
               n_grouped++;
            }
            n_contig++;
            ctg_list[n_contig]=n_reads;
          }
//          printf("\n");
        }
     }


     if(outmod==0)
     {
       printf("All jobs finished\n");
       fflush(stdout);
       system("date");
       return;
     }
     free(*r_matrix);
     free(r_matrix);
     printf("Relation matrix finished\n");
     fflush(stdout);
     system("date");
     if((num_pair= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - num_pair\n");
       exit(1);
     }
     if((mat_pair= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - mat_pair\n");
       exit(1);
     }
     if((ins_pair= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - ins_pair\n");
       exit(1);
     }
     if((dev_pair= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - dev_pair\n");
       exit(1);
     }

/*   input mates information    */
     if((fpMate = fopen(argv[args],"r")) == NULL)
     {
       printf("Error ssaha: mate input\n");
       exit(1);
     }
     seqcharc=0;
     while(!feof(fpMate))
     {
       fgets(line,500,fpMate);
       seqcharc += strlen(line);
     }
     fclose(fpMate); 

     mates = (char **) malloc(nseq*sizeof(char *));
     mates[0] = (char *) calloc(seqcharc+nseq,1);     
     if((fpMate = fopen(argv[args],"r")) == NULL)
     {
       printf("Error ssaha: mate input\n");
       exit(1);
     }
     seqc = 0;
     while(!feof(fpMate))
     {
       int id_pt,nPair;

       fgets(line,500,fpMate);
       strcpy(mates[seqc],line);
       mates[seqc+1] = mates[seqc]+strlen(line)+1;
       nPair=0;
       strcpy(line2,line);
       for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
       {
       }
       i=0;
       for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
       {
          if(i==0)
          {
            strcpy(base,zero);
            strcat(base,ptr);
            id_pt = atoi(base);
          }
          else if(i==1)
          {
            strcpy(base,zero);
            strcat(base,ptr);
            num_pair[id_pt] = atoi(base);
          }
          else if(i==2)
          {
            strcpy(base,zero);
            strcat(base,ptr);
            mat_pair[id_pt] = atoi(base);
          }
          else if((nPair>2)&&(i==(nPair-2)))
          {
            strcpy(base,zero);
            strcat(base,ptr);
            ins_pair[id_pt] = atoi(base);
          }
          else if((nPair>2)&&(i==(nPair-1)))
          {
            strcpy(base,zero);
            strcat(base,ptr);
            dev_pair[id_pt] = atoi(base);
          }
       }
       seqc++;
     }
     fclose(fpMate); 

//     for(i=0;i<nseq;i++)
//        printf("%d %d %d %d %d\n",i,num_pair[i],mat_pair[i],ins_pair[i],dev_pair[i]);
//     for(i=0;i<nseq;i++)
//        printf("len: %d %s\n",i,mates[i]);

/*   output read names, fasta and fastq dat in every contig */ 
     nseq=0;
     n_input=0;
     for(ac=(args+1);ac<argc;ac++)
     {
	seq = decodeFastq ( argv[ac], &nSeq, &totalBases, qthresh);
	if(seq == NULL)
	{
		printf("ERROR ssaha: no data found\n");
		exit(1);
	}
	fastaUC(seq,nSeq);
        expp[n_input]=seq;
/*      rd_used and rd_rept are used again to store the read index */
        for(i=0;i<nSeq;i++)
        {
           rd_used[i+nseq]=n_input;
           rd_rept[i+nseq]=i;
        }
        nseq=nseq+nSeq;
        n_input++;
     }

     ctg_list[0]=0;
     for(i=1;i<=n_contig;i++)
        ctg_list[i]=ctg_list[i]+ctg_list[i-1]; 
     n_patch=0;
     n_outreads=0;
     new_patch=0;
     rd_st=0;
     rd_ed=0;
     n_ctg=0;
     for(i=1;i<=n_contig;i++)
     {
        n_outreads=n_outreads+ctg_list[i]-ctg_list[i-1];
        if(n_outreads>N_SET)
        {
          new_patch=1;
          rd_st=rd_ed+1;
          rd_ed=i;
          n_ctg=0;
          n_outreads=0;
        }       
        n_ctg++;

        if(new_patch)
        {
          sprintf(outName,"readname%04d",n_patch);
          if((fpOutname = fopen(outName,"w")) == NULL)
          {
            printf("Unable to open file for name out\n");
            exit(1);
          }
          sprintf(outFast,"%04d.fastq",n_patch);
          if((fpOutfast = fopen(outFast,"w")) == NULL)
          {
            printf("Unable to open file for fastq out\n");
            exit(1);
          }
          for(k=rd_st;k<=rd_ed;k++)
          {
             int id=out_list[ctg_list[k-1]];
           
             n_reads=ctg_list[k]-ctg_list[k-1];
             seqp=expp[rd_used[id]]+rd_rept[id];
             fprintf(fpOutname,"readnames for contig %d %s %d \n",k-1,seqp->name,n_reads);
             for(j=0;j<n_reads;j++)
             {
                int rdrc,rd=out_list[ctg_list[k-1]+j],rc,i_pair,np;
                int nPair,iPair;
 
                seqp=expp[rd_used[rd]]+rd_rept[rd];
                fprintf(fpOutfast,"@%s %s\n",seqp->name,seqp->name2);
                for(rc=0;rc<seqp->length;rc++)
                   fprintf(fpOutfast,"%c",seqp->data[rc]);
                fprintf(fpOutfast,"\n");
                fprintf(fpOutfast,"+%s %s\n",seqp->name,seqp->name2);
                putc(041,fpOutfast);
                for(rc=1;rc<seqp->length;rc++)
                {
                   if(seqp->finished)
                     putc(40+041,fpOutfast); 
                   else
                     putc(seqp->qual[rc]+041,fpOutfast); 
                }
                fprintf(fpOutfast,"\n");
                if(num_pair[rd]==0)
                {
                  fprintf(fpOutname,"%s\n",seqp->name);
                }
                else if(num_pair[rd]>0)
                {
                  i_pair=0;
                  if(num_pair[rd]==1)
                  {
                    read_pair[i_pair]=mat_pair[rd];
                    i_pair++;
                  }
                  else
                  {
                    strcpy(line,mates[rd]);
                    strcpy(line2,mates[rd]);
                    nPair=0;
                    iPair=0;
                    for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
                    {
                    }
                    for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),iPair++)
                    {
                       if((iPair>1)&&(iPair<(nPair-2)))
                       {
                         strcpy(base,zero);
                         strcat(base,ptr);
                         read_pair[i_pair]=atoi(base);
                         i_pair++;
                       }
                    }
                    
                  }
                  fprintf(fpOutname,"%s ",seqp->name);
                  for(np=0;np<i_pair;np++)
                  {
                     rdrc=read_pair[np];
                     seqr=expp[rd_used[rdrc]]+rd_rept[rdrc];
                     fprintf(fpOutfast,"@%s %s\n",seqr->name,seqr->name2);
                     for(rc=0;rc<seqr->length;rc++)
                        fprintf(fpOutfast,"%c",seqr->data[rc]);
                     fprintf(fpOutfast,"\n");
                     fprintf(fpOutfast,"+%s %s\n",seqr->name,seqr->name2);
                     putc(041,fpOutfast);
                     for(rc=1;rc<seqr->length;rc++)
                     {
                        if(seqr->finished)
                          putc(40+041,fpOutfast); 
                        else
                          putc(seqr->qual[rc]+041,fpOutfast); 
                     }
                     fprintf(fpOutfast,"\n");
                     fprintf(fpOutname,"%s ",seqr->name);
                  }
                  fprintf(fpOutname,"%d %d\n",ins_pair[rd],dev_pair[rd]);
                }
             }
          }  
          n_patch++;
          fclose(fpOutname);
          fclose(fpOutfast);
          sprintf(syscmd,"gzip %s",outName);
          system (syscmd);
          sprintf(syscmd,"gzip %s",outFast);
          system (syscmd);
          new_patch=0;
        }
     }
     if(rd_ed<n_contig) 
     {
          sprintf(outName,"readname%04d",n_patch);
          if((fpOutname = fopen(outName,"w")) == NULL)
          {
            printf("Unable to open file for name out\n");
            exit(1);
          }
          sprintf(outFast,"%04d.fastq",n_patch);
          if((fpOutfast = fopen(outFast,"w")) == NULL)
          {
            printf("Unable to open file for fastq out\n");
            exit(1);
          }
          for(k=rd_ed+1;k<=n_contig;k++)
          {
             int id=out_list[ctg_list[k-1]];
           
             n_reads=ctg_list[k]-ctg_list[k-1];
             seqp=expp[rd_used[id]]+rd_rept[id];
             fprintf(fpOutname,"readnames for contig %d %s %d \n",k-1,seqp->name,n_reads);
             for(j=0;j<n_reads;j++)
             {
                int rdrc,rd=out_list[ctg_list[k-1]+j],rc,i_pair,np;
                int nPair,iPair;
 
                seqp=expp[rd_used[rd]]+rd_rept[rd];
                fprintf(fpOutfast,"@%s %s\n",seqp->name,seqp->name2);
                for(rc=0;rc<seqp->length;rc++)
                   fprintf(fpOutfast,"%c",seqp->data[rc]);
                fprintf(fpOutfast,"\n");
                fprintf(fpOutfast,"+%s %s\n",seqp->name,seqp->name2);
                putc(041,fpOutfast);
                for(rc=1;rc<seqp->length;rc++)
                {
                   if(seqp->finished)
                     putc(40+041,fpOutfast); 
                   else
                     putc(seqp->qual[rc]+041,fpOutfast); 
                }
                fprintf(fpOutfast,"\n");
                if(num_pair[rd]==0)
                {
                  fprintf(fpOutname,"%s\n",seqp->name);
                }
                else if(num_pair[rd]>0)
                {
                  i_pair=0;
                  if(num_pair[rd]==1)
                  {
                    read_pair[i_pair]=mat_pair[rd];
                    i_pair++;
                  }
                  else
                  {
                    strcpy(line,mates[rd]);
                    strcpy(line2,mates[rd]);
                    nPair=0;
                    iPair=0;
                    for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
                    {
                    }
                    for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),iPair++)
                    {
                       if((iPair>1)&&(iPair<(nPair-2)))
                       {
                         strcpy(base,zero);
                         strcat(base,ptr);
                         read_pair[i_pair]=atoi(base);
                         i_pair++;
                       }
                    }
                  }
                  fprintf(fpOutname,"%s ",seqp->name);
                  for(np=0;np<i_pair;np++)
                  {
                     rdrc=mat_pair[rd];
                     seqr=expp[rd_used[rdrc]]+rd_rept[rdrc];
                     fprintf(fpOutfast,"@%s %s\n",seqr->name,seqr->name2);
                     for(rc=0;rc<seqr->length;rc++)
                        fprintf(fpOutfast,"%c",seqr->data[rc]);
                     fprintf(fpOutfast,"\n");
                     fprintf(fpOutfast,"+%s %s\n",seqr->name,seqr->name2);
                     putc(041,fpOutfast);
                     for(rc=1;rc<seqr->length;rc++)
                     {
                        if(seqr->finished)
                          putc(40+041,fpOutfast); 
                        else
                          putc(seqr->qual[rc]+041,fpOutfast); 
                     }
                     fprintf(fpOutfast,"\n");
                     fprintf(fpOutname,"%s ",seqr->name);
                  }
                  fprintf(fpOutname,"%d %d\n",ins_pair[rd],dev_pair[rd]);
                }
             }
          }  
          n_patch++;
          fclose(fpOutname);
          fclose(fpOutfast);
          sprintf(syscmd,"gzip %s",outName);
          system (syscmd);
          sprintf(syscmd,"gzip %s",outFast);
          system (syscmd);
     }
     printf("All jobs finished\n");
     fflush(stdout);
     system("date");
     return;
     {
       int j = (int)sarr[i]&nmask;
       B64_long ii=(B64_long)sarr[i]>>sshift;
       B64_long mask=3;
       B64_long jj;
       char *base = "ACGT";
       for(jj=n_Sbase;--jj>=0;)
   	  printf("%c",base[(ii>>(jj<<1))&mask]);
       printf(" %-30s\n",names[j]);
     }
     seqp = seq;
     for(i=0;i<nSeq;i++,seqp++)
     {
        b = (char*)seqp->data;
        ns=(seqp->length)-n_Sbase;
        for(j=0;j<ns;j++)
        {
           IntSeg=0;
           ch=b;
           for(k=n_Sbase;k>0;k--)
           {
              if     (*ch == 'A') IntBase=0;
              else if(*ch == 'C') IntBase=1;
              else if(*ch == 'G') IntBase=2;
              else if(*ch == 'T') IntBase=3;
              else
              {
                IntBase=1;
              }
              IntBase=IntBase<<((k-1)<<1);
              IntSeg=IntSeg+IntBase;
              ch++;
           }
		   printf("%c %d\n",*b,head[IntSeg]);
           b++;
        }
     }

}


int ssaha_init( char **argv, int argc, B64_long dataSize, int IMOD, int ISUB, int NCUT, int NHIT, int SEG_LEN, int mhisti)
{
    FILE *fpSM;
    int i,i_Fasta,sizesm;
    int ac = 3;

    mhist = mhisti;
    n_Entry=1L<<(SEG_LEN<<1);

/*  allocate memory for arrays HEAD and LIST   */
	
    if((head = (signed char *)calloc((n_Entry),sizeof(signed char))) == NULL)
    {
      printf("ERROR ssaha_init: calloc - head\n");
      exit(1);
    }
    if((hist = (B64_long *)calloc(400,sizeof(B64_long))) == NULL)
    {
      printf("ERROR ssaha_init: calloc - hist\n");
      exit(1);
    }

    if(IMOD==0)
    {
      /* define the sequence matrix SM     */

      i_Fasta=1;
      fastaSort(argv,argc,i_Fasta, SEG_LEN);
    }
  return(1);   
}

/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
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
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

