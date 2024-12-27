#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>


#include "fasta.h"
#include "ssaha.h"
         
int ssaha(char *queryfile, char *subjectFile, int IMOD, int ISUB, int NCUT, int NHIT, int SEG_LEN)
{
    FILE *fpFasta,*fpSM;
    char dname[60]; 
    int nSeq,nSeg,batch=0;
    fasta *seg, *seq;
    DIR   *dp;
    struct dirent *dep;
    HitInfo *hip;
    int n_hip;
    int n_hip_max = 2000;
	long totalBases;

	seg = decodeFastq(subjectFile, &nSeg, &totalBases);
    fastaUC(seg,nSeg);
    ssaha_init(seg, nSeg, totalBases, IMOD, ISUB, NCUT, NHIT, SEG_LEN);
    if((hip = (HitInfo *) calloc(n_hip_max,sizeof(HitInfo))) == NULL)
    {
		printf("ssaha: calloc error\n");
		exit(0);
    }

    if (IMOD==0 || IMOD==1 || IMOD==3 || IMOD==6)
    {
/*    read the test sequence and search the whole data base of SM   */
      seq = decodeFastq(queryfile,&nSeq,&totalBases);
      fastaUC(seq,nSeq);
      SM_Info(NCUT);
      for(batch=0;batch<nSeq;batch++) 
      {
		n_hip = n_hip_max;
		Search_SM(seq+batch,IMOD,ISUB,NCUT,NHIT,SEG_LEN,NULL,hip,&n_hip);
/*      print out the gene names which have been found    */
  printf("Input Query Sequence Name: %s %-35.35s\n ",seq[batch].name,seq[batch].name2);
  printf("----------------------------------------------------------------\n");
         HitPrint(n_hip, hip);
      }

    }
    else if (IMOD==2||IMOD==4)
    {
/*    read the test sequence and search the whole data base of SM   */
      seq = decodeFastq(queryfile,&nSeq,&totalBases);
      fastaUC(seq,nSeq);
      SM_Info(NCUT);
      if((fpSM=fopen("sm.out","rb"))==NULL)
      {
		printf("Cannot open file: sm.out\n");
		exit(0);
      }
      for(batch=0;batch<nSeq;batch++) 
      {
		n_hip = n_hip_max;
		Search_SM(seq+batch,IMOD,ISUB,NCUT,NHIT,SEG_LEN,NULL,hip,&n_hip);
/*      print out the gene names which have been found    */
  printf("Input Query Sequence Name: %s %-35.35s\n ",seq[batch].name,seq[batch].name2);
  printf("----------------------------------------------------------------\n");
         HitPrint(n_hip, hip);
      }
      fclose(fpSM);
    }
    else
    {
      printf("ERROR ssaha: input operation mode?\n");
      exit(1);
    }
   
    printf("Job finished for %d input sequences!\n",nSeq);
  return(0);
}

int main(int argc, char **argv)
{
    int i;

    /* SSAHA default parameters    */
    int SEG_LEN=12;
    int IMOD=0;
    int ISUB=0;
    int NCUT=200;
    int NHIT=18;

    if(argc < 3)
    {
      printf("Usage: Not enough premeters to run SSAHA\n");
      exit(1);
    }

    for(i=3;i<argc;i++)
    {
       if(!strcmp(argv[i],"-mod"))
         sscanf(argv[++i],"%d",&IMOD); 
       else if(!strcmp(argv[i],"-sub"))
         sscanf(argv[++i],"%d",&ISUB); 
       else if(!strcmp(argv[i],"-cut"))
         sscanf(argv[++i],"%d",&NCUT); 
       else if(!strcmp(argv[i],"-hit"))
         sscanf(argv[++i],"%d",&NHIT); 
       else if(!strcmp(argv[i],"-len"))
         sscanf(argv[++i],"%d",&SEG_LEN);
       else
       {
         printf("input parameter error\n");
         exit(1);
       }
    }
    ssaha(argv[1], argv[2], IMOD, ISUB, NCUT, NHIT, SEG_LEN);
    return(0);
}
/* end of the main */
