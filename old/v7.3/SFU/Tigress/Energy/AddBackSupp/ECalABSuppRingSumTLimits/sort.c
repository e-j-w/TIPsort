#include "sort.h"

int analyze_data(raw_event *data)
{
  cal_event* cev;
  double eAddBack=0,tAddBack=0,tForGating=0;
  int suppFlag=0;
  
  cev=(cal_event*)malloc(sizeof(cal_event));
  memset(cev,0,sizeof(cal_event));
  calibrate_TIGRESS(data,&cal_par->tg,&cev->tg);
  
  if(cev->tg.h.FH>0)
    for(pos=1;pos<NPOSTIGR;pos++)
      if(((cev->tg.h.HHP&(1<<(pos-1)))!=0) && (cev->tg.det[pos].hge.FH>0))
	if((cev->tg.h.AHP&(1<<(pos-1)))!=0)
	  {
	    //if(cev->tg.det[pos].suppress==0)
	    //time add back is average of cores participating in event
	    eAddBack   = cev->tg.det[pos].addback.E/cal_par->tg.contr_e;;
	    tAddBack   = cev->tg.det[pos].addback.T;
	    suppFlag   = cev->tg.det[pos].suppress;
	    col        = cev->tg.det[pos].addbackC;
	    ring       = cev->tg.det[pos].ge[col].ring;

	    //for the core of first interaction (the one with highest energy)
	    //make sure the time is within time gating
	    tForGating = cev->tg.det[pos].ge[col].seg[0].T;
	
	    if((tForGating>cal_par->tg.ctlow[pos][col]) && (tForGating<cal_par->tg.cthigh[pos][col]))
	      {
		if((eAddBack>0) && (eAddBack<S32K))
		  {
		    if((ring>0) && (ring<NRING))
		      {
			ring += NRING*suppFlag;
			hist[ring][(int)rint(eAddBack)]++;
		      }}
		else hist[ring][S32K-1000];
	      }
	    else hist[ring][S32K-2000];
	  }
  free(cev);
  return SEPARATOR_DISCARD;
}
/*====================================================================================*/
int main(int argc, char *argv[])
{
  input_names_type* name;
  FILE* output,*cluster;
  char DataFile[132];
  
  if(argc!=2)
    {
      printf("Tigress_ECalABSuppRingSumTLimits master_file_name\n");
      exit(-1);
    }
    
  printf("Program sorts TIGRESS ECalABSuppRingSum spectra gated on times.\n");

  name=(input_names_type*)malloc(sizeof(input_names_type));
  memset(name,0,sizeof(input_names_type));

  cal_par=(calibration_parameters*)malloc(sizeof(calibration_parameters));
  memset(cal_par,0,sizeof(calibration_parameters));

  memset(hist,0,sizeof(hist));
  read_master(argv[1],name);
  
  if(name->flag.cluster_file==1)
    {
      printf("Sorting data from cluster file: %s\n",name->fname.cluster_file);
      if((cluster=fopen(name->fname.cluster_file,"r"))==NULL)
	{
	  printf("ERROR!!! I cannot open input file %s!\n",name->fname.cluster_file);
	  exit(-2);
	}}
  else
    {
      printf("ERROR!!! Cluster file not defined!\n");
      exit(-1);
    }
  
  if(name->flag.TIGRESS_cal_par==1)
    {
      printf("TIGRESS calpar read from: %s\n",name->fname.TIGRESS_cal_par);
      initialize_TIGRESS_calibration(&cal_par->tg,name->fname.TIGRESS_cal_par);
    }
  else
    {
      printf("ERROR!!! TIGRESS calibration parameters not defined\n");
      exit(EXIT_FAILURE);
    }
  
  while(fscanf(cluster,"%s",DataFile) != EOF)
    {
      memset(name,0,sizeof(input_names_type));
      strcpy(name->fname.inp_data,DataFile);
      
      printf("Sorting data from file %s\n", name);
      sort(name);
    }
  
  if((output=fopen("Ring_ECalABSuppSumTLimits.mca","w"))==NULL)
    {
      printf("ERROR!!! I cannot open the mca file!\n");
      exit(EXIT_FAILURE);
    }
  fwrite(hist,2*NRING*S32K*sizeof(int),1,output);
  fclose(output);  
}

