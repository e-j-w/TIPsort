#include "sort.h"

int analyze_data(raw_event *data)
{
  cal_event* cev;
  double t;
  int tt;
  
  cev=(cal_event*)malloc(sizeof(cal_event));
  memset(cev,0,sizeof(cal_event));
  calibrate_TIGRESS(data,&cal_par->tg,&cev->tg);
  
  if(cev->tg.h.FT>0)
    for(pos=1;pos<NPOSTIGR;pos++)
      if((cev->tg.h.THP&(1<<(pos-1)))!=0)
	if(cev->tg.det[pos].hge.FT>0)
	  for(col=0;col<NCOL;col++)
	    if((cev->tg.det[pos].hge.THP&(1<<col))!=0)
	      if(cev->tg.det[pos].ge[col].h.FT>0)
		if((cev->tg.det[pos].ge[col].h.THP&1)!=0)
		  {
		    //In CFD units, then scaled by contraction
		    t=cev->tg.det[pos].ge[col].seg[0].T/cal_par->tg.contr_t;
		    //In ns
		    t*=0.625;
		    
		    tt=(int)rint(t);
		    if(tt<0) tt=S32K-1000;
		    if(tt>S32K) tt=S32K-2000;
		    hist[pos][col][tt]++;
		  }
  free(cev);
  return SEPARATOR_DISCARD;
}
/*=============================================================*/
int main(int argc, char *argv[])
{
  FILE *output, *cluster;
  input_names_type *name;
  char n[132], DataFile[256];
  int pos,col,stop;

  if(argc!=2)
    {
      printf("Tigress_TCalSum master_file_name\n");
      exit(-1);
    }
  
  printf("Program sorts summed TCal histograms for Tigress crystal cores from a cluster file\n");
  
  name=(input_names_type*)malloc(sizeof(input_names_type));
  cal_par=(calibration_parameters*)malloc(sizeof(calibration_parameters));
  memset(cal_par,0,sizeof(calibration_parameters));
  memset(hist,0,sizeof(hist));
  read_master(argv[1],name);
  
  if(name->flag.cluster_file==1)
    {
      printf("Sorting calibrated time histograms for TIGRESS clovers and cores based upon the cluster file: %s\n",name->fname.cluster_file);
      if((cluster=fopen(name->fname.cluster_file,"r"))==NULL)
	{
	  printf("ERROR! I can't open input file %s\n",name->fname.cluster_file);
	  exit(-2);
	}}
  else
    {
      printf("ERROR! Cluster file not defined\n");
      exit(-1);
    }

  if(name->flag.TIGRESS_cal_par==1)
    {
      printf("TIGRESS calibration read from the file: %s\n",name->fname.TIGRESS_cal_par);
      initialize_TIGRESS_calibration(&cal_par->tg,name->fname.TIGRESS_cal_par);	  
    }
  else
    {
      printf("ERROR! TIGRESS calibration parameters not defined\n");
      exit(EXIT_FAILURE);
    }
  
 while(fscanf(cluster,"%s",DataFile) != EOF)
    {
      memset(name,0,sizeof(input_names_type));
      strcpy(name->fname.inp_data,DataFile);
      
      printf("Sorting timing (CFD algorithm) data from file %s\n",DataFile);
      sort(name);
    }
  
  for(pos=1;pos<NPOSTIGR;pos++)
    {
      stop=0;
      for(col=0;col<NCOL;col++)
	stop+=cal_par->tg.ceflag[pos][col];
      if(stop>0)
	{
	  sprintf(n,"pos%1d%1d_CalSumCoreTime.mca",pos/10,pos-(pos/10)*10);
	  if((output=fopen(n,"w"))==NULL)
	    {
	      printf("ERROR! I can't open file %s\n",n);
	      exit(EXIT_FAILURE);
	    }
	  for(col=0;col<NCOL;col++)
	    fwrite(hist[pos][col],S32K*sizeof(int),1,output);
	  fclose(output);
	}}
 
  printf("The event-dropper increments are %d %d %d %d %d %d.\n",count_1,count_2,count_3,count_4,count_5,count_6);
}
