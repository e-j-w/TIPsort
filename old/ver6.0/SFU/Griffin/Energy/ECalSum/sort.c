#include "sort.h"

int analyze_data(raw_event *data)
{
  cal_event* cev;
  int pos,col;
  double e;
  
  cev=(cal_event*)malloc(sizeof(cal_event));
  memset(cev,0,sizeof(cal_event));
  calibrate_GRIFFIN(data,&cal_par->gr,&cev->gr);
  if(cev->gr.h.FE>0)
    for(pos=1;pos<NPOSGRIF;pos++)
      if((cev->gr.h.EHP&(1<<(pos-1)))!=0)
  	if(cev->gr.det[pos].hge.FE>0)
  	  for(col=0;col<NCOL;col++)
  	    if((cev->gr.det[pos].hge.EHP&(1<<col))!=0)
  	      if(cev->gr.det[pos].ge[col].h.FE>0)
  		if((cev->gr.det[pos].ge[col].h.EHP&1)!=0)
  		  {
  		    e=cev->gr.det[pos].ge[col].seg[0].E/cal_par->gr.contr_e;
		    if(e>=0 && e<S32K)
		      hist[pos][col][(int)rint(e)]++;
		    else hist[pos][col][S32K-1000];
  		  }
  free(cev);
  return SEPARATOR_DISCARD;
}
/*====================================================================================*/
int main(int argc, char *argv[])
{
  FILE *output, *cluster;
  input_names_type* name;
  int stop,pos,col;
  char n[132], DataFile[132];
  
  if(argc!=2)
    {
      printf("Griffin_ECalSum master_file_name\n");
      exit(-1);
    }
  
  printf("Program sorts ECalSum energy histogram for GRIFFIN cores.\n");

  name=(input_names_type*)malloc(sizeof(input_names_type));
  memset(name,0,sizeof(input_names_type));
  cal_par=(calibration_parameters*)malloc(sizeof(calibration_parameters));
  memset(cal_par,0,sizeof(calibration_parameters));
  memset(hist,0,sizeof(hist));
  read_master(argv[1],name);

  if(name->flag.cluster_file==1)
    {
      printf("Sorting ECalSum energy histograms for GRIFFIN cores based upon the cluster file: %s\n",name->fname.cluster_file);
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
  
  if(name->flag.GRIFFIN_cal_par==1)
    {
      printf("GRIFFIN calibration read from: %s\n",name->fname.GRIFFIN_cal_par);
      initialize_GRIFFIN_calibration(&cal_par->gr,name->fname.GRIFFIN_cal_par);
    }
  else
    {
      printf("ERROR! GRIFFIN calibration parameters not defined!\n");
      exit(EXIT_FAILURE);
    }
  
  summarize_GRIFFIN_calibration(&cal_par->gr,"in.in");  
  
  while(fscanf(cluster,"%s",DataFile) != EOF)
    {
      memset(name,0,sizeof(input_names_type));
      strcpy(name->fname.inp_data,DataFile);
      
      printf("Sorting gamma-ray energy data from %s\n", name);
      sort(name);
    }
  
  for(pos=1;pos<NPOSGRIF;pos++)
    {
      stop=0;
      for(col=0;col<NCOL;col++)
	stop+=cal_par->gr.ceflag[pos][col];
      if(stop>0)
	{
	  sprintf(n,"pos%1d%1d_ECalSum.mca",pos/10,pos-(pos/10)*10);
	  if((output=fopen(n,"w"))==NULL)
	    {
	      printf("ERROR! I can't open %s!\n",n);
	      exit(EXIT_FAILURE);
	    }
	  
	  for(col=0;col<NCOL;col++)
	    fwrite(hist[pos][col],S32K*sizeof(int),1,output);
	    
	  fclose(output);
	}}}
