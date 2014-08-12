#include "sort.h"

int analyze_data(raw_event *data)
{
  cal_event* cev;
  int eBGO,suppFlag,eAddBack;
  
  cev=(cal_event*)malloc(sizeof(cal_event));
  memset(cev,0,sizeof(cal_event));
  calibrate_TIGRESS(data,&cal_par->tg,&cev->tg);

  if(cev->tg.h.FH>0)
    for(pos=1;pos<NPOSTIGR;pos++)
      {
	suppFlag=0;
	if(((cev->tg.h.HHP&(1<<(pos-1)))!=0) && (cev->tg.det[pos].hge.FE>0))
	  if((cev->tg.h.AHP&(1<<(pos-1)))!=0)
	    {
	      eAddBack = cev->tg.det[pos].addback.E/cal_par->tg.contr_e;
	      if(data->tg.h.BGOfold>0) /* BGO in the event */
		if(((data->tg.h.BGOHP&(1<<(pos-1))!=0)) && ((data->tg.det[pos].h.BGOfold)>0)) //At the right position in the event
		  for(col=0;col<NCOL;col++) /* loop over colors */
		    if(((data->tg.det[pos].h.BGOHP&(1<<col))!=0) && (data->tg.det[pos].bgo[col].h.Efold>0)) //At the right pos/col with energy
		      for(sup=0;sup<NSUP;sup++) /* loop over BGO suppressors */
			if(((data->tg.det[pos].bgo[col].h.EHP&(1<<sup))!=0) && ((data->tg.det[pos].bgo[col].h.Efold)>0))
			  {
			    eBGO=(int)rint(data->tg.det[pos].bgo[col].sup[sup].charge);
			    if((eBGO>elow) && (eBGO<ehigh)) suppFlag=1;
			  }
	      if((eAddBack>=0) && (eAddBack<S32K)) hist[pos][suppFlag][eAddBack]++;
	    }}
  free(cev);
  return SEPARATOR_DISCARD;
}
/*================================================================================*/
int main(int argc, char *argv[])
{
  FILE * output;
  input_names_type* name;
  char mcaFile[132];

  if(argc!=4)
    {
      printf("Tigress_ECalABSupp_OnBGOE master_file_name elow ehigh\n");
      exit(-1);
    }
  
  printf("Program sorts ECalABSupp spectra for TIGRESS using BGO energy regions.\n");
  
  name=(input_names_type*)malloc(sizeof(input_names_type));
  memset(name,0,sizeof(input_names_type));
  
  cal_par=(calibration_parameters*)malloc(sizeof(calibration_parameters));
  memset(cal_par,0,sizeof(calibration_parameters));
  
  memset(hist,0,sizeof(hist));
  read_master(argv[1],name);
  
  elow=atoi(argv[2]);
  ehigh=atoi(argv[3]);

  if(name->flag.inp_data!=1)
    {
      printf("ERROR!!! Input data file not defined!\n");
      exit(EXIT_FAILURE);
    }
  
  if(name->flag.TIGRESS_cal_par==1)
    {
      printf("TIGRESS calibration read from %s\n",name->fname.TIGRESS_cal_par);
      initialize_TIGRESS_calibration(&cal_par->tg,name->fname.TIGRESS_cal_par);
    }
  else
    {
      printf("ERROR!!! TIGRESS calibration parameters not defined!\n");
      exit(EXIT_FAILURE);
    }
  
  sort(name);

  for(pos=1;pos<NPOSTIGR;pos++)
    {
      sprintf(mcaFile,"pos%1d%1d_ECalABSupp_OnBGOE.mca",pos/10,pos-(pos/10)*10);
      if((output=fopen(mcaFile,"w"))==NULL)
	{
	  printf("ERROR!!! I cannot open the mca file!\n");
	  exit(EXIT_FAILURE);
	}
      
      for(int k=0;k<2;k++) fwrite(hist[pos][k],S32K*sizeof(int),1,output);
      fclose(output);
    }}
