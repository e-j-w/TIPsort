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
		    if(e>0)
		      if(e<S32K)
			hist[pos][col][(int)rint(e)]++;
  		  }
  free(cev);
  return SEPARATOR_DISCARD;
}
/*====================================================================================*/
int main(int argc, char *argv[])
{
  FILE * output;
  input_names_type* name;
  int stop,pos,col;
  char n[132];

  if(argc!=2)
    {
      printf("\n ./GRIFFIN_ECal master_file_name\n");
      exit(-1);
    }
  
  printf("Program sorts calibrated energy histogram for GRIFFIN cores \n");
  name=(input_names_type*)malloc(sizeof(input_names_type));
  memset(name,0,sizeof(input_names_type));
  cal_par=(calibration_parameters*)malloc(sizeof(calibration_parameters));
  memset(cal_par,0,sizeof(calibration_parameters));
  memset(hist,0,sizeof(hist));
  read_master(argv[1],name);

  if(name->flag.inp_data!=1)
    {
      printf("\nInput data file not defined\n");
      exit(EXIT_FAILURE);
    }

  if(name->flag.GRIFFIN_cal_par==1)
        {
          printf("\nGRIFFIN calibration read from the file:\n %s\n",name->fname.GRIFFIN_cal_par);
          initialize_GRIFFIN_calibration(&cal_par->gr,name->fname.GRIFFIN_cal_par);
	  
        }
      else
        {
          printf("\nGRIFFIN calibration parameters not defined\n");
          exit(EXIT_FAILURE);
        }


  summarize_GRIFFIN_calibration(&cal_par->gr,"in.in");
  sort(name);

  summarize_GRIFFIN_calibration(&cal_par->gr,"out.out");
  for(pos=1;pos<NPOSGRIF;pos++)
    {
      stop=0;
      for(col=0;col<NCOL;col++)
	stop+=cal_par->gr.ceflag[pos][col];
      if(stop>0)
	{
	  sprintf(n,"pos%1d%1d_cal_core_energy.mca",pos/10,pos-(pos/10)*10);
	  if((output=fopen(n,"w"))==NULL)
	    {
	      printf("\nI can't open file %s\n",n);
	      exit(EXIT_FAILURE);
	    }
	  
	  for(col=0;col<NCOL;col++)
	    {
	      fwrite(hist[pos][col],S32K*sizeof(int),1,output);
	    }
	  fclose(output);
	}
    }

 
}

