#include "sort.h"

int analyze_data(raw_event *data)
{
  cal_event* cev;
  int ene;
  cev=(cal_event*)malloc(sizeof(cal_event));
  memset(cev,0,sizeof(cal_event));
  calibrate_TIGRESS(data,&cal_par->tg,&cev->tg);

  if(cev->tg.h.FE>0)
    for(pos=1;pos<NPOSTIGR;pos++)
      if((cev->tg.h.AHP&(1<<(pos-1)))!=0)
	{	  
	  ene=(int)rint(cev->tg.det[pos].addback.E/cal_par->tg.contr_e);
	  if(ene>0)
	    if(ene<S32K)
	      hist[pos][ene]++;
	}
  free(cev);

  return SEPARATOR_DISCARD;
}
/*====================================================================================*/
int main(int argc, char *argv[])
{
  FILE * output;
  input_names_type* name;
  char n[132];
  int k;

  if(argc!=2)
    {
      printf("\n ./Tigress_SuppRFAddBack master_file_name\n");
      exit(-1);
    }
  
  printf("Program sorts BGO suppressed Tigress addback spectra \n");

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

  if(name->flag.TIGRESS_cal_par==1)
        {
          printf("\nTIGRESS calibration read from the file:\n %s\n",name->fname.TIGRESS_cal_par);
          initialize_TIGRESS_calibration(&cal_par->tg,name->fname.TIGRESS_cal_par);
	  
        }
      else
        {
          printf("\nTIGRESS calibration parameters not defined\n");
          exit(EXIT_FAILURE);
        }


 
  sort(name);
 
  sprintf(n,"SuppRFAddBack.mca");

  if((output=fopen(n,"w"))==NULL)
    {
      printf("\nI can't open file %s\n",n);
      exit(EXIT_FAILURE);
    }
  for(k=0;k<NPOSTIGR;k++)
    fwrite(hist[k],S32K*sizeof(int),1,output);
  fclose(output);




}

