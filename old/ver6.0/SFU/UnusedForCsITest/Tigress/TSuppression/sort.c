#include "sort.h"

int analyze_data(raw_event *data)
{
  cal_event* cev;
  int col,sup;
  double t;
  cev=(cal_event*)malloc(sizeof(cal_event));
  memset(cev,0,sizeof(cal_event));
  calibrate_TIGRESS(data,&cal_par->tg,&cev->tg);
  if(cev->tg.h.FE>0)
    if((cev->tg.h.EHP&(1<<(pos-1)))!=0)
      if(cev->tg.det[pos].hge.FE>0)
	if(data->tg.h.BGOfold>0) /* BGO in the event */
	  if((data->tg.h.BGOHP&(1<<(pos-1)))!=0) /* BGO at the right position in the event */
	    if(data->tg.det[pos].h.BGOfold>0) /* BGO at the right position in the event */
	      for(col=0;col<NCOL;col++) /* loop over colors */
		if((data->tg.det[pos].h.BGOHP&(1<<col))!=0)  /* BGO at the right position and color in the event */
		  if(data->tg.det[pos].bgo[col].h.Tfold>0)   /* BGO time information at the right position and color in the event */
		    for(sup=0;sup<NSUP;sup++) /* loop over BGO suppressors */
		      if((data->tg.det[pos].bgo[col].h.THP&(1<<sup))!=0) 
			{
			  /* BGO time */
			  t=data->tg.det[pos].bgo[col].sup[sup].cfd&0x00ffffff;
			  t-=(data->tg.det[pos].bgo[col].sup[sup].timestamp*16)&0x00ffffff;
			  if((data->h.setupHP&RF_BIT)!=0)
			    {
			      t-=(int)data->rf.sin.t0;
			      t+=S16K;
			      if(t<0)t=S32K-2;
			      if(t>S32K) t=S32K-3;
			    }
			  else
			    t=S32K-4;

			  if(t>=0)
			    if(t<S32K)
			      hist[col*NSUP+sup][(int)rint(t)]++;
			
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

  if(argc!=5)
    {
      printf("\n ./Tigress_Suppression master_file_name position tlow thigh\n");
      exit(-1);
    }
  
  printf("Program sorts raw time histogram for Tigress suppressors \n");

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

  pos=atoi(argv[2]);
  tlow=atoi(argv[3]);
  thigh=atoi(argv[4]);
  sort(name);


  sprintf(n,"pos%1d%1d_suppresion.mca",pos/10,pos-(pos/10)*10);

  if((output=fopen(n,"w"))==NULL)
    {
      printf("\nI can't open file %s\n",n);
      exit(EXIT_FAILURE);
    }
  for(k=0;k<NCOL*NSUP;k++)
    fwrite(hist[k],S32K*sizeof(int),1,output);
  fclose(output);




}

