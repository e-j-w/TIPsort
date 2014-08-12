#include "sort.h"

int analyze_data(raw_event *data)
{
  cal_event* cev;
  unsigned long long int one=1;
  int pos,col,pin,ring;
  double ttg,tpin,epin,etg;
  
  cev=(cal_event*)malloc(sizeof(cal_event));
  memset(cev,0,sizeof(cal_event));
  calibrate_TIGRESS(data,&cal_par->tg,&cev->tg);
  calibrate_PINARRAY(data,&cal_par->pinarray,&cev->pinarray);

  if(cev->tg.h.FH>0)
    if(cev->pinarray.h.FH>0)
      for(pos=1;pos<NPOSTIGR;pos++)
	if((cev->tg.h.HHP&(1<<(pos-1)))!=0)
	  if(cev->tg.det[pos].hge.FH>0)
	    if(cev->tg.det[pos].suppress==0)
	    for(col=0;col<NCOL;col++)
	      if((cev->tg.det[pos].hge.HHP&(1<<col))!=0)
		if(cev->tg.det[pos].ge[col].h.FH>0)
		  if((cev->tg.det[pos].ge[col].h.HHP&1)!=0)
		    {
		      ttg=cev->tg.det[pos].ge[col].seg[0].T;
		      etg=cev->tg.det[pos].ge[col].seg[0].E;
		      ring=cev->tg.det[pos].ge[col].ring;
		      if(ttg>cal_par->tg.ctlow[pos][col])
			if(ttg<cal_par->tg.cthigh[pos][col])
			  if(etg>cal_par->tg.celow[pos][col])
			    if(etg<cal_par->tg.cehigh[pos][col])
			      for(pin=1;pin<NPIN;pin++)
				if((cev->pinarray.h.HHP&(one<<pin))!=0)
				  {
				    tpin=cev->pinarray.pin[pin].T;
				    epin=cev->pinarray.pin[pin].E;
				    if(tpin>cal_par->pinarray.tlow[pin])
				      if(tpin<cal_par->pinarray.thigh[pin])
					if(epin>cal_par->pinarray.elow[pin])
					  if(epin<cal_par->pinarray.ehigh[pin])
					    {
					      etg/=cal_par->tg.contr_e;
					      if(etg>0)
						if(etg<S4K)
						  if(ring>0)
						    if(ring<NRING)
						      hist[ring][(int)rint(etg)]++;

					    }
				       
			      }
		    }
  free(cev);
  return SEPARATOR_DISCARD;
}
/*====================================================================================*/
int main(int argc, char *argv[])
{
  input_names_type* name;
  FILE* output;
 
  if(argc!=2)
    {
      printf("\n ./Tigress_PINARRAY_TigressRing_AllLim master_file_name\n");
      exit(-1);
    }
  

  printf("Program sorts calibrated TIGRESS ring spectra gated on PINARRAY and TIGRESS energies and times \n");
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
          exit(-1);
        }

  if(name->flag.PINARRAY_cal_par==1)
        {
          printf("\nPINARRAY calibration read from the file:\n %s\n",name->fname.PINARRAY_cal_par);
          initialize_PINARRAY_calibration(&cal_par->pinarray,name->fname.PINARRAY_cal_par);
	  
        }
      else
        {
          printf("\nPINARRAY calibration parameters not defined\n");
          exit(-1);
        }



  sort(name);

  if((output=fopen("gated_rings.spn","w"))==NULL)
    {
      printf("\nI can't open file ring_cal_core_energy.spn\n");
      exit(EXIT_FAILURE);
    }
	  
  fwrite(hist,NRING*S4K*sizeof(int),1,output);
	 
  fclose(output);


}

