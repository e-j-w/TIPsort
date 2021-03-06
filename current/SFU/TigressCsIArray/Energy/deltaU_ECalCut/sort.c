#include "sort.h"

int analyze_data(raw_event *data)
{
  cal_event* cev;
  double eAddBack=0;
  int suppFlag=0;
  int take=0;
  double u=-1.;
  
  cev=(cal_event*)malloc(sizeof(cal_event));
  memset(cev,0,sizeof(cal_event));

  //calibrate Tigress here - earlier is a waste of CPU
  calibrate_TIGRESS(data,&cal_par->tg,&cev->tg);
  
  //check the Ge fold
  if(cev->tg.h.FH>0)
    //look through each Tigress position
    for(pos=1;pos<NPOSTIGR;pos++)
      {
	//reset suppression flag for each position
	suppFlag=0;
	//check if the position is in the hit pattern
	if((cev->tg.h.HHP&(1<<(pos-1)))!=0) 
	  //check the energy fold
	  if(cev->tg.det[pos].hge.FH>0)
	    //check if the position is in the addback hit pattern
	    if((cev->tg.h.AHP&(1<<(pos-1)))!=0)
	      {
		/* set take=0 for each position for suppression 
		   of all possible clover hits in add-back mode*/
		take=0;
		//Run through four cores for each position for suppression
		for(col=0;col<NCOL;col++)
		  {
		    //Check if this color is indicated in the hit pattern
		    if((cev->tg.det[pos].hge.HHP&(1<<col))!=0)
		      //Check that this combination has a fold great than zero
		      if(cev->tg.det[pos].ge[col].h.FH>0)
			//Check if this combination is indicated in the hit pattern
			if((cev->tg.det[pos].ge[col].h.HHP&1)!=0)
			  //suppress if the position is in the map and has not yet been suppressed
			  if(cev->tg.det[pos].ge[col].suppress>=940 && cev->tg.det[pos].ge[col].suppress<=1060 && take==0)
			    { 
			      /* once suppression flag is set
				 do not reset it, could remove the take bit
				 and keep resetting suppFlag, but this
				 is nicer */   
			      suppFlag=1;
			      take=1;
			    }
		  }
		
		eAddBack = cev->tg.det[pos].addback.E/cal_par->tg.contr_e;  
		colAddBack = cev->tg.det[pos].addbackC;
		ring = cev->tg.det[pos].ge[colAddBack].ring;
		
		if(eAddBack>=elow && eAddBack<ehigh) 
		  {
		    calibrate_CSIARRAY(data,&cal_par->csiarray,&cev->csiarray);
		    u=cev->csiarray.U;
		    //printf("CsIArray ene = %f\n",u);
		    //getc(stdin);
		    if(u>=0 && u<S32K)
		      hist[(int)rint(u)]++;
		    else 
		      hist[S32K-1000]++;
		  }
		//else hist[S32K-1000]++;
	      }
      }
  
  
  free(cev);
  return SEPARATOR_DISCARD;
}
/*=========================================================================*/
int main(int argc, char *argv[])
{
  FILE *output, *cluster;
  input_names_type* name;
  char DataFile[132];

  elow=0.;
  ehigh=0.;  

  if(argc!=4)
    {
      printf("TigressCsI_deltaU_ECalCut master_file_name E_gamma_low E_gamma_high\n");
      exit(-1);
    }
  
  printf("Program sorts deltaU histograms for the CsI array with gamma ray energy cuts imposed by TIGRESS.\n");
  
  name=(input_names_type*)malloc(sizeof(input_names_type));
  memset(name,0,sizeof(input_names_type));
  
  cal_par=(calibration_parameters*)malloc(sizeof(calibration_parameters));
  memset(cal_par,0,sizeof(calibration_parameters));
  
  memset(hist,0,sizeof(hist));
  read_master(argv[1],name);

  elow=atof(argv[2]);
  ehigh=atof(argv[3]);

  //printf("elow: %f ehigh %f\n",elow,ehigh);
  //getc(stdin);
  
  if(name->flag.cluster_file==1)
    {
      printf("Sorting deltaU histograms for the CsI Array based upon the cluster file: %s\n",name->fname.cluster_file);
      if((cluster=fopen(name->fname.cluster_file,"r"))==NULL)
	{
	  printf("ERROR!!! I can't open input file %s\n",name->fname.cluster_file);
	  exit(-2);
	}
    }
  else
    {
      printf("ERROR!!! Cluster file not defined\n");
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

  if(name->flag.CSIARRAY_cal_par==1)
    {
      printf("CsIArray calpar read from: %s\n",name->fname.CSIARRAY_cal_par);
      initialize_CSIARRAY_calibration(&cal_par->csiarray,name->fname.CSIARRAY_cal_par);
    }
  else
    {
      printf("ERROR! CsIArray calibration parameters not defined!\n");
      exit(EXIT_FAILURE);
    }
  
  while(fscanf(cluster,"%s",DataFile) != EOF)
    {
      memset(name,0,sizeof(input_names_type));
      strcpy(name->fname.inp_data,DataFile);
      
      printf("Sorting data from file %s\n", name);
      sort(name);
    }
  
  if((output=fopen("deltaU_ECalCut.mca","w"))==NULL)
    {
      printf("ERROR!!! I cannot open the mca file!\n");
      exit(EXIT_FAILURE);
    }
  fwrite(hist,S32K*sizeof(int),1,output);
  fclose(output);
}
