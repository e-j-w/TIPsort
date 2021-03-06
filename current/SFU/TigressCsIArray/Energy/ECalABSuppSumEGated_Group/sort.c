#include "sort.h"

int analyze_data(raw_event *data)
{
  cal_event* cev;
  int suppFlag=0;
  int take=0;
  int i=0;
  int j=0;
  unsigned long long int one=1;
  int csi1,csi2;
  double* energy;
  int* group;

  if((data->h.setupHP&TIGRESS_BIT)==0)
    return SEPARATOR_DISCARD;

  /*if((data->h.setupHP&RF_BIT)==0)
    return SEPARATOR_DISCARD;*/

  if((data->h.setupHP&CsIArray_BIT)==0)
    return SEPARATOR_DISCARD;

  cev=(cal_event*)malloc(sizeof(cal_event));
  memset(cev,0,sizeof(cal_event));
  calibrate_CSIARRAY(data,&cal_par->csiarray,&cev->csiarray);
  calibrate_TIGRESS(data,&cal_par->tg,&cev->tg);
  
  energy=(double*)calloc(cev->tg.h.FA,sizeof(double));
  group=(int*)calloc(cev->tg.h.FA,sizeof(int));
  
  if(cev->tg.h.FA>0) //addback fold>0
    //look through each Tigress position
    for(pos=1;pos<NPOSTIGR;pos++)
      {
      	//reset suppression flag for every position
      	suppFlag=0;
      	//check if the position is in the hit pattern
      	if((cev->tg.h.HHP&(1<<(pos-1)))!=0)
			    //check the position fold
			    if(cev->tg.det[pos].hge.FH>0)
			      //check if the position is in the addback hit pattern
			      if((cev->tg.h.AHP&(1<<(pos-1)))!=0)
			        {
				        //reset take for add-back suppression
				        take=0;
				        //Run through four cores for each position
				        for(col=0;col<NCOL;col++)
				          {
				            //check if the position and color is in the hit pattern
				            if((cev->tg.det[pos].hge.HHP&(1<<col))!=0)
				              //check the fold
				              if(cev->tg.det[pos].ge[col].h.FH>0)
					              {
		              			  //suppress if the position is in the map and has not yet been suppressed
		              			  if(cev->tg.det[pos].ge[col].suppress>=supLow && cev->tg.det[pos].ge[col].suppress<=supHigh && take==0)
		              			    {
		                  			  /* once suppression flag is set
		                  				do not reset it, could remove the take bit
		                  				and keep resetting suppFlag, but this
		                  				is nicer */
		              			      suppFlag=1;
		              			      take=1;
					                    //printf("event at pos %d col %d suppressed\n",pos,col);
		              			    }
					              }
				          }
					   
		            energy[i]=cev->tg.det[pos].addback.E/cal_par->tg.contr_e;
			          colAddBack=cev->tg.det[pos].addbackC;

                for(csi1=1;csi1<NCSI;csi1++)
                  if((cev->csiarray.h.HHP[csi1/64]&(one<<csi1%64))!=0)
                    {
                      if(cal_par->tg.groupflag[pos][colAddBack][csi1][0]==1)
                        {
                          //groups defined by one csi hit
                          group[i] = cal_par->tg.group_map[pos][colAddBack][csi1][0]+NGROUP*suppFlag;
                          //printf("one hit - group %i\n",group[i]);
                        }
                      else
                        {
                          //groups defined by two csi hits
                          for(csi2=csi1+1;csi2<NCSI;csi2++)
                            if((cev->csiarray.h.HHP[csi2/64]&(one<<csi2%64))!=0)
                              if(cal_par->tg.groupflag[pos][colAddBack][csi1][csi2]==2)
                                {
                                  group[i] = cal_par->tg.group_map[pos][colAddBack][csi1][csi2]+NGROUP*suppFlag;
                                  //printf("two hits - group %i\n",group[i]);
                                }
                        }
                    }
			          
			          if(energy[i]>=0 && energy[i]<S32K)
                	projhist[group[i]][(int)(energy[i])]++;
			          
			          i++;//increment event counter
			          
			        }
      }
  
  if(i!=cev->tg.h.FA)
    printf("WARNING: Addback fold not equal to the number of events seen!\n");
  
  for(i=0;i<cev->tg.h.FA;i++)
    {
      //look for a gamma that falls into the gate
      if(energy[i]>=0)
	      if(energy[i]<S32K)
	        if(energy[i]>=gateELow/cal_par->tg.contr_e)
		        if(energy[i]<=gateEHigh/cal_par->tg.contr_e)
		          if(group[i]>0)
		            if(group[i]<NGROUP)
		              {
		                //add all gammas in the event that aren't the gamma that fell into the gate
		                for(j=0;j<cev->tg.h.FA;j++)
		                	{
				                if(j!=i)
				                  if(energy[j]>=0)
			                      if(energy[j]<S32K)
				                      if(group[j]>0)
				                        if(group[j]<NGROUP)
				                          hist[group[j]][(int)(energy[j])]++;
		                  }
		                gatehist[group[i]][(int)(energy[i])]++;
		                break;//don't double count
		              }
    }
  
  //printf("***** END OF EVENT *****\n");
  //getc(stdin);
  
  free(cev);
  free(energy);  
  free(group);
  return SEPARATOR_DISCARD;
}
/*====================================================================================*/
int main(int argc, char *argv[])
{
  FILE * output;
  input_names_type* name;
  FILE *cluster;
  char n[132],FileName[132];
  double avgGateE;
  double gateWidth;

  if(argc!=6)
    {
      printf("TigressCsIArray_ECalABSuppSumEGated_Group master_file_name supLow supHigh gateELow gateEHigh\n");
      printf("Program attempts to generate energy gated gamma ray spectra (not background subtracted).  Energy gate values (gateEHigh, gateELow) should be specified in keV.  Sorted spectrum will be mapped by group rather than Tigress ring.\n");
      exit(-1);
    }
  
  printf("Program attempts to generate energy gated gamma ray spectra.\n");
  
  name=(input_names_type*)malloc(sizeof(input_names_type));
  memset(name,0,sizeof(input_names_type));
  
  cal_par=(calibration_parameters*)malloc(sizeof(calibration_parameters));
  memset(cal_par,0,sizeof(calibration_parameters));
  
  memset(hist,0,sizeof(hist));
  memset(gatehist,0,sizeof(gatehist));
  memset(projhist,0,sizeof(projhist));
  read_master(argv[1],name);
  
  supLow=atof(argv[2]);
  supHigh=atof(argv[3]);
  gateELow = atof(argv[4]);
  gateEHigh = atof(argv[5]);

  if(gateELow>gateEHigh)
  	{
  		printf("ERROR: cannot have lower gate bound larger than upper gate bound.\n");
  		exit(-1);
  	}
  	
  avgGateE=(gateELow+gateEHigh)/2.;
  gateWidth=gateEHigh-gateELow;

  if(name->flag.cluster_file==1)
    {
      printf("\nSorting data from cluster file:\n %s\n",name->fname.cluster_file);
      if((cluster=fopen(name->fname.cluster_file,"r"))==NULL)
				{
					printf("\nI can't open input file %s\n",name->fname.cluster_file);
					exit(-2);
				}
    }
  else
    {
      printf("\nCluster file not defined\n");
      exit(-1);
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
  if(name->flag.CSIARRAY_cal_par==1)
    {
      printf("CsIArray calpar read from: %s\n",name->fname.CSIARRAY_cal_par);
      initialize_CSIARRAY_calibration(&cal_par->csiarray,name->fname.CSIARRAY_cal_par);
    }
  else
    {
      printf("\nCSIARRAY calibration parameters not defined\n");
      exit(EXIT_FAILURE);
    }
  
  printf("\n");
  printf("  Using gamma energy gate between %f and %f keV.\n",gateELow,gateEHigh);
  printf("\n");
  
  name->flag.inp_data=1; 
  while(fscanf(cluster,"%s",n)!=EOF)
    {
      strcpy(name->fname.inp_data,n);
      printf("Sorting data from file %s\n", name->fname.inp_data);
      sort(name);
    }

  fclose(cluster);
  
  printf("\n");

  sprintf(FileName,"ECalABSuppGroupSum_c%.0f_w%.0fgated.mca",avgGateE,gateWidth);
  if((output=fopen(FileName,"w"))==NULL)
    {
      printf("ERROR!!! I cannot open the mca file!\n");
      exit(EXIT_FAILURE);
    }
  fwrite(hist,2*NGROUP*S32K*sizeof(int),1,output);
  fclose(output);
  printf("Gated spectrum saved to file: %s\n",FileName);
  
  sprintf(FileName,"ECalABSuppGroupSum_c%.0f_w%.0fgate.mca",avgGateE,gateWidth);
  if((output=fopen(FileName,"w"))==NULL)
    {
      printf("ERROR!!! I cannot open the mca file!\n");
      exit(EXIT_FAILURE);
    }
  fwrite(gatehist,2*NGROUP*S32K*sizeof(int),1,output);
  fclose(output);
  printf("Gate spectrum saved to file: %s\n",FileName);
  
  sprintf(FileName,"ECalABSuppGroupSum_c%.0f_w%.0fproj.mca",avgGateE,gateWidth);
  if((output=fopen(FileName,"w"))==NULL)
    {
      printf("ERROR!!! I cannot open the mca file!\n");
      exit(EXIT_FAILURE);
    }
  fwrite(projhist,2*NGROUP*S32K*sizeof(int),1,output);
  fclose(output);
  printf("Projection spectrum saved to file: %s\n",FileName);

}

  

