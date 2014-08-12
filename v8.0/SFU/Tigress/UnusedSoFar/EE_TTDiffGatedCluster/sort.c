#include "sort.h"

int analyze_data(raw_event *data)
{
  cal_event* cev;
  int pos1,pos2;
  double t1,t2,t;
  int  e1,e2;
  
  cev=(cal_event*)malloc(sizeof(cal_event));
  memset(cev,0,sizeof(cal_event));
  calibrate_TIGRESS(data,&cal_par->tg,&cev->tg);

  if(cev->tg.h.FA>1)
    for(pos1=1;pos1<NPOSTIGR;pos1++)
      if((cev->tg.h.AHP&(1<<(pos1-1)))!=0)
	{
	  t1=cev->tg.det[pos1].addback.T;
	  e1=(int)rint(cev->tg.det[pos1].addback.E/cal_par->tg.contr_e);
	  for(pos2=pos1+1;pos2<NPOSTIGR;pos2++)
	    if((cev->tg.h.AHP&(1<<(pos2-1)))!=0)
	      {
		t2=cev->tg.det[pos2].addback.T;
		e2=(int)rint(cev->tg.det[pos2].addback.E/cal_par->tg.contr_e);
		t=t1-t2+S2K;
		if(t>=tlow)
		  if(t<=thigh)
		    if(e1>=0)
		      if(e1<S4K)
			if(e2>=0)
			  if(e2<S4K)
			    {
			      mat[e1][e2]++;
			      mat[e2][e1]++;
			    }

	      }
	}


  free(cev);
  return SEPARATOR_DISCARD;
}
/*====================================================================================*/
int main(int argc, char *argv[])
{
  FILE * output;
  input_names_type* name;
  FILE *cluster;
  char n[132];



  if(argc!=4)
    {
      printf("\n ./TIGRESS_EE_TTDiffGatedCluster master_file_name tlow thigh\n");
      exit(-1);
    }
  
  printf("Program sorts calibrated time histogram for TIGRESS cores \n");
  name=(input_names_type*)malloc(sizeof(input_names_type));
  memset(name,0,sizeof(input_names_type));
  cal_par=(calibration_parameters*)malloc(sizeof(calibration_parameters));
  memset(cal_par,0,sizeof(calibration_parameters));
  memset(hist,0,sizeof(hist));
  memset(mat,0,sizeof(mat));
  read_master(argv[1],name);
  tlow=atof(argv[2]);
  thigh=atof(argv[3]);

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

  name->flag.inp_data=1; 
  while(fscanf(cluster,"%s",n)!=EOF)
    {
      strcpy(name->fname.inp_data,n);
      sort(name);
    }

  fclose(cluster);


 
  if((output=fopen("matrix.mat","w"))==NULL)
    {
      printf("\nI can't open file matrix.mat\n");
      exit(EXIT_FAILURE);
    }
	  
  fwrite(mat,S4K*S4K*sizeof(short int),1,output);
  fclose(output);

  for(int i=0;i<S4K;i++)
    for(int j=0;j<S4K;j++)
      hist[i]+=mat[i][j];
    

  if((output=fopen("projection.spn","w"))==NULL)
    {
      printf("\nI can't open file projection.spn\n");
      exit(EXIT_FAILURE);
    }
	  
  fwrite(hist,S4K*sizeof(int),1,output);
  fclose(output);
  
}

  

