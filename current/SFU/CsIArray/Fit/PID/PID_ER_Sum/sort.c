#include "sort.h"

int analyze_data(raw_event *data)
{
  uint64_t one=1;
  int type;
  double chi,s,f,r,e;

  if((data->h.setupHP&RF_BIT)==0) return SEPARATOR_DISCARD;
  if((data->h.setupHP&CsIArray_BIT)==0) return SEPARATOR_DISCARD;

  for(pos=1;pos<NCSI;pos++)
    if((data->csiarray.h.THP[pos/64]&(one<<pos%64))!=0)
      {
	type=data->csiarray.wfit[pos].type;
	chi=data->csiarray.wfit[pos].chisq;
	chi/=data->csiarray.wfit[pos].ndf;
	
	if(type>=idmin && type<=idmax)
	  if(chi>=chimin && chi<=chimax)
	    {
	      e=data->csiarray.wfit[pos].am[1];
	      s=data->csiarray.wfit[pos].am[3];
	      f=data->csiarray.wfit[pos].am[2];
	      if(f==0)
		r=s*100;
	      else
		r=s/f*100;
	      r+=100;

	      if(r>S32K-4) r=S32K-4;
	      h[pos]->Fill(e,r);
	      h_e[pos]->Fill(e);
	      	      
	      hist1[(int)rint(r)]++;
	      hist2[(int)rint(e)]++;
	    }}
  
  return SEPARATOR_DISCARD;
}
/*=================================================================*/
int main(int argc, char *argv[])
{
  FILE  *cluster,*output;
  input_names_type* name;
  char DataFile[132];
  char HistName1[132],HistName2[132];

  if(argc!=6)
    {
      printf("CsIArray_PID_ER_Sum cluster_file idmin idmax chimin chimax\n");
      exit(-1);
    }

  for(pos=1;pos<NCSI;pos++)
    {
      sprintf(HistName1,"PID_ER_%02d",pos);
      h[pos] = new TH2D(HistName1,HistName1,S8K,0,S32K,S1K/2,0,S1K/2);
      h[pos]->Reset();
      sprintf(HistName2,"PID_E_%02d",pos);
      h_e[pos] = new TH1D(HistName2,HistName2,S8K,0,S32K);
      h_e[pos]->Reset();
    }

  idmin=atoi(argv[2]);
  idmax=atoi(argv[3]);
  chimin=atof(argv[4]);
  chimax=atof(argv[5]);
  
  printf("Program sorts CsIArray_PID_ER_Sum.\n");
  
  memset(hist1,0,sizeof(hist1));
  memset(hist2,0,sizeof(hist2));

  name=(input_names_type*)malloc(sizeof(input_names_type));
  memset(name,0,sizeof(input_names_type));
  
  if((cluster=fopen(argv[1],"r"))==NULL)
    {
      printf("ERROR!!! I cannot open the cluster file!\n");
      exit(-2);
    }
  
  while(fscanf(cluster,"%s",DataFile)!=EOF)
    {
      printf("Reading data from file %s.\n",DataFile);
      strcpy(name->fname.inp_data,DataFile);
      sort(name);
    }
  
  fclose(cluster);

  if((output=fopen("CSIArray_PID_ratio.mca","w"))==NULL)
    {
      printf("ERROR!!! I cannot open the mca file!\n");
      exit(EXIT_FAILURE);
    }
  //for(int pos=0;pos<NCSI;pos++)
  fwrite(hist1,S32K*sizeof(int),1,output);
  fclose(output);

  if((output=fopen("CSIArray_PID_energy.mca","w"))==NULL)
    {
      printf("ERROR!!! I cannot open the mca file!\n");
      exit(EXIT_FAILURE);
    }
  //for(int pos=0;pos<NCSI;pos++)
  fwrite(hist2,S32K*sizeof(int),1,output);
  fclose(output);

  TFile *f= new TFile("PID.root","RECREATE");
  f->cd();
  for(pos=1;pos<NCSI;pos++)
    {
      sprintf(HistName1, "PID_ER_%02d", pos);
      sprintf(HistName2, "PID_E_%02d", pos);
      h[pos]->Draw("COLZ");
      h[pos]->Write(HistName1);
      h_e[pos]->Draw("");
      h_e[pos]->Write(HistName2);
    }
  f->Close();
  
}
