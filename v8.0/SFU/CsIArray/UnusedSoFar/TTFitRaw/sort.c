#include "sort.h"

int analyze_data(raw_event *data)
{
  unsigned long long int one=1;
  int pos1,pos2,type;
  double t1,t2,chi;


  for(pos1=1;pos1<NCSI;pos1++)
    if((data->csiarray.h.THP&(one<<pos1))!=0)
      {
	type=data->csiarray.wfit[pos1].type;
	chi=data->csiarray.wfit[pos1].chisq;
	chi/=data->csiarray.wfit[pos1].ndf;
	if(type>=idmin)
	  if(type<=idmax)
	    if(chi>=chimin)
	      if(chi<=chimax)
		{
		  t1=16.*data->csiarray.wfit[pos1].t[0]; 
		
		  if((data->h.setupHP&RF_BIT)!=0)
		    {
		      t1-=data->rf.sin.t0;
		   
		      t1+=S4K;
		      t1/=16;
		      if(t1<0)t1=S32K-2;
		      if(t1>S32K) t1=S32K-3;
		    }
		  else
		    t1=S32K-4;

		  for(pos2=pos1+1;pos2<NCSI;pos2++)
		    if((data->csiarray.h.THP&(one<<pos2))!=0)
		      {
			type=data->csiarray.wfit[pos2].type;
			chi=data->csiarray.wfit[pos2].chisq;
			chi/=data->csiarray.wfit[pos2].ndf;
			if(type>=idmin)
			  if(type<=idmax)
			    if(chi>=chimin)
			      if(chi<=chimax)
				{
				  t2=16.*data->csiarray.wfit[pos2].t[0]; 
				  if((data->h.setupHP&RF_BIT)!=0)
				    {
				      t2-=data->rf.sin.t0;		   
				      t2+=S4K;
				      t2/=16;
				      if(t2<0)t2=S32K-2;
				      if(t2>S32K) t2=S32K-3;
				    }
				  else
				    t2=S32K-4;
				  if(t1==t2)
				    h->Fill(t1,t2);
				  else
				    {
				      h->Fill(t1,t2);
				      h->Fill(t2,t1);
				    }
				}
		      }
		  
		}
      }
  return SEPARATOR_DISCARD;
}
/*====================================================================================*/
int main(int argc, char *argv[])
{

  input_names_type* name;
  TCanvas *canvas;
  TApplication *theApp;


  if(argc!=6)
    {
      printf("\n ./CSIARRAY_TTFitRaw sfu_input_data_file_name idmin idmax chimin chimax\n");
      exit(-1);
    }
  h = new TH2D("TTFit time","TTFit time",S1K,0,S1K,S1K,0,S1K);
  h->Reset();
  idmin=atoi(argv[2]);
  idmax=atoi(argv[3]);
  chimin=atof(argv[4]);
  chimax=atof(argv[5]);
  printf("Program sorts raw time/time matrix from the CSIARRAY fits \n");

  name=(input_names_type*)malloc(sizeof(input_names_type));
  memset(name,0,sizeof(input_names_type));
  strcpy(name->fname.inp_data,argv[1]);
  sort(name);

  
 theApp=new TApplication("App", &argc, argv);
  canvas = new TCanvas("TTFit time", "TTFit time",10,10, 500, 500);
  gPad->SetLogz(1);
  gStyle->SetPalette(1);
  h->Draw("COLZ");
  
  theApp->Run(kTRUE);

}

