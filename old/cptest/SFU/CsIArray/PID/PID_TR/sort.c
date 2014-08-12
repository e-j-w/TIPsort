#include "sort.h"

int analyze_data(raw_event *data)
{
  unsigned long long int one=1;
  int pos,type;
  double chi,s,f,r,t;
  
  pos=det_id;
  if((data->csiarray.h.THP&(one<<pos))!=0)
    //if((data->csiarray.h.TSfold)==1)
    {
      type=data->csiarray.wfit[pos].type;
      chi=data->csiarray.wfit[pos].chisq;
      chi/=data->csiarray.wfit[pos].ndf;
      
      if((type>=idmin) && (type<=idmax))
	if((chi>=chimin) && (chi<=chimax))
	  {
	    t=16.*data->csiarray.wfit[pos].t[0]; 
	    //t=data->csiarray.csi[pos].cfd&0x00ffffff;
	    //t-=(data->csiarray.csi[pos].timestamp*16)&0x00ffffff;
	    
	    if((data->h.setupHP&RF_BIT)!=0)
	      {
		t-=data->rf.sin.t0;
		t+=S8K;
		
		if(t<0)t=S32K-2;
		if(t>S32K) t=S32K-3;
	      }
	    else
	      t=S32K-4;
	    
	    t/=16;
	    s=data->csiarray.wfit[pos].am[3];
	    f=data->csiarray.wfit[pos].am[2];
	    r=s/f*100;
	    
	    h->Fill(t,100+r);
	    if(r>S32K-4) r=S32K-4;
	    
	    hist1[(int)rint(r)]++;
	    hist2[(int)rint(t)]++;
	  }}
  
  return SEPARATOR_DISCARD;
}
/*====================================================================================*/
int main(int argc, char *argv[])
{
  FILE * output;
  input_names_type* name;
  TCanvas *canvas;
  TApplication *theApp;
  
  if(argc!=7)
    {
      printf("CsIArray_PID_TR sfu_input_file idmin idmax chimin chimax det_id\n");
      exit(-1);
    }
  h = new TH2D("PID","PID",S1K,0,S1K,S1K/2,0,S1K/2);
  h->Reset();
  
  idmin=atoi(argv[2]);
  idmax=atoi(argv[3]);
  chimin=atof(argv[4]);
  chimax=atof(argv[5]);
  det_id=atoi(argv[6]);
  
  printf("Program sorts PID histogram for the CsIArray \n");

  memset(hist1,0,sizeof(hist1));
  memset(hist2,0,sizeof(hist2));
  name=(input_names_type*)malloc(sizeof(input_names_type));
  memset(name,0,sizeof(input_names_type));
  strcpy(name->fname.inp_data,argv[1]);

  sort(name);
  
  if((output=fopen("CsIArray_PID_Ratio.mca","w"))==NULL)
    {
      printf("ERROR! I can't open the mca file!\n"); 
      exit(EXIT_FAILURE);
    }
 
  fwrite(hist1,S32K*sizeof(int),1,output);
  fclose(output);

  if((output=fopen("CsIArray_PID_Time.mca","w"))==NULL)
    {
      printf("ERROR! I can't open the mca file!\n"); 
      exit(EXIT_FAILURE);
    }
 
  fwrite(hist2,S32K*sizeof(int),1,output);
  fclose(output);

  theApp=new TApplication("App", &argc, argv);
  canvas = new TCanvas("PID_TR", "PID_TR",10,10, 500, 300);
  gPad->SetLogz(1);
  gStyle->SetPalette(1);
  h->Draw("COLZ");
  
  theApp->Run(kTRUE);
}
