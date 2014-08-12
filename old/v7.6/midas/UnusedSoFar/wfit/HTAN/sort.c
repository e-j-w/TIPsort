#include "sort.h"

/*================================================================*/
Double_t myfunction(Double_t *y, Double_t *par)
{
  Double_t x=(y[0]-par[3])/par[2];

  Double_t f=par[0]+par[1]*((exp(x)-exp(-x))/(exp(x)+exp(-x))+1);

  return f;
}

/*================================================================*/
int analyze_fragment(Tig10_event* ptr,short* waveform)
{
  double b;
 if(ptr->channel==chn)
  if(ptr->waveform_length!=0)
    {
     
      if(h!=NULL) delete h;
      h=new TH1D("Waveform","Waveform",ptr->waveform_length,0,ptr->waveform_length);

      h->Reset();
      for(Int_t i=0;i<ptr->waveform_length;i++)
	{
	  h->Fill(i,waveform[i]);
	}
      b=0.;
      for(Int_t i=0;i<base;i++)
	{
	  b+=waveform[i];
	}      
      b/=base;
    
      f->FixParameter(0,b);
      f->SetParameter(1,900);
      f->SetParameter(2,150);
      f->SetParameter(3,30);

  
      h->Fit("wfit","","",0,(Double_t)ptr->waveform_length-1);
      h->Draw();
  
      theApp->Run(kTRUE);
     
    }

  return 0;
}
/*================================================================*/
int main(int argc, char *argv[])
{

  int ac;
  char *av[10];
  if(argc!=4)
    {
      printf("\n ./wfit_HTAN midas_input_data_file_name channel baseline_range\n");
      exit(-1);
    }
  chn=atoi(argv[2]);
  base=atoi(argv[3]);
  theApp=new TApplication("App", &ac, av);
  f=new TF1("wfit",myfunction,0,MAX_SAMPLES,4);
  f->SetParName(0,"C");
  f->SetParName(1,"A");
  f->SetParName(2,"tau");
  f->SetParName(3,"t0");
  f->SetLineColor(kRed);
  f->SetLineWidth(0.3);
 
/* do sorting */
  sort_but_not_assemble(argv[1]);
 /* display results */
  printf("\n"); 
  printf("Program shows waveform fit for selected channel\n");


}
