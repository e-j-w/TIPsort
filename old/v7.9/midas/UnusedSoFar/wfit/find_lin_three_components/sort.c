#include "sort.h"
/*================================================================*/
Double_t myfunction(Int_t i)
{
 
  Double_t x,s;


  x=i-par->t[0];
  if(x<=0) 
    return par->am[0];
  else
    {   
      s=par->am[0];
      s+=par->am[1]*exp(-x/par->t[1]);
      for(int i=2;i<NSHAPE;i++)
	s-=par->am[i]*exp(-x/par->t[i]);
   
      return s;
    }
}
/*================================================================*/
int analyze_fragment(Tig10_event* ptr,short* waveform)
{
 

  if(ptr->channel==chn)
    if((ptr->trigger_num&0x00ffffff)==trig_num)
      if(ptr->waveform_length!=0)
	{	  
	  print_fragment_info(ptr,S16K);
	  if(h!=NULL) delete h;		
	  h=new TH1D(label,label,ptr->waveform_length,0,ptr->waveform_length);
	  
	  if(g!=NULL) delete g;
	  g=new TH1D("Fit","Fit",ptr->waveform_length,0,ptr->waveform_length);
	  if(c!=NULL) delete c;
	  c = new TCanvas("Waveform", "Waveform",10,10, 700, 500);
	  for(Int_t i=0;i<ptr->waveform_length;i++)
	    {
	      h->Fill(i,waveform[i]);
	      g->Fill(i,myfunction(i));
	    }
	  h->GetXaxis()->SetTitle("Time [10ns]");
	  h->GetYaxis()->SetTitle("Amplitude [arb.]");
	  h->GetYaxis()->SetRangeUser(rlow,rhigh);
	  h->GetYaxis()->SetTitleOffset(1.25);
	  h->SetStats(0);
	  h->Draw();
	  g->SetLineColor(kRed);
	  g->Draw("same");
	  theApp->Run(kTRUE);
      }

  return 0;
}
/*================================================================*/
int main(int argc, char *argv[])
{

  
 
  int ac;
  char *av[10];

 if(argc!=17)
    {
      printf("\n ./wfit_find_lin_three_components\n midas_input_data_file_name channel trig_number\n t0 trc tf ts tr\n a0 arc af as ar label\n");
      exit(-1);
    }
  par=(ShapePar*)malloc(sizeof(ShapePar));
  memset(par,0,sizeof(ShapePar));
  chn=atoi(argv[2]);
  trig_num=atoi(argv[3]);
  par->t[0]=atof(argv[4]);
  par->t[1]=atof(argv[5]);
  par->t[2]=atof(argv[6]);
  par->t[3]=atof(argv[7]);
  par->t[4]=atof(argv[8]);
  par->am[0]=atof(argv[9]);
  par->am[1]=atof(argv[10]);
  par->am[2]=atof(argv[11]);
  par->am[3]=atof(argv[12]);
  par->am[4]=atof(argv[13]);
  label=argv[14];
  rlow=atof(argv[15]);
  rhigh=atof(argv[16]);
  theApp=new TApplication("App", &ac, av);

/* do sorting */
  sort_but_not_assemble(argv[1]);
 /* display results */

 


}
