#include "sort.h"
/*================================================================*/
Double_t myfunction(Int_t i)
{ 
  Double_t x,s;

  x=i-par->t[0];

  if(x<=0) return par->am[0];
  else
    {   
      s=par->am[0];
      s+=par->am[1]*exp(-x/par->t[1]);
      s-=par->am[2]*exp(-x/par->t[2]);
      s-=par->am[3]*exp(-x/par->t[3]);
      return s;
    }
}

/*================================================================*/
int analyze_fragment(Tig10_event* ptr,short* waveform)
{
  Int_t d;
  WaveFormPar wpar;
  double ch;
  
  if(ptr->channel==chn)
    if((d=ptr->waveform_length)!=0)
      {
	wpar.baseline_range=CSI_BASELINE_RANGE;
	get_t10t90(d,waveform,&wpar);
	get_exclusion_zone_for_CsI(d,waveform,&wpar);
	if(wpar.t10t90_flag==1)
	  {
	    get_shape(4,d,waveform,par,&wpar);
	    get_CsI_t0(d,waveform,par,&wpar);
	    ch=par->chisq/par->ndf;
	    
	    if((ch>=chmin) && (ch<=chmax))
	      {
		print_fragment_info(ptr,S16K);
		printf("   Waveform    d: %10d\n",d);
		printf("   Waveform    b: %10.1f\n",wpar.b);
		printf("   Waveform  max: %10.1f\n",wpar.max);
		printf("   Waveform tmax: %10.4f\n",wpar.tmax);
		printf("  Waveform tlmin: %10d\n",wpar.tlmin);
		printf("  Waveform tlmax: %10d\n",wpar.tlmax);
		printf("   Waveform  t30: %10.4f\n",wpar.t30);
		printf("   Waveform  t10: %10.4f\n",wpar.t10);
		printf("   Waveform  t90: %10.4f\n",wpar.t90);
		printf(" Risetime T10T90: %10.4f\n",wpar.t10t90);
		printf("   Waveform   b0: %10.4f\n",wpar.b0);
		printf("   Waveform   b1: %10.4f\n",wpar.b1);
		printf("   Waveform   s0: %10.4Lf\n",wpar.s0);
		printf("   Waveform   s1: %10.4Lf\n",wpar.s1);
		printf("   Waveform   s2: %10.4Lf\n",wpar.s2);
		printf("   Waveform   t0: %10.4f\n",wpar.t00);
		printf("------ Chisq/NDF: %10.4f\n",ch);		

		for(Int_t i=0;i<5;i++)
		  printf("Id %1d Amp. %10.3f  T %10.3Lf.\n",i,par->am[i],par->t[i]);

		if(h!=NULL) delete h;		
		h=new TH1D("Waveform","Waveform",ptr->waveform_length,0,ptr->waveform_length);

		if(g!=NULL) delete g;
		g=new TH1D("Fit","Fit",ptr->waveform_length,0,ptr->waveform_length);

		if(c!=NULL) delete c;
		c = new TCanvas("Waveform", "Waveform",10,10, 700, 500);

		for(Int_t i=0;i<ptr->waveform_length;i++)
		  {
		    h->Fill(i,waveform[i]);
		    g->Fill(i,myfunction(i));
		  }
		h->SetLineColor(kBlack);
		h->Draw();
		g->SetLineColor(kCyan);
		g->Draw("same");
		
		if(wpar.t00>0)
		  {
		    if(base!=NULL) delete base;
		    base=new TH1D("Baseline","Baseline",ptr->waveform_length,0,ptr->waveform_length);
		    
		    if(fit!=NULL) delete fit;
		    fit=new TH1D("Local fit","Local fit",ptr->waveform_length,0,ptr->waveform_length);

		    if(t0Spot!=NULL) delete t0Spot;
		    t0Spot=new TH1D("Location of t0Spot","Location of t0Spot",ptr->waveform_length,0,ptr->waveform_length);
		    
		    for(Int_t i=0;i<=(int)(wpar.t00+5);i++)
		      base->Fill(i,wpar.b0+wpar.b1*i);
		    
		    for(Int_t i=(int)(wpar.t00-5);i<=wpar.tlmax;i++)
		      fit->Fill(i,wpar.s0+wpar.s1*i+wpar.s2*i*i);

		    //Fill the location of t0, plus 2 bins before and after
		    for(Int_t i=(int)(wpar.t00-2);i<=(int)(wpar.t00+2);i++)
		      t0Spot->Fill(i,waveform[i]);

		    t0Spot->SetFillColor(kGreen);
		    t0Spot->Draw("same");
		    base->SetLineColor(kRed);
		    base->Draw("same");
		    fit->SetLineColor(kBlue);
		    fit->Draw("same");
		  }
		theApp->Run(kTRUE);
	      }}}
  return 0;
}
/*================================================================*/
int main(int argc, char *argv[])
{
  char runFile[132];

  if(argc!=8)
    {
      printf("wfit_CsI_to midas_input_data channel trc tf ts min max\n");
      exit(-1);
    }
  
  par=(ShapePar*)malloc(sizeof(ShapePar));
  memset(par,0,sizeof(ShapePar));
  
  strcpy(runFile,argv[1]);
  chn=atoi(argv[2]);
  par->t[1]=atof(argv[3]);
  par->t[2]=atof(argv[4]);
  par->t[3]=atof(argv[5]);
  chmin=atof(argv[6]);
  chmax=atof(argv[7]);

  theApp=new TApplication("App", &argc, argv);
  sort_but_not_assemble(runFile);
}
