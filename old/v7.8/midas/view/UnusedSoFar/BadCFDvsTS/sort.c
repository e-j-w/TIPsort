#include "sort.h"

/*================================================================*/
int analyze_fragment(Tig10_event* ptr, short* waveform)
{
  double time;
  if(ptr->channel==chn)
    {
      
      time=ptr->cfd&0x00ffffff;
      time-=(ptr->timestamp*16)&0x00ffffff;
      time*=0.625;
      time/=10.;
      time+=S1K;
      if(time>=low)
	if(time<=high)
	  {
	    printf("=========================================================\n");
	    printf("Trigger  number : %8d 0x%8.8x\n",ptr->trigger_num&0x0fffffff,ptr->trigger_num);
	    printf("           Port : %8d 0x%8.8x\n",ptr->port,ptr->port);
	    printf("          Tig10 : %8d 0x%8.8x\n",ptr->tig10,ptr->tig10);
	    printf("      Collector : %8d 0x%8.8x\n",ptr->collector,ptr->collector);
	    printf("        Channel : %8d 0x%8.8x\n",ptr->channel,ptr->channel);
	    printf("         Charge : %8d 0x%8.8x\n",ptr->charge,ptr->charge);
	    printf("            CFD : %8d 0x%8.8x\n",ptr->cfd,ptr->cfd);
	    printf("            LED : %8d 0x%8.8x\n",ptr->led,ptr->led);
	    printf("     Time Stamp : %8d 0x%8.8x\n",ptr->timestamp,ptr->timestamp);
	    printf("  Time Stamp Up : %8d 0x%8.8x\n",ptr->timestamp_up,ptr->timestamp_up);
	    printf("           Time : %8.0f \n",time, time);
	    printf("Waveform length : %8d 0x%8.8x\n",ptr->waveform_length,ptr->waveform_length);
	    
	    if(ptr->waveform_length!=0)
	      {
		if(h!=NULL) delete h;
		h=new TH1D("Waveform","Waveform",ptr->waveform_length,0,ptr->waveform_length);
		if(c!=NULL) delete c;
		c = new TCanvas("Waveform", "Waveform",10,10, 700, 500);
		
		for(Int_t i=0;i<ptr->waveform_length;i++)
		  h->Fill(i,waveform[i]);
		
		h->Draw();
		theApp->Run(kTRUE);
	      }
	  }
    
    }


  return 0;
}
/*================================================================*/
int main(int argc, char *argv[])
{

  int ac;
  char* av[10];
 if(argc!=5)
    {
      printf("\n ./view_BadCFDvsTS midas_input_data_file_name channel low high\n");
      exit(-1);
    }
  chn=atoi(argv[2]);
  low=atoi(argv[3]);
  high=atoi(argv[4]);
  theApp=new TApplication("App", &ac, av);
/* do sorting */
  sort_but_not_assemble(argv[1]);
 /* display results */
  printf("\n"); 

printf("Program provides information on selected channel\n");
 
}
