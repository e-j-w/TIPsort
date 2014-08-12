#include "sort.h"

/*================================================================*/
int analyze_fragment(Tig10_event* ptr, short* waveform)
{
  int time,i;
  int length;
  long sup,sdown,s;
  double v;
  int   TLE,TCFD;

  if(ptr->channel==chn)
    {
      time=ptr->cfd&0x00ffffff;
      time-=(ptr->timestamp*16)&0x00ffffff;
      printf("=========================================================\n");
      printf(" Trigger  number : %8d 0x%8.8x\n",ptr->trigger_num&0x0fffffff,ptr->trigger_num);
      printf("            Port : %8d 0x%8.8x\n",ptr->port,ptr->port);
      printf("           Tig10 : %8d 0x%8.8x\n",ptr->tig10,ptr->tig10);
      printf("       Collector : %8d 0x%8.8x\n",ptr->collector,ptr->collector);
      printf("         Channel : %8d 0x%8.8x\n",ptr->channel,ptr->channel);
      printf("          Charge : %8d 0x%8.8x\n",ptr->charge,ptr->charge);
      printf("             CFD : %8d 0x%8.8x\n",ptr->cfd,ptr->cfd);
      printf("      Time Stamp : %8d 0x%8.8x\n",ptr->timestamp,ptr->timestamp);
      printf("   Time Stamp Up : %8d 0x%8.8x\n",ptr->timestamp_up,ptr->timestamp_up);
      printf("            Time : %8d 0x%8.8x\n",time+S16K, time+S16K);
      printf(" Waveform length : %8d 0x%8.8x\n",ptr->waveform_length,ptr->waveform_length);

      if(ptr->waveform_length!=0)
	{
	  if(h!=NULL) delete h;
	  if(d!=NULL) delete d;
	  if(cfd!=NULL) delete cfd;
	  length=ptr->waveform_length;
	  h=new TH1D("Waveform","Waveform",length,0,length);
	  d=new TH1D("Diff","Diff",length,0,length);
	  cfd=new TH1D("CFD","CFD",length,0,length);	  
	  if(c!=NULL) delete c;
	  c = new TCanvas("Waveform", "Waveform",10,10, 700, 500);
	  c->Divide(1,3);
          sdown=0;
	  sup=0;
	  h->Reset();
	  d->Reset();
	  cfd->Reset();
	  for(i=0;i<window;i++)
	    {
	      sdown+=waveform[i];
	      h->Fill(i,waveform[i]);
	    }
	  for(i=window;i<window+gap;i++)
	    {
	      h->Fill(i,waveform[i]);
	    }

	  for(i=window+gap;i<2*window+gap;i++)
	    {
	      sup+=waveform[i];
	      h->Fill(i,waveform[i]);
	    }

	  for(i=2*window+gap;i<ptr->waveform_length;i++)
	    {
	      sup+=waveform[i];
	      h->Fill(i,waveform[i]);
	      sup-=waveform[i-window];

	      sdown+=waveform[i-window-gap];
	      sdown-=waveform[i-2*window-gap];
	      s=sup-sdown;
	      d->Fill(i,(int)s);


	    }
	  for(i=2*window+gap;i<length;i++)
	    {
	      v=-fraction*d->GetBinContent(i)+d->GetBinContent(i-delay);
	      cfd->Fill(i,v);
	    }

	  TLE=0;
	  while(d->GetBinContent(TLE)<hit_threshold)
	    {
	      TLE++;
	      if(TLE>ptr->waveform_length)
		{
		  TLE=-1;
		  break;
		}
	    }
	  
	  TCFD=-1;
	  if(TLE>0&&TLE<ptr->waveform_length)
	    {
	      TCFD=TLE;
		while(1)
		  {
		    TCFD++;
		    if(TCFD>ptr->waveform_length)
		      {
			TCFD=-1;
			break;
		      }
		    else
		      {
			if(cfd->GetBinContent(TCFD)<=0)
			  if(cfd->GetBinContent(TCFD+1)>0)
			    break;
		      }
		  }
	    }

	  printf("Leading Edge Time: %8d 0x%8.8x\n",TLE, TLE);
	  printf("         CFD Time: %8d 0x%8.8x\n",16*TCFD, 16*TCFD);
	  c->cd(1);
	  h->Draw();
	  c->cd(2);
	  d->Draw();
	  c->cd(3);
	  cfd->Draw();
	  theApp->Run(kTRUE);
	}
    }


  return 0;
}
/*================================================================*/
int main(int argc, char *argv[])
{

  int ac;
  char* av[10];
 if(argc!=7)
    {
      printf("\n ./view_CFD midas_input_data_file_name channel clip_delay hit_threshold fraction delay\n");
      exit(-1);
    }

 window=1;
 chn=atoi(argv[2]);
 gap=atof(argv[3]);
 hit_threshold=atoi(argv[4]);
 fraction=atof(argv[5]);
 delay=atoi(argv[6]);
 theApp=new TApplication("App", &ac, av);
 /* do sorting */
 sort_but_not_assemble(argv[1]);
 /* display results */
 printf("\n"); 
 
 printf("Program provides information on selected channel\n");
 
}
