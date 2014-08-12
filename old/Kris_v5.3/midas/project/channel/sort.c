#include "sort.h"

/*================================================================*/
int analyze_fragment(Tig10_event* ptr, short* waveform)
{
  
  if(ptr->channel==channel)
    if(ptr->charge_flag!=0)
      h->Fill(ptr->charge/K);

  return 0;
}
/*================================================================*/
int main(int argc, char *argv[])
{
  TCanvas *c;
  TApplication *theApp;
  int ac;
  char *av[10];
  char name[132];

  if(argc!=4)
    {
      printf("\n ./project_channel midas_input_data_file_name channel K\n");
      exit(-1);
    }
  channel=atoi(argv[2]);
  K=atof(argv[3]);
  h=new TH1D("","",S65K,0,S65K);
  sprintf(name,"ch. %s",argv[2]);
  h->Reset();
  h->SetName(name);
  h->SetTitle(name);

/* do sorting */
  sort_but_not_assemble(argv[1]);

 /* display results */
  printf("\n"); 
  printf("Program provides a 1D projection for a selected DAQ channel\n");

  theApp=new TApplication("App", &ac, av);
  c = new TCanvas(name, name,10,10, 700, 500);
  h->GetXaxis()->SetTitle("Charge");
  h->GetYaxis()->SetTitle("Counts");
  h->Draw();
  theApp->Run(kTRUE);


}
