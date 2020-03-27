/* ECalABSuppSum */

#include "sort.h"

int analyze_data(raw_event *data)
{
  cal_event* cev;
  double eAddBack=-1.;
  double phaseRF=-1.;
  int suppFlag=0;
  int take=0;

  memcpy(&phaseRF,&data->rf.sin.t0,sizeof(double)); //get the RF phase for this raw event (in CFD units)
  phaseRF *= 0.625; //convert to ns

  cev=(cal_event*)malloc(sizeof(cal_event));
  memset(cev,0,sizeof(cal_event));
  calibrate_TIGRESS(data,&cal_par->tg,&cev->tg);
 
  //check the Ge fold
  if(cev->tg.h.FA>0)
    //look through each Tigress position
    for(pos=1;pos<NPOSTIGR;pos++)
      {
        //reset suppression flag for every position
        suppFlag=0;
        //check if the position is in the hit pattern
        //if((cev->tg.h.HHP&(1<<(pos-1)))!=0)
        if((cev->tg.h.HHP&(1<<(pos-1)))!=0)
          if(cev->tg.det[pos].hge.FH>0)
            //check the fold
            if(cev->tg.det[pos].hge.FH>0)
              //check if the position is in the addback hit pattern
              if((cev->tg.h.AHP&(1<<(pos-1)))!=0)
                {
                  //printf("AHP = %d\n",cev->tg.h.AHP);
                  //reset take for add-back suppression
                  take=0;
                  //Run through four cores for each position
                  for(col=0;col<NCOL;col++)
                    {
                    //Check if this color is indicated in the hit pattern
                    if((cev->tg.det[pos].hge.HHP&(1<<col))!=0)
                      //Check that this combination has a fold great than zero
                      if(cev->tg.det[pos].ge[col].h.FH>0)
                        //Check if this combination is indicated in the hit pattern
                        if((cev->tg.det[pos].ge[col].h.HHP&1)!=0)
                          //suppress if the position is in the map and has not yet been suppressed
                          if(cev->tg.det[pos].ge[col].suppress>=supLow && cev->tg.det[pos].ge[col].suppress<=supHigh && take==0)
                            {
                              /* once suppression flag is set
                              do not reset it, could remove the take bit
                              and keep resetting suppFlag, but this
                              is nicer */
                              suppFlag=1;
                              take=1;
                            }
                    }
                    if(suppFlag){
                      //eAddBack = cev->tg.det[pos].addback.E/cal_par->tg.contr_e;
                      eAddBack = cev->tg.det[pos].addback.E;
                      h->Fill(eAddBack,phaseRF);
                      if(rf_period > 0.){
                        hrad->Fill(eAddBack,TWOPI*phaseRF/rf_period);
                      }
                      
                    }
                    
                }
      }
  
  free(cev);
  return SEPARATOR_DISCARD;
}
/*=========================================================================*/
int main(int argc, char *argv[])
{
  FILE *cluster;
  input_names_type* name;
  char title[132];
  char DataFile[132];
  
  if((argc!=4)&&(argc!=5))
    {
      printf("TigressRF_ECalABSuppSum master_file_name supLow supHigh rf_period_ns\n");
      printf("Program sorts RF phase vs. ECalABSuppSum histograms for TIGRESS.\n");
      printf("The RF period may be omitted, in which case the RF phase in radians will not be sorted, only the t0 in ns.\n");
      exit(-1);
    }
  
  printf("Program sorts RF phase vs. ECalABSuppSum histograms for TIGRESS.\n");
  
  h = new TH2D("Tigress vs RF ns","Tigress vs RF ns",S4K,0,S8K-1,90,0,90-1);
  h->Reset();
  hrad = new TH2D("Tigress vs RF rad","Tigress vs RF rad",S4K,0,S8K-1,90,0,TWOPI);
  hrad->Reset();

  name=(input_names_type*)malloc(sizeof(input_names_type));
  memset(name,0,sizeof(input_names_type));
  
  cal_par=(calibration_parameters*)malloc(sizeof(calibration_parameters));
  memset(cal_par,0,sizeof(calibration_parameters));
  
  read_master(argv[1],name);

  supLow = atof(argv[2]);
  supHigh = atof(argv[3]);
  if(argc==5){
    rf_period = atof(argv[4]);
    if(rf_period <= 0){
      printf("ERROR: RF period must be a positive value!\n");
      exit(-1);
    }  
  }else{
    rf_period = -1.; //will be used to skip making an RF phase histogram
  }
  

  if(name->flag.cluster_file==1)
    {
      printf("Sorting ECalABSuppSum histograms for TIGRESS clovers and cores based upon the cluster file: %s\n",name->fname.cluster_file);
      if((cluster=fopen(name->fname.cluster_file,"r"))==NULL)
	{
	  printf("ERROR!!! I can't open input file %s\n",name->fname.cluster_file);
	  exit(-2);
	}
    }
  else
    {
      printf("ERROR!!! Cluster file not defined\n");
      exit(-1);
    }
  
  if(name->flag.TIGRESS_cal_par==1)
    {
      printf("TIGRESS calpar read from: %s\n",name->fname.TIGRESS_cal_par);
      initialize_TIGRESS_calibration(&cal_par->tg,name->fname.TIGRESS_cal_par);
    }
  else
    {
      printf("ERROR!!! TIGRESS calibration parameters not defined\n");
      exit(EXIT_FAILURE);
    }
  
  while(fscanf(cluster,"%s",DataFile) != EOF)
    {
      memset(name,0,sizeof(input_names_type));
      strcpy(name->fname.inp_data,DataFile);
      
      printf("Sorting data from file %s\n", name->fname.inp_data);
      sort(name);
    }
  

  sprintf(title,"TigressRF_ECal.root");
  TFile f(title, "recreate");
  h->GetXaxis()->SetTitle("TIGRESS Energy (keV)");
  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->SetTitle("Tigress RF t0 (ns)");
  h->GetYaxis()->CenterTitle(true);
  h->GetYaxis()->SetTitleOffset(1.5);
  h->SetOption("COLZ");
  gStyle->SetPalette(1);
  h->Write();

  if(rf_period > 0.){
    hrad->GetXaxis()->SetTitle("TIGRESS Energy (keV)");
    hrad->GetXaxis()->CenterTitle(true);
    hrad->GetYaxis()->SetTitle("Tigress RF phase (rad)");
    hrad->GetYaxis()->CenterTitle(true);
    hrad->GetYaxis()->SetTitleOffset(1.5);
    hrad->SetOption("COLZ");
    hrad->Write();
  }

  printf("Output data written to: TigressRF_ECal.root\n");

}
