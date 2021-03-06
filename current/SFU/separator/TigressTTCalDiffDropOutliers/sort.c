#include "sort.h"

int analyze_data(raw_event *data)
{
  cal_event* cev;
  int pos,col;
  int pos1,col1,pos2,col2;
  double thit;
  double t1=0;
  bool keepHit;
  
  long long int one=1,none=-1,kill;
  int id,id_ge;
  long long int flag_ge,drop;

  if((data->h.setupHP&TIGRESS_BIT)==0)
    return SEPARATOR_DISCARD;

  //if((data->h.setupHP&RF_BIT)==0)
  //  return SEPARATOR_DISCARD;

  cev=(cal_event*)malloc(sizeof(cal_event));
  memset(cev,0,sizeof(cal_event));
  calibrate_TIGRESS(data,&cal_par->tg,&cev->tg);

  flag_ge=0;

  if(cev->tg.h.FT>0)
    for(pos1=1;pos1<NPOSTIGR;pos1++)
      if((cev->tg.h.THP&(1<<(pos1-1)))!=0)
	      if(cev->tg.det[pos1].hge.FT>0)
	        for(col1=0;col1<NCOL;col1++)
	          if((cev->tg.det[pos1].hge.THP&(1<<col1))!=0)
	            if(cev->tg.det[pos1].ge[col1].h.FT>0)
		            if((cev->tg.det[pos1].ge[col1].h.THP&1)!=0) //is there a hit in the detector?
		              {
                    t1=cev->tg.det[pos1].ge[col1].seg[0].T/cal_par->tg.contr_t;

                    //now look and see if the hit is within the time window wrt any other hits
                    keepHit = false;

                    for(pos2=1;pos2<NPOSTIGR;pos2++)
                      if(pos2!=pos1)
                        if((cev->tg.h.THP&(1<<(pos2-1)))!=0)
                          if(cev->tg.det[pos2].hge.FT>0)
                            {
                              for(col2=0;col2<NCOL;col2++)
                              {
                                if(col2!=col1)
                                if((cev->tg.det[pos2].hge.THP&(1<<col2))!=0)
                                  if(cev->tg.det[pos2].ge[col2].h.FT>0)
                                    if((cev->tg.det[pos2].ge[col2].h.THP&1)!=0) //is there a hit in the detector?
                                      {
                                        thit=cev->tg.det[pos2].ge[col2].seg[0].T/cal_par->tg.contr_t;
                                        if(fabs(thit-t1) < gate_length)
                                          {
                                            keepHit = true;
                                            break;
                                          }
                                      }
                              }
                              if(keepHit)
                                break;
                            }
                    
                    if(keepHit){
                      id=pos1-1;
                      id_ge=id*NCOL+col1;
                      flag_ge|=(one<<id_ge); //flag the correlated hit for preservation
                    }
                    
                  }
  //printf("t1=%f\n",t1);

    free(cev);
    
    //drop TIGRESS data out of the time limits
    for(pos=1;pos<NPOSTIGR;pos++)
      {
	      id=pos-1;
	      for(col=0;col<NCOL;col++)
	        {	
	          id_ge=id*NCOL+col;
	          drop=(one<<id_ge);
	          drop&=data->tg.g.GeHP;
	          if(drop!=0)
	            {
		            drop&=flag_ge;
		            if(drop==0)
		              {
		                //drop this crystal
		                memset(&data->tg.det[pos].ge[col],0,sizeof(SegTIGR));
		                kill=none-(one<<col);
		                data->tg.det[pos].h.GeHP&=kill;
		                data->tg.det[pos].h.Gefold--;
		                kill=none-(one<<id_ge);
		                data->tg.g.GeHP&=kill;
		                data->tg.g.Gefold--;
		                data->tg.g.THP&=kill;
		                data->tg.g.Tfold--;
		              }
	            }
	        }
      }
    
    for(pos=1;pos<NPOSTIGR;pos++)
      {
	      id=pos-1;
	      if((data->tg.h.GeHP&(1<<id))!=0)
	        if(data->tg.det[pos].h.Gefold<=0)
	          {
	            //drop this position
	            memset(&data->tg.det[pos],0,sizeof(CssTIGR));
	            kill=none-(one<<id);
	            data->tg.h.GeHP&=kill;
	            data->tg.h.Gefold--;
	            data->tg.g.PosHP&=kill;
	            data->tg.g.Posfold--;
	          }
      }    
    
    if(data->tg.h.Gefold<=0)
      {
	      memset(&data->tg,0,sizeof(Tigress));
	      kill=none-TIGRESS_BIT;
	      data->h.setupHP&=kill;
      }
    
    
    if((data->h.setupHP&TIGRESS_BIT)==0)
      return SEPARATOR_DISCARD;
    
    //if((data->h.setupHP&RF_BIT)==0)
    //  return SEPARATOR_DISCARD; 
    
    //check that there are at least 2 Tigress events 
		if(data->tg.h.Gefold<2)
			return SEPARATOR_DISCARD;
    
    encode(data,output,enb);
    
    return SEPARATOR_DISCARD;
}
/*====================================================================================*/
int main(int argc, char *argv[])
{
  input_names_type* name;
  //TCanvas *canvas;
  //TApplication *theApp;

  if(argc!=2)
    {
      printf("separate_Tigress_TTCalDiffDropOutliers master_file_name\n");
      printf("\n Separates out gamma-gamma coincidences where the gammas arrive within (calibrated) time gate of gate_length.  Discards any gamma hits which arrive outside of the time gate (with respect to all other hits).\nThe time gate length is specified in the Tigress array calibration parameters file (under 'TIGRESS_TTCal_gate_length').\n");
      exit(-1);
    }
    
  printf("Program sorts data separated on Tigress timing.\n");
  
  gate_length=0;
  
  name=(input_names_type*)malloc(sizeof(input_names_type));
  memset(name,0,sizeof(input_names_type));
  cal_par=(calibration_parameters*)malloc(sizeof(calibration_parameters));
  memset(cal_par,0,sizeof(calibration_parameters));

  read_master(argv[1],name);

  if(name->flag.inp_data!=1)
    {
      printf("\nInput data file not defined\n");
      exit(EXIT_FAILURE);
    }

  if(name->flag.TIGRESS_cal_par==1)
        {
          printf("\nTIGRESS calibration read from the file:\n %s\n",name->fname.TIGRESS_cal_par);
          initialize_TIGRESS_calibration(&cal_par->tg,name->fname.TIGRESS_cal_par);
	  
        }
      else
        {
          printf("\nTIGRESS calibration parameters not defined\n");
          exit(-1);
        }

  if(name->flag.out_data!=1)
    {
      printf("\nOutput data file not defined\n");
      exit(EXIT_FAILURE);
    }
  
  if((output=fopen(name->fname.out_data,"w"))==NULL)
    {
      printf("\nI can't open output file %s for writing\n",name->fname.out_data);
      exit(-2);
    }
  memset(enb,0,sizeof(enb));
  enb[0]=BUFFER_TAG;
  enb[1]++;
  enb[1]++;
  
  gate_length=cal_par->tg.TTCal_gate_length;
  printf("Gate length is %f\n",gate_length);

  sort(name);
  
	//save the last buffer which will be dropped otherwise
	//if enb[1]==2 then the buffer contains no data, only the header
	if(enb[1]>2)
		fwrite(enb,sizeof(int),BUFFSIZE,output);

  /*theApp=new TApplication("App", &argc, argv); 
  canvas = new TCanvas("TTCalDiff", "TTCalDiff",10,10, 500, 300); 
  gPad->SetLogz(1); 
  g->SetLineColor(2); 
  h->Draw(); 
  g->Draw("same"); 
  theApp->Run(kTRUE);*/ 
}

