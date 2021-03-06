#include "sort.h"

int analyze_data(raw_event *data)
{
  cal_event* cev;
  int pos,col,csi;
  double tdiff;
  bool keepHit;
  
  int64_t one=1,none=-1,kill;
  int id,id_ge;
  long long int flag_ge,drop;
  int64_t flag_csi[4];

  double ttg=10E10;
  double tcsi=-10E10;

  /* printf("init bits:\n  Tigress: %d\n  CsI: %d\n  RF: %d\n",data->h.setupHP&TIGRESS_BIT,data->h.setupHP&CsIArray_BIT,data->h.setupHP&RF_BIT); */
  /* getc(stdin); */

  if((data->h.setupHP&TIGRESS_BIT)==0)
    return SEPARATOR_DISCARD;

  //if((data->h.setupHP&RF_BIT)==0)
  //  return SEPARATOR_DISCARD;

  if((data->h.setupHP&CsIArray_BIT)==0)
    return SEPARATOR_DISCARD;

  cev=(cal_event*)malloc(sizeof(cal_event));
  memset(cev,0,sizeof(cal_event));
  calibrate_TIGRESS(data,&cal_par->tg,&cev->tg);
  calibrate_CSIARRAY(data,&cal_par->csiarray,&cev->csiarray);

  flag_ge=0;
  memset(flag_csi,0,sizeof(flag_csi));

  if(high>low)
    {

      if(cev->tg.h.FT>0)
        for(pos=1;pos<NPOSTIGR;pos++) //look at each Tigress position
        	{
		        if((cev->tg.h.THP&(1<<(pos-1)))!=0)
		          if(cev->tg.det[pos].hge.FT>0)
		            for(col=0;col<NCOL;col++)
		              if((cev->tg.det[pos].hge.THP&(1<<col))!=0)
		                if(cev->tg.det[pos].ge[col].h.FT>0)
											if((cev->tg.det[pos].ge[col].h.THP&1)!=0) //is there a hit in the detector?
												{
                          ttg=cev->tg.det[pos].ge[col].seg[0].T/cal_par->tg.contr_t; //get hit time

                          //now look and see if the hit is within the time window wrt ALL CsI hits
                          keepHit = true;
                          if(cev->csiarray.h.FT>0)
                            {
                              for(csi=1;csi<NCSI;csi++) //look at each CsI position
                                {
                                  if((cev->csiarray.h.THP[csi/64]&(one<<csi%64))!=0) //is there a hit in the detector?
                                    {
                                      flag_csi[csi/64]|=(one<<csi%64); //flag every CsI hit for preservation
                                      tcsi=cev->csiarray.csi[csi].T/cal_par->csiarray.contr_t;
                                      tdiff = fabs(tcsi - ttg);
                                      if((tdiff > high)||(tdiff < low))
                                        keepHit = false;
                                    }
                                    
                                }
                            }
                          else
                            {
                              keepHit = false;
                            }
                            

                          if(keepHit){
                            id=pos-1;
                            id_ge=id*NCOL+col;
                            flag_ge|=(one<<id_ge); //flag the hit for preservation
                          }
													
												}
					}
    }

  // cev malloc
  // disasterous memory leak!
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
  
  //drop csi out of the time limits
  for(csi=1;csi<NCSI;csi++)
    if((data->csiarray.h.TSHP[csi/64]&(one<<csi%64))!=0)
      if((flag_csi[csi/64]&(one<<csi%64))==0)
	      {
	        memset(&data->csiarray.csi[csi],0,sizeof(channel));
	        memset(&data->csiarray.wfit[csi],0,sizeof(ShapePar));
	        memset(&data->csiarray.t0[csi],0,sizeof(double));
	        data->csiarray.h.Efold--;
	        data->csiarray.h.Tfold--;	  
	        data->csiarray.h.TSfold--;
	        kill=none-(one<<csi%64);
	        data->csiarray.h.TSHP[csi/64]&=kill;
	        data->csiarray.h.EHP[csi/64]&=kill;
	        data->csiarray.h.THP[csi/64]&=kill;
	      }
  
  if(data->csiarray.h.TSfold<=0)
    {
      kill=none-CsIArray_BIT; 
      data->h.setupHP&=kill;
      memset(&data->csiarray,0,sizeof(CsIArray));
    } 

  /* printf("final bits:\n  Tigress: %d\n  CsI: %d\n  RF: %d\n",data->h.setupHP&TIGRESS_BIT,data->h.setupHP&CsIArray_BIT,data->h.setupHP&RF_BIT); */

  if((data->h.setupHP&TIGRESS_BIT)==0)
    return SEPARATOR_DISCARD;
  
  /* if((data->h.setupHP&RF_BIT)==0) */
   /*   return SEPARATOR_DISCARD; */
  
  if((data->h.setupHP&CsIArray_BIT)==0)
      return SEPARATOR_DISCARD;
      
  /* printf("pre encode:\n  Tigress: %d\n  CsI: %d\n  RF: %d\n",data->h.setupHP&TIGRESS_BIT,data->h.setupHP&CsIArray_BIT,data->h.setupHP&RF_BIT); */
  /* getc(stdin); */
  
  encode(data,output,enb);
  
  return SEPARATOR_DISCARD;
}
/*====================================================================================*/
int main(int argc, char *argv[])
{
  input_names_type* name;

  if(argc!=4)
    {
      printf("\n ./separate_TigressCsIArray_TTCalDropOutliers master_file_name low high\n");
      printf("\n Sorts events where the all gammas come within the time window defined by 'low' and 'high' with respect to all CsI hits.\n");
      exit(-1);
    }

  low=atof(argv[2]);
  high=atof(argv[3]);
  
  printf("Program sorts data separated on Tigress-CsIArray timing.\n");
	if(low>high)
		{
			//swap values
			double tmp=low;
			low=high;
			high=tmp;
		}
  

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

  if(name->flag.CSIARRAY_cal_par==1)
    {
      printf("\nCSIARRAY calibration read from the file:\n %s\n",name->fname.CSIARRAY_cal_par);
      initialize_CSIARRAY_calibration(&cal_par->csiarray,name->fname.CSIARRAY_cal_par);
    }
  else
    {
      printf("\nCSIARRAY calibration parameters not defined\n");
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

  sort(name);
  
	//save the last buffer which will be dropped otherwise
	//if enb[1]==2 then the buffer contains no data, only the header
	if(enb[1]>2)
		fwrite(enb,sizeof(int),BUFFSIZE,output);

  free(name);
  free(cal_par);
 
}

