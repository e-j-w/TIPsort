#include "SFU-cal-CSIARRAY.h"
/*******************************************************************/
void initialize_CSIARRAY_calibration(CSIARRAY_calibration_parameters *CSIARRAY_cal_par, char *name)
{
  FILE *inp;
  char str1[256],str2[256];

  memset(CSIARRAY_cal_par,0,sizeof(*CSIARRAY_cal_par));
  CSIARRAY_cal_par->contr_e=1.;
  CSIARRAY_cal_par->contr_t=1.;
  CSIARRAY_cal_par->use_time_fit=0;
  //For RFUnwrapping                                                                                                                                                                                                                        
  CSIARRAY_cal_par->DoRFUnwrapping=0;

  if((inp=fopen(name,"r"))==NULL)
      {
         printf("\nFile %s can not be opened\n",name);
         exit(EXIT_FAILURE);
      }
  else
    printf("CSIARRAY calibration parameters read from file %s\n",name);

 
  if(fgets(str1,256,inp)!=NULL)
    {
      if(fgets(str1,256,inp)!=NULL)
	{
	  while(fscanf(inp,"%s %s",str1,str2)!=EOF)
	    {

	      //printf("%s %s\n",str1,str2);
	      //getc(stdin);
	      if(strcmp(str1,"Energy_spectra_contraction_factor")==0)
		CSIARRAY_cal_par->contr_e=atof(str2);	      
	      
	      if(strcmp(str1,"Time_spectra_contraction_factor")==0)
		CSIARRAY_cal_par->contr_t=atof(str2);
	      
	      if(strcmp(str1,"Energy_calibration_parameters")==0)
		read_CSIARRAY_e_cal_par(CSIARRAY_cal_par,str2);
	      
	      if(strcmp(str1,"Time_calibration_parameters")==0)
		read_CSIARRAY_t_cal_par(CSIARRAY_cal_par,str2);

	      if(strcmp(str1,"Energy_limits")==0)
		read_CSIARRAY_energy_limits(CSIARRAY_cal_par,str2);

	      if(strcmp(str1,"Time_limits")==0)
		read_CSIARRAY_time_limits(CSIARRAY_cal_par,str2);

	      if(strcmp(str1,"Detector_position_file")==0)
		read_CSIARRAY_detector_positions(CSIARRAY_cal_par,str2);
	      
	      if(strcmp(str1,"DeltaU_file")==0)
		read_CSIARRAY_deltaU_parameters(CSIARRAY_cal_par,str2);
	      
	      if(strcmp(str1,"Timing_from")==0)
		{
		  if(strcmp(str2,"fit_function")==0)
		    CSIARRAY_cal_par->use_time_fit=1;
		  if(strcmp(str2,"FIT_FUNCTION")==0)
		    CSIARRAY_cal_par->use_time_fit=1;
		  if(strcmp(str2,"Fit_function")==0)
		    CSIARRAY_cal_par->use_time_fit=1;
		  if(strcmp(str2,"Fit_Function")==0)
		    CSIARRAY_cal_par->use_time_fit=1;

		  if(strcmp(str2,"fit_t0")==0)
		    CSIARRAY_cal_par->use_time_fit=2;
		  if(strcmp(str2,"FIT_T0")==0)
		    CSIARRAY_cal_par->use_time_fit=2;
		  if(strcmp(str2,"Fit_t0")==0)
		    CSIARRAY_cal_par->use_time_fit=2;
		  if(strcmp(str2,"Fit_T0")==0)
		    CSIARRAY_cal_par->use_time_fit=2;
		}

	      //For RFUnwrapping                                                                                                                                                                                                            
              if(strcmp(str1,"RFUnwrapping_option")==0)
                {
                  if(strcmp(str2,"YES")==0)
                    CSIARRAY_cal_par->DoRFUnwrapping=1;
                  if(strcmp(str2,"Yes")==0)
                    CSIARRAY_cal_par->DoRFUnwrapping=1;
                  if(strcmp(str2,"yes")==0)
                    CSIARRAY_cal_par->DoRFUnwrapping=1;
                  if(strcmp(str2,"y")==0)
                    CSIARRAY_cal_par->DoRFUnwrapping=1;
                  if(strcmp(str2,"Y")==0)
                    CSIARRAY_cal_par->DoRFUnwrapping=1;
                }

              if(strcmp(str1,"CSIARRAY_RFunwrapping_offset")==0)
                {
                  CSIARRAY_cal_par->offset=atoi(str2);
                  printf("RFunwrapping CSIARRAY_offset is %d\n",CSIARRAY_cal_par->offset);
                }

              if(strcmp(str1,"CSIARRAY_RFunwrapping_shift")==0)
                {
                  CSIARRAY_cal_par->shift=atoi(str2);
                  printf("RFunwrapping CSIARRAY_shift is %d\n",CSIARRAY_cal_par->shift);
                }
	    }
	}
    }
  else
    {
      printf("Wrong structure of file %s\n",name);
      printf("Aborting sort\n");
      exit(1);
    }    
  fclose(inp);
    
}
/***********************************************************************/
void read_CSIARRAY_e_cal_par(CSIARRAY_calibration_parameters *CSIARRAY_cal_par, char *filename)
{
  FILE *inp;
  char line[132];
  int  pos,i;
  double a[3];

  if((inp=fopen(filename,"r"))==NULL)
      {
         printf("\nI can't open file %s\n",filename);
         exit(EXIT_FAILURE);
      }
  printf("\nCSIARRAY energy calibration parameters read from the file:\n %s\n",filename);

  if(fgets(line,132,inp)!=NULL)
    {
      if(fgets(line,132,inp)!=NULL)
	//while(fscanf(inp,"%d %lf %lf",&pos,&a[0],&a[1])!=EOF)
	while(fscanf(inp,"%d %lf %lf %lf",&pos,&a[0],&a[1],&a[2])!=EOF)
	  if(pos>0&&pos<NCSI)
	    {
	      CSIARRAY_cal_par->ceflag[pos]=1;
	      for(i=0;i<3;i++) CSIARRAY_cal_par->ce[pos][i]=a[i];
	    }
    }
  else
    {
      printf("Wrong structure of file %s\n",filename);
      printf("Aborting sort\n");
      exit(1);
    }
  fclose(inp);  
}
/***********************************************************************/
void read_CSIARRAY_t_cal_par(CSIARRAY_calibration_parameters *CSIARRAY_cal_par, char *filename)
{
  FILE *inp;
  char line[132];
  int  pos,i;
  double a[2];

  if((inp=fopen(filename,"r"))==NULL)
      {
         printf("\nI can't open file %s\n",filename);
         exit(EXIT_FAILURE);
      }
  printf("\nCSIARRAY time calibration parameters read from the file:\n %s\n",filename);

  if(fgets(line,132,inp)!=NULL)
    {
      if(fgets(line,132,inp)!=NULL)
	while(fscanf(inp,"%d %lf %lf",&pos,&a[0],&a[1])!=EOF)
	  if(pos>0&&pos<NCSI)
	    {
	      CSIARRAY_cal_par->ctflag[pos]=1;
	      for(i=0;i<2;i++) CSIARRAY_cal_par->ct[pos][i]=a[i];
	    }
    }
  else
    {
      printf("Wrong structure of file %s\n",filename);
      printf("Aborting sort\n");
      exit(1);
    }
  fclose(inp);
}
/***********************************************************************/
void read_CSIARRAY_energy_limits(CSIARRAY_calibration_parameters *CSIARRAY_cal_par, char *filename)
{
  FILE *inp;
  char line[132];
  int  pos;
  float low,high;
 

  if((inp=fopen(filename,"r"))==NULL)
      {
         printf("\nI can't open file %s\n",filename);
         exit(EXIT_FAILURE);
      }
  printf("\nCSIARRAY energy limits read from the file:\n %s\n",filename);

  if(fgets(line,132,inp)!=NULL)
    {
      if(fgets(line,132,inp)!=NULL)
	while(fscanf(inp,"%d %f %f",&pos,&low,&high)!=EOF)
	  if(pos>0&&pos<NCSI)
	      {
		CSIARRAY_cal_par->elow[pos]=low;
		CSIARRAY_cal_par->ehigh[pos]=high;
	      }
    }
  else
    {
      printf("Wrong structure of file %s\n",filename);
      printf("Aborting sort\n");
      exit(1);
    }
  fclose(inp);
  
}
/***********************************************************************/
void read_CSIARRAY_time_limits(CSIARRAY_calibration_parameters *CSIARRAY_cal_par, char *filename)
{
  FILE *inp;
  char line[132];
  int  pos;
  float low,high;
 

  if((inp=fopen(filename,"r"))==NULL)
      {
         printf("\nI can't open file %s\n",filename);
         exit(EXIT_FAILURE);
      }
  printf("\nCSIARRAY time limits read from the file:\n %s\n",filename);

  if(fgets(line,132,inp)!=NULL)
    {
      if(fgets(line,132,inp)!=NULL)
	while(fscanf(inp,"%d %f %f",&pos,&low,&high)!=EOF)
	  if(pos>0&&pos<NCSI)
	      {
		CSIARRAY_cal_par->tlow[pos]=low;
		CSIARRAY_cal_par->thigh[pos]=high;
	      }
    }
  else
    {
      printf("Wrong structure of file %s\n",filename);
      printf("Aborting sort\n");
      exit(1);
    }
  fclose(inp);
  
}

/***********************************************************************/
void read_CSIARRAY_detector_positions(CSIARRAY_calibration_parameters *CSIARRAY_cal_par, char *filename)
{
  FILE *inp;
  char line[132];
  int  pos;
  double x,y,z; 
  double R,theta,phi;

  if((inp=fopen(filename,"r"))==NULL)
      {
         printf("\nI can't open file %s\n",filename);
         exit(EXIT_FAILURE);
      }
  printf("\nCSIARRAY detector positions read from the file:\n %s\n",filename);

  if(fgets(line,132,inp)!=NULL)
    {
      if(fgets(line,132,inp)!=NULL)
	while(fscanf(inp,"%d %lf %lf %lf",&pos,&x,&y,&z)!=EOF)
	  if(pos>0&&pos<NCSI)
	    {
	      CSIARRAY_cal_par->cposflag[pos]=1;

	      R = sqrt(x*x+y*y+z*z);
	      theta = acos(z/R);
	      phi = atan2(y,x);
	      
	      CSIARRAY_cal_par->cpos[pos][0]=R;
	      CSIARRAY_cal_par->cpos[pos][1]=theta;
	      CSIARRAY_cal_par->cpos[pos][2]=phi;
	    }
    }
  else
    {
      printf("Wrong structure of file %s\n",filename);
      printf("Aborting sort\n");
      exit(1);
    }
  fclose(inp);
}

/***********************************************************************/
void read_CSIARRAY_deltaU_parameters(CSIARRAY_calibration_parameters *CSIARRAY_cal_par, char *filename)
{
  FILE *inp;
  char line[132];
  char str1[256],str2[256];
  
  if((inp=fopen(filename,"r"))==NULL)
      {
         printf("\nI can't open file %s\n",filename);
         exit(EXIT_FAILURE);
      }
  printf("\nDeltaU calculation parameters read from the file:\n %s\n",filename);

  if(fgets(line,132,inp)!=NULL)
    {
      if(fgets(line,132,inp)!=NULL)
	while(fscanf(inp,"%s %s",&str1,&str2)!=EOF)
	  {
	    if(strcmp(str1,"Beam_energy_for_deltaU_in_MeV")==0)
	      CSIARRAY_cal_par->Ebeam=atof(str2);
	    
	    if(strcmp(str1,"Projectile_mass_in_amu")==0)
	      CSIARRAY_cal_par->mproj=amuMeV*atoi(str2);
	    
	    if(strcmp(str1,"Target_mass_in_amu")==0)
	      CSIARRAY_cal_par->mt=amuMeV*atoi(str2);
	    
	    if(strcmp(str1,"Particle_mass_in_amu")==0)
	      CSIARRAY_cal_par->mp=amuMeV*atoi(str2);
	    
	    if(strcmp(str1,"Daughter_mass_in_amu")==0)
	      CSIARRAY_cal_par->md=amuMeV*atoi(str2);
	  }
    }
  else
    {
      printf("Wrong structure of file %s\n",filename);
      printf("Aborting sort\n");
      exit(1);
    }

  /* printf("Ebeam %f mproj %f mt %f mp %f md %f\n",CSIARRAY_cal_par->Ebeam,CSIARRAY_cal_par->mproj,CSIARRAY_cal_par->mt,CSIARRAY_cal_par->mp,CSIARRAY_cal_par->md); */
  /* getc(stdin); */

  fclose(inp);
}

/*******************************************************************/
void summarize_CSIARRAY_calibration(CSIARRAY_calibration_parameters *CSIARRAY_cal_par, char *filename)
{
  FILE *out;
  int  pos;

  if((out=fopen(filename,"w"))==NULL)
      {
         printf("\nFile %s can not be opened\n",filename);
         exit(EXIT_FAILURE);
      }

  fprintf(out,"CSIARRAY energy calibration parameters\n");
  fprintf(out,"  pos#       a0             a1              a2\n");
  for(pos=1;pos<NCSI;pos++)
      if(CSIARRAY_cal_par->ceflag[pos]==1)
	fprintf(out,"   %02d  %12.5e   %12.5e   %12.5e\n",pos,CSIARRAY_cal_par->ce[pos][0],CSIARRAY_cal_par->ce[pos][1],CSIARRAY_cal_par->ce[pos][2]);
     
  fprintf(out,"CSIARRAY time calibration parameters\n");
  fprintf(out,"  pos#  a0             a1 \n");
  for(pos=1;pos<NCSI;pos++)
    if(CSIARRAY_cal_par->ctflag[pos]==1)
      fprintf(out,"   %02d  %12.5e   %12.5e\n",pos,CSIARRAY_cal_par->ct[pos][0],CSIARRAY_cal_par->ct[pos][1]);
  fclose(out);

 fprintf(out,"CSIARRAY detector positions\n");
  fprintf(out,"  pos#       x             y              z\n");
  for(pos=1;pos<NCSI;pos++)
      if(CSIARRAY_cal_par->cposflag[pos]==1)
	fprintf(out,"   %02d  %12.1f   %12.1f   %12.1f\n",pos,CSIARRAY_cal_par->cpos[pos][0],CSIARRAY_cal_par->cpos[pos][1],CSIARRAY_cal_par->cpos[pos][2]);
    
}
/********************************************************************/
void calibrate_CSIARRAY(raw_event* rev, CSIARRAY_calibration_parameters *CSIARRAY_cal_par, cCSIARRAY *cp)
{
  unsigned long long int one=1; 
  int pos,pos1,pos2;
  double ran,ren,e,t;

  //For RFUnwrapping
  double trf,func;
  long int a;

  //For excitation energy change
  double theta1,phi1,theta2,phi2;
  double T1,T2,Vcm;
  double mproj,mt,ma,md;
  double Ebeam;
  //double R1,R2,dU;

  memset(cp,0,sizeof(cCSIARRAY));
  /* energy calibration */
  if(rev->csiarray.h.Efold>0)
    for(pos=1;pos<NCSI;pos++)
      if(CSIARRAY_cal_par->ceflag[pos]==1)
	if((rev->csiarray.h.THP&(one<<pos))!=0)
	  if(rev->csiarray.wfit[pos].am[1]>0)
  	    {
  	      ran=(double)rand()/(double)RAND_MAX-0.5;
  	      ren=rev->csiarray.csi[pos].charge+ran;
  	      //e=CSIARRAY_cal_par->ce[pos][0]+CSIARRAY_cal_par->ce[pos][1]*ren;
	      //1000* converts to 1keV/ch
  	      e=1000*(CSIARRAY_cal_par->ce[pos][0]*ren + CSIARRAY_cal_par->ce[pos][2]*(log(CSIARRAY_cal_par->ce[pos][1]*ren) + 1)); 
	      
	      if(e>0)
		if(e<S65K)
		  {
		    cp->csi[pos].E=e;
		    cp->h.EHP|=(one<<pos);
		    cp->h.FE++;

		    //printf("for pos %d parameters g = %.4f b = %.4f a = %.4f\n",pos,CSIARRAY_cal_par->ce[pos][0],CSIARRAY_cal_par->ce[pos][1], CSIARRAY_cal_par->ce[pos][2]);
		    //printf("energy is %f for charge %f\n",e,ren);
		    //getc(stdin);
		  }
 	    }

  /* DAQ CFD time calibration */
  if(CSIARRAY_cal_par->use_time_fit==0)
      if(rev->csiarray.h.Tfold>0)
	for(pos=1;pos<NCSI;pos++)
	  if(CSIARRAY_cal_par->ctflag[pos]==1)
	    if((rev->csiarray.h.THP&(one<<pos))!=0)
	      if(rev->csiarray.csi[pos].cfd>0)
		{
		  t=rev->csiarray.csi[pos].cfd&0x00ffffff;
		  t-=(rev->csiarray.csi[pos].timestamp*16)&0x00ffffff;

		  if((rev->h.setupHP&RF_BIT)!=0)
		    {
		      //For RFUnwrapping
		      trf=rev->rf.sin.t0;
		      if(CSIARRAY_cal_par->DoRFUnwrapping==1)
			{
			  
			  func=(1.000*trf)-(RFphase+CSIARRAY_cal_par->offset);
			  a=fmod(t+CSIARRAY_cal_par->shift,RFphase);
			  if(a>func) trf+=RFphase;
			}
		    }
		  else
		    trf=0.;

		  t-=trf;
		  //The line below was replaced by the line above and below the above if statement
		  //t-=rev->rf.sin.t0;
		  t+=S16K;

		  ran=(double)rand()/(double)RAND_MAX-0.5;
		  ren=t+ran;		
		  t=CSIARRAY_cal_par->ct[pos][0]+CSIARRAY_cal_par->ct[pos][1]*ren;
		  if(t>0)
		    if(t<S65K)
		      {
			cp->csi[pos].T=t;
			cp->h.THP|=(one<<pos);
			cp->h.FT++;
		      }
		}

  /* fit function time calibration */
  if(CSIARRAY_cal_par->use_time_fit==1)
    for(pos=1;pos<NCSI;pos++)
      if(CSIARRAY_cal_par->ctflag[pos]==1)
	if((rev->csiarray.h.TSHP&(one<<pos))!=0)
	  {
	    t=rev->csiarray.wfit[pos].t[0]*16;
	    
	    if((rev->h.setupHP&RF_BIT)!=0)
	      {
		//For RFUnwrapping
		trf=rev->rf.sin.t0;
		if(CSIARRAY_cal_par->DoRFUnwrapping==1)
		  {
		    func=(1.000*trf)-(RFphase+CSIARRAY_cal_par->offset);
		    a=fmod(t+CSIARRAY_cal_par->shift,RFphase);
		    if(a>func) trf+=RFphase;
		  }
	      }
	    else
	      trf=0.;

	    t-=trf;
	    //The line below was replaced by the line above and below the above if statement
	    //t-=rev->rf.sin.t0;
	    t+=S16K;
	    
	    ran=(double)rand()/(double)RAND_MAX-0.5;
	    ren=t+ran;		
	    t=CSIARRAY_cal_par->ct[pos][0]+CSIARRAY_cal_par->ct[pos][1]*ren;
	    //printf("calibrator CsI time: %.4f for raw time %.4f\n",t,ren);
	    //getc(stdin);
	    if(t>0)
	      if(t<S65K)
		{
		  cp->csi[pos].T=t;
		  cp->h.THP|=(one<<pos);
		  cp->h.FT++;
		}
	  }
  
  /* fit t0 time calibration */
  if(CSIARRAY_cal_par->use_time_fit==2)
      for(pos=1;pos<NCSI;pos++)
	if(CSIARRAY_cal_par->ctflag[pos]==1)
	  if((rev->csiarray.h.TSHP&(one<<pos))!=0)
	    {
	      t=rev->csiarray.t0[pos]*16;

	      if((rev->h.setupHP&RF_BIT)!=0)
		{
		  //For RFUnwrapping
		  trf=rev->rf.sin.t0;
		  if(CSIARRAY_cal_par->DoRFUnwrapping==1)
		    {
		      func=(1.000*trf)-(RFphase+CSIARRAY_cal_par->offset);
		      a=fmod(t+CSIARRAY_cal_par->shift,RFphase);
		      if(a>func) trf+=RFphase;
		    }
		}
	      else
		trf=0.;
	      
              t-=trf;
              //The line below was replaced by the line above and below the above if statement
	      //t-=rev->rf.sin.t0;
	      t+=S16K;

	      ran=(double)rand()/(double)RAND_MAX-0.5;
	      ren=t+ran;		
	      t=CSIARRAY_cal_par->ct[pos][0]+CSIARRAY_cal_par->ct[pos][1]*ren;
	      if(t>0)
		if(t<S65K)
		  {
		    cp->csi[pos].T=t;
		    cp->h.THP|=(one<<pos);
		    cp->h.FT++;
		  }
	    }
  
  if(cp->h.FE>0)
    if(cp->h.FT>0)
      for(pos=1;pos<NCSI;pos++)
	if((cp->h.EHP&(one<<pos))!=0)
	  if((cp->h.THP&(one<<pos))!=0)
	    {
	      cp->h.HHP|=(one<<pos);
	      cp->h.FH++;
	    }

  //initialize the excitation energy change U = -1
  cp->U=-1.;

  //initialize angles and energies
  T1=0.;
  T2=0.;
  theta1=0.;
  theta2=0.;
  phi1=0.;
  phi2=0.;

  /* if(cp->h.FH!=2) */
  /*   { */
  /*     printf("FE %d FT %d FH %d\n",cp->h.FE,cp->h.FT,cp->h.FH); */
  /*     getc(stdin); */
  /*   } */

  //change U for fold 2 events based on energies and positions
  if(cp->h.FH==2)
      for(pos1=1;pos1<NCSI;pos1++)
	if((cp->h.THP&(one<<pos1))!=0)
	  for(pos2=pos1+1;pos2<NCSI;pos2++)
	    if((cp->h.THP&(one<<pos2))!=0)
	      {
	      if(CSIARRAY_cal_par->cposflag[pos1]==1)
		if(CSIARRAY_cal_par->cposflag[pos2]==1)
		  {
		    //theta and phi of the emitted particles
		    theta1 = CSIARRAY_cal_par->cpos[pos1][1];
		    phi1 = CSIARRAY_cal_par->cpos[pos1][2];
		    theta2 = CSIARRAY_cal_par->cpos[pos2][1];
		    phi2 = CSIARRAY_cal_par->cpos[pos2][2];
		    
		    //converts to MeV
		    T1=cp->csi[pos1].E/1000.;
		    T2=cp->csi[pos2].E/1000.;
		  }

	      Ebeam=CSIARRAY_cal_par->Ebeam;
	      mproj=CSIARRAY_cal_par->mproj;
	      mt=CSIARRAY_cal_par->mt;
	      ma=CSIARRAY_cal_par->mp; //alpha mass is particle mass mp
	      md=CSIARRAY_cal_par->md;

	      Vcm = sqrt(2*mproj*Ebeam)/(mproj+mt);
	      
	      /* printf("Ebeam %f Vcm %f mproj %f mt %f ma %f md %f T1 %f T2 %f\n",Ebeam,Vcm,mproj,mt,ma,md,T1,T2); */

	      /* dU = (T1+T2)*(1+(ma/md)) + (1+2*(ma/md))*ma*Vcm*Vcm + 2*(ma/md)*sqrt(T1*T2)*(sin(theta1)*sin(theta2)*cos(phi1-phi2)+cos(theta1)*cos(theta2)) - (1+2*(ma/md))*sqrt(2*ma)*Vcm*(sqrt(T1)*cos(theta1)+sqrt(T2)*cos(theta2)); */

	      //converts back to keV (2 keV/ch contraction)
	      cp->U=500*((T1+T2)*(1+(ma/md)) + (1+2*(ma/md))*ma*Vcm*Vcm + 2*(ma/md)*sqrt(T1*T2)*(sin(theta1)*sin(theta2)*cos(phi1-phi2)+cos(theta1)*cos(theta2)) - (1+2*(ma/md))*sqrt(2*ma)*Vcm*(sqrt(T1)*cos(theta1)+sqrt(T2)*cos(theta2)));   
	      }

  /* printf("dU [MeV] %f cp->U [2 keV/ch] %f\n",dU,cp->U); */
  /* getc(stdin);   */
}
