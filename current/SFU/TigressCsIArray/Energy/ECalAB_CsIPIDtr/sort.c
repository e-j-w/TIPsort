#include "sort.h"

int analyze_data(raw_event *data)
{
  cal_event* cev;
  uint64_t one=1;
  int type;
  double eAddBack=0;
  double chi,s,f,r,t;
  int na,np,sp;
  
  
  //check if all the elements of the setup are present; discard if not
  if((data->h.setupHP&RF_BIT)==0) return SEPARATOR_DISCARD;
  if((data->h.setupHP&CsIArray_BIT)==0) return SEPARATOR_DISCARD;
  if((data->h.setupHP&TIGRESS_BIT)==0) return SEPARATOR_DISCARD;
  
  ev++;
  na=0;
  np=0;
  
  /* work out number and type of particles identified in the event */
  for(pos=1;pos<NCSI;pos++)
    if((aGateFlag[pos]==1)||(pGateFlag[pos]==1))
      if((data->csiarray.h.TSHP[pos/64]&(one<<pos%64))!=0)
	{
	  type=data->csiarray.wfit[pos].type;
	  chi=data->csiarray.wfit[pos].chisq;
	  chi/=data->csiarray.wfit[pos].ndf;
	  if((type>=idmin) && (type<=idmax))
	    if((chi>=chimin) && (chi<=chimax))
	      {
		t=16.*data->csiarray.wfit[pos].t[0]; 
		t-=data->rf.sin.t0;		      
		t+=S4K;
		t/=16;
		s=data->csiarray.wfit[pos].am[3];
		f=data->csiarray.wfit[pos].am[2];
		r=100+s/f*100;
		if(pGateFlag[pos]==1)
		  {
		    if(pGate[pos]->IsInside(t,r)) np++;
		  }
		if(aGateFlag[pos]==1)
		  {
		    if(aGate[pos]->IsInside(t,r)) na++;
		  }}}

  if((np<pMax) && (na<aMax)) ch[np][na]++;
  sp=id[np][na];

  /* now calibrate TIGRESS (earlier is just a waste of CPU) */
  cev=(cal_event*)malloc(sizeof(cal_event));
  memset(cev,0,sizeof(cal_event));
  calibrate_TIGRESS(data,&cal_par->tg,&cev->tg);
  
  if(cev->tg.h.FH>0)
    for(pos=1;pos<NPOSTIGR;pos++)
      if(((cev->tg.h.HHP&(one<<(pos-1)))!=0) && (cev->tg.det[pos].hge.FE>0))
	if((cev->tg.h.AHP&(1<<(pos-1)))!=0)
	  {
	    eAddBack = cev->tg.det[pos].addback.E/cal_par->tg.contr_e;
	    if(eAddBack<0 || eAddBack>S32K-10) eAddBack=S32K-10;
	    hist[sp][(int)(eAddBack)]++;
	  }
  free(cev);
  return SEPARATOR_DISCARD;
}
/*=======================================================================*/
int main(int argc, char *argv[])
{
  input_names_type *name;
  FILE *out,*gateNameFile;
  TFile *f;
  char aGateName[132],pGateName[132],det[132];
  int np,na,sum;
  
  if(argc!=6)
    {
      printf("Tigress_ECalAB_CsIPIDtr master_file_name idmin idmax chimin chimax\n");
      exit(-1);
    }
  
  printf("Program generates Tigress ECalAB histograms separated on CsIPIDtr.\n");
  
  memset(ch,0,sizeof(ch));
  ev=0;
  memset(id,0,sizeof(id));
  id[0][0]=1;//random
  id[1][0]=2;//1p
  id[2][0]=3;//2p
  id[0][1]=4;//1a
  id[1][1]=5;//1a1p
  id[0][2]=6;//2a
  //everything else goes to spectrum 0
  name=(input_names_type*)malloc(sizeof(input_names_type));
  memset(name,0,sizeof(input_names_type));
  
  cal_par=(calibration_parameters*)malloc(sizeof(calibration_parameters));
  memset(cal_par,0,sizeof(calibration_parameters));
  
  memset(hist,0,sizeof(hist));
  
  read_master(argv[1],name);
  
  idmin=atoi(argv[2]);
  idmax=atoi(argv[3]);
  chimin=atof(argv[4]);
  chimax=atof(argv[5]);
  
  if(name->flag.inp_data!=1)
    {
      printf("ERROR!!! Input data file not defined!\n");
      exit(EXIT_FAILURE);
    }
  
  if(name->flag.TIGRESS_cal_par==1)
    {
      printf("Tigress calibration read from %s.\n",name->fname.TIGRESS_cal_par);
      initialize_TIGRESS_calibration(&cal_par->tg,name->fname.TIGRESS_cal_par);
    }
  else
    {
      printf("ERROR!!! Tigress calibration parameters not defined!\n");
      exit(EXIT_FAILURE);
    }
  
  if(name->flag.root_gate_file==1)
    {
      printf("Using root gate file: %s\n",name->fname.root_gate_file);
      f=new TFile(name->fname.root_gate_file,"READ");
    }
  else
    {
      printf("ERROR!!! Root gate file not defined!\n");
      exit(-1);
    }
  
  if(name->flag.gate_name_file==1)
    {
      printf("Using gate name file: %s\n",name->fname.gate_name_file);
      gateNameFile=fopen(name->fname.gate_name_file,"r");
      
      while(fscanf(gateNameFile,"%s %s %s",det,aGateName,pGateName)!=EOF)
	{
	  pos=atoi(det);
	  if(strcmp(aGateName,"null"))
	    {
	      aGate[pos] = (TCutG *) gROOT->FindObject(aGateName);
	      aGateFlag[pos]=1;
	    }
	  if(strcmp(pGateName,"null"))
	    {
	      pGate[pos] = (TCutG *) gROOT->FindObject(pGateName);
	      pGateFlag[pos]=1;
	    }}
      fclose(gateNameFile);
    }
  else
    {
      printf("ERROR!!! Gate name file not defined!\n");
      exit(-1);
    }
  f->Close();
  
  for(pos=1;pos<NCSI;pos++)    
    if((pGateFlag[pos]+aGateFlag[pos])!=0)
      {
	printf("CsI %2d ",pos);
	if(aGateFlag[pos]==1)
	  printf("  alpha gate: %10s ",aGate[pos]->GetName());
	else
	  printf("%25s"," ");
	if(pGateFlag[pos]==1)
	  printf(" proton gate: %10s ",pGate[pos]->GetName());
	else
	  printf("%25s"," ");
	printf("\n");
      }
  
  sort(name); 

  printf("Number of PID events is     %10d\n",ev);
  printf(" Np/Na");
  for(np=0;np<pMax;np++)
    printf(" %8d ",np);
  printf("\n");
  for(na=0;na<aMax;na++)
    {
      printf("   %2d ",na);
      
      for(int np=0;np<pMax;np++)
	printf(" %8d ",ch[np][na]);
      printf("\n");
    }
  sum=0;
  for(np=0;np<pMax;np++)
    for(na=0;na<aMax;na++)
      sum+=ch[np][na];

  printf("Fraction of PID events is   %10.1f\n",100.*sum/ev);
  printf(" Np/Na");
  for(int np=0;np<pMax;np++)
    printf(" %8d ",np);
  printf("\n");
  for(int na=0;na<aMax;na++)
    {
      printf("   %2d ",na);
      
      for(int np=0;np<pMax;np++)
	printf(" %8.1f ",100.*ch[np][na]/ev);
      printf("\n");
    }
 
  if((out=fopen("ECal.mca","w"))==NULL)
    {
      printf("ERROR!!! I cannot open the mca file\n");
      exit(EXIT_FAILURE);
    }
	  
  for(np=0;np<NSP;np++)
    fwrite(hist[np],S32K*sizeof(int),1,out);
 
  fclose(out);


}
