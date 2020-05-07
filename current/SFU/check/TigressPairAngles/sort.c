#include "sort.h"
/*=========================================================================*/
int getAngBin(double angleDeg){

  if((angleDeg > 181.)||(angleDeg < -0.1))
    return 99;

  int i;
  for(i=0;i<numAngBins;i++){
    if(i<99){
      if(fabs(angleDeg-angleList[i]) < 0.001){
        return i;
      }
    }else{
      return 99;
    }
  }
  //if we get here, then the bin doesn't exist yet
  if(numAngBins<99){
    //make a new bin
    angleList[numAngBins]=angleDeg;
    numAngBins++;
    return numAngBins-1;
  }else{
    return 99;
  }
}

double getAngDeg(int pos1, int pos2, int col1, int col2, int useAddBack){

  double pVec1[3], pVec2[3];

  if(useAddBack){
    //angles are between clovers, not cores
    pVec1[0] = (cal_par->tg.tpos_xyz[pos1][0][0] + cal_par->tg.tpos_xyz[pos1][2][0])/2.;
    pVec1[1] = (cal_par->tg.tpos_xyz[pos1][0][1] + cal_par->tg.tpos_xyz[pos1][2][1])/2.;
    pVec1[2] = (cal_par->tg.tpos_xyz[pos1][0][2] + cal_par->tg.tpos_xyz[pos1][2][2])/2.;
    pVec2[0] = (cal_par->tg.tpos_xyz[pos2][0][0] + cal_par->tg.tpos_xyz[pos2][2][0])/2.;
    pVec2[1] = (cal_par->tg.tpos_xyz[pos2][0][1] + cal_par->tg.tpos_xyz[pos2][2][1])/2.;
    pVec2[2] = (cal_par->tg.tpos_xyz[pos2][0][2] + cal_par->tg.tpos_xyz[pos2][2][2])/2.;
  }else{
    pVec1[0] = cal_par->tg.tpos_xyz[pos1][col1][0];
    pVec1[1] = cal_par->tg.tpos_xyz[pos1][col1][1];
    pVec1[2] = cal_par->tg.tpos_xyz[pos1][col1][2];
    pVec2[0] = cal_par->tg.tpos_xyz[pos2][col2][0];
    pVec2[1] = cal_par->tg.tpos_xyz[pos2][col2][1];
    pVec2[2] = cal_par->tg.tpos_xyz[pos2][col2][2];
  }

  double dp = pVec1[0]*pVec2[0];
  dp += pVec1[1]*pVec2[1];
  dp += pVec1[2]*pVec2[2];
  double mag1 = pVec1[0]*pVec1[0];
  mag1 += pVec1[1]*pVec1[1];
  mag1 += pVec1[2]*pVec1[2];
  mag1 = sqrt(mag1);
  double mag2 = pVec2[0]*pVec2[0];
  mag2 += pVec2[1]*pVec2[1];
  mag2 += pVec2[2]*pVec2[2];
  mag2 = sqrt(mag2);
  
  double ang = acos(dp/(mag1*mag2)); //from 0 to pi radians
  ang = ang*360.0/TWOPI;
  
  return ang;
}

/*=========================================================================*/
int main(int argc, char *argv[])
{
  input_names_type* name;
  char inp[256];
  int useAddBack = 0;
  
  if(argc!=3)
    {
      printf("check_TigressPairAngles master_file_name useAddBack\n");
      printf("Specifies the number of unique crystal pairs at a given angle in TIGRESS.  You will be prompted to specify any missing positions/cores in the array.\n");
      exit(-1);
    }
  
  numAngBins = 0;
  memset(ignoredPositions,0,sizeof(ignoredPositions));

  if(atoi(argv[2])==1){
    useAddBack = 1;
  }

  int posIgn,coreIgn;
  while(true){
    printf("Enter any positions and cores to ignore (position and core number separated by a space, or just a position number, or -1 if done):\n");
    scanf("%[^\n]%*c",inp); //get string with spaces, until newline
    //printf("input: %s\n",inp);
    if((sscanf(inp,"%i %i",&posIgn,&coreIgn)) == 2){
      if(((posIgn>0)&&(posIgn<NPOSTIGR))&&((coreIgn>=0)&&(coreIgn<NCOL))){
        if(useAddBack){
          printf("Ignoring position %i.\n",posIgn);
          for(col=0;col<NCOL;col++){
            ignoredPositions[posIgn][col]=1;
          }
        }else{
          printf("Ignoring position %i, core %i.\n",posIgn,coreIgn);
          ignoredPositions[posIgn][coreIgn]=1;
        }
      }else{
        printf("Invalid position and/or core values.  Try again:\n");
      }
    }else if((sscanf(inp,"%i",&posIgn)) == 1){
      if(posIgn < 0){
        break;
      }else if((posIgn>0)&&(posIgn<NPOSTIGR)){
        printf("Ignoring position %i.\n",posIgn);
        for(col=0;col<NCOL;col++){
          ignoredPositions[posIgn][col]=1;
        }
      }else{
        printf("Invalid position value.  Try again:\n");
      }
    }else{
      printf("Didn't understand input.  Try again:\n");
    }
  }
  
  name=(input_names_type*)malloc(sizeof(input_names_type));
  memset(name,0,sizeof(input_names_type));
  
  cal_par=(calibration_parameters*)malloc(sizeof(calibration_parameters));
  memset(cal_par,0,sizeof(calibration_parameters));

  read_master(argv[1],name);
  
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

  if(useAddBack){
    for(pos=1;pos<NPOSTIGR;pos++){
      if(ignoredPositions[pos][0]==0){
        for(pos2=pos+1;pos2<NPOSTIGR;pos2++){
          if(ignoredPositions[pos2][0]==0){
            double corrAng = getAngDeg(pos, pos2, 0, 1,useAddBack);
            numPairs[getAngBin(corrAng)]++;
          }
        }
      }
    }
  }else{
    for(pos=1;pos<NPOSTIGR;pos++){
      for(col=0;col<NCOL;col++){
        if(ignoredPositions[pos][col]==0){
          for(pos2=pos+1;pos2<NPOSTIGR;pos2++){
            for(col2=0;col2<NCOL;col2++){
              if(ignoredPositions[pos2][col2]==0){
                double corrAng = getAngDeg(pos, pos2, col, col2,useAddBack);
                numPairs[getAngBin(corrAng)]++;
              }
            }
          }
          //count pairs from the same position/different core, can only do this if no addback
          for(col2=col+1;col2<NCOL;col2++){
            if(ignoredPositions[pos][col2]==0){
              double corrAng = getAngDeg(pos, pos, col, col2,useAddBack);
              numPairs[getAngBin(corrAng)]++;
            }
          }  
        }
      }
    }
  }

  

  int i;
  printf("Angle      Number of Unique Pairs\n");
  for(i=0;i<numAngBins;i++){
    if(i<100){
      printf("%6.2f  %10i\n",angleList[i],numPairs[i]);
    }
  }
  
}
