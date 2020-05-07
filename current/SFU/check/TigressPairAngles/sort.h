#ifndef SORT_H
#define SORT_H

#include "SFU-common.h"
#include "SFU-format.h"
#include "SFU-decoder.h"
#include "SFU-cal.h"

#include "TCutG.h"
#include "TFile.h"
#include "TROOT.h"

#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TApplication.h"

#define MAXNPART 10
#define amuMeV 931.494061

int pos,col,pos2,col2;

calibration_parameters* cal_par;

//quantities for calculation of doppler shift
char str[256];

//angular bins
double angleList[100];
int numPairs[100];
int numAngBins;

int ignoredPositions[NPOSTIGR][NCOL];

#endif
