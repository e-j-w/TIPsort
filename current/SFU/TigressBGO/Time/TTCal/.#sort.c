#ifndef SORT_H
#define SORT_H


#include "SFU-common.h"
#include "SFU-format.h"
#include "SFU-decoder.h"
#include "SFU-cal.h"

#include "TH2D.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TApplication.h"

int  hist[NRING*2][S32K];
calibration_parameters* cal_par;
double supLow,supHigh;
int pos1,pos2,col1,col2,ring1,ring2;

TH2D *h;

#endif
