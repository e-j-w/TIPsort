#ifndef SORT_H
#define SORT_H

//Prior to putting in 1ns units!!!
#define RFphase 135.74544

#include "SFU-common.h"
#include "SFU-format.h"
#include "SFU-decoder.h"

#include "TH2D.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TApplication.h"

int pos,col,sup;

TH2D *hFITlocal,*hCFD;
TH1D *h_tFITlocal,*h_tCFD,*h_trf,*h_tdiffFITlocal,*h_tdiffCFD;

double chimin,chimax;

#endif
