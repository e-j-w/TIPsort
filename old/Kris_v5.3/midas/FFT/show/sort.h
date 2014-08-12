#ifndef SORT_H
#define SORT_H

#include "SFU-common.h"
#include "sort_but_not_assemble.h"
#include "TH1D.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TApplication.h"


TH1D *h;
TH1 *fourier_transform;
//TFile fi("30MHzrefFFTnocage.root", "new"); //use if one wants to write histogram to a file
TCanvas *c;  
TApplication *theApp;
int  chn;

#endif
