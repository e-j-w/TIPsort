#ifndef SORT_H
#define SORT_H

#include "SFU-common.h"
#include "SFU-format.h"
#include "SFU-decoder.h"
#include "SFU-cal.h"

#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TApplication.h"

int pos,col;
double supLow,supHigh,rf_period;

TH2D *h, *hrf, *hrad;

calibration_parameters* cal_par;
#endif
