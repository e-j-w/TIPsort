#ifndef SORT_H
#define SORT_H

#define NRING 7

#include "SFU-common.h"
#include "SFU-format.h"
#include "SFU-decoder.h"
#include "SFU-cal.h"

calibration_parameters* cal_par;
int hist[2*NRING][S32K];
int pos,col,pin,ring;

#endif
