#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

void resample(atmos_data_struct *atmos, int n,double dt_rs, double dt_bs,int ipos){

  /** Declare local variables **/
  //double dt_rs = 72.00; resampled time step
  //double dt_bs = 3.00; baseline time step
  //int ipos; initial position
  double data_rs[n];
  int i,lpos,rpos,ncount;
  int nt = (int)dt_rs/dt_bs; //Number of baseline time steps in a resampled one
  /** Fill up the temporary array **/
  for (i=0;i<n;i++){
    data_rs[i] = *atmos[i].prec;
  }
  /** Resample the precipitation **/
  lpos = ipos;
  rpos = lpos + nt;
  ncount = 0;
  for (i=0;i<n;i++){
    ncount = ncount + 1;
    //Determine what happens depending on the distance
    if(abs(rpos - i) > abs(i - lpos)){
        data_rs[i] = *atmos[lpos].prec;
    }
    else if(abs(rpos - i) < abs(i - lpos)){
        data_rs[i] = *atmos[rpos].prec;
    }
    else{
        data_rs[i] = (*atmos[lpos].prec + *atmos[rpos].prec)/2;
    }
    //Update lpos and rpos after reaching nt
    if (ncount == nt){
        lpos = lpos + nt;
        rpos = rpos + nt;
        ncount = 0;
    }
    //If the rpos is out of bounds quit the do loop
    if (rpos >= n)break;
  }
  double sum0,sum1;
  sum0 = 0.0;
  sum1 = 0.0;
  for (i=0;i<n;i++){
    sum0 = sum0 + data_rs[i];
    sum1 = sum1 + *atmos[i].prec;
  }
  /** Set the original precipitation to the resampled one **/
  for (i=0;i<n;i++){
    *atmos[i].prec = data_rs[i];
  }
}
