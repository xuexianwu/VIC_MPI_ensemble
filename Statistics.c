#include <math.h>

void Comparision_Statistics(double *bs_output,double *rs_output,int n,
			   double *rrmse){
  /** Calculate RRMSE **/
  int nvars = 7;
  double rmse[nvars],bs_mean[nvars],rs_mean[nvars];
  double min[nvars],max[nvars];
  int i,j,pos;
  /** Initialize all the variables to zero **/
  for (j = 0; j < nvars; j++){
    rmse[j] = 0.0;
    bs_mean[j] = 0.0;
    rs_mean[j] = 0.0;
  }
  /** Sum up all the squared errors **/
  for (i = 0; i < n; i++){
    for (j = 0; j < nvars ; j++){          
      pos = i*nvars + j;
      rmse[j] = rmse[j] + pow(rs_output[pos]-bs_output[pos],2);
      bs_mean[j] = bs_mean[j] + bs_output[pos];
      rs_mean[j] = rs_mean[j] + rs_output[pos]; 
    }
  }
  /** Calculate the maximum and minimum values **/
  for (j = 0; j < nvars; j++){
    min[j] = 1000000.0;
    max[j] = -1000000.0;
  }
  for (i = 0; i < n; i++){
    for (j = 0; j < nvars; j++){          
      pos = i*nvars + j;
      if (min[j] > bs_output[pos])min[j] = bs_output[pos];
      if (max[j] < bs_output[pos])max[j] = bs_output[pos];
    }
  }
  /** Calculate the relative root mean squared error**/
  for (j = 0; j < nvars ; j++){
    bs_mean[j] = bs_mean[j]/n;
    rs_mean[j] = rs_mean[j]/n;
    rmse[j] = sqrt(rmse[j]/n);
    if (max[j]-min[j] == 0){rrmse[j] = rrmse[j];}
    else{rrmse[j] = rrmse[j] + rmse[j]/(max[j]-min[j]);}
    //printf("bs_mean: %f rs_mean: %f min: %f max: %f rmse: %f rrmse: %f\n", 
    //       bs_mean[j],rs_mean[j],min[j],max[j],rmse[j],rrmse[j]);
    //printf("rrmse: %f ",rrmse[j]);
  }
  //printf("\n");

}
