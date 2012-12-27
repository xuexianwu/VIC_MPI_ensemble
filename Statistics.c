#include <math.h>
#include <time.h>

void Comparision_Statistics(double *bs_output,double *rs_output,int n,
			   double *rrmse,int year, int month, int day, int hour,
			   double dt, int current_month){  
  time_t t;
  struct tm gtime;
  struct tm gtime_original;
  gtime.tm_sec = 0;
  gtime.tm_min = 0;
  gtime.tm_hour = hour;
  gtime.tm_mday = day;
  gtime.tm_mon = month - 1;
  gtime.tm_year = year - 1900;
  gtime_original = gtime;

  /* Convert it to local time representation. */
  int i;
  //for (i = 0; i < n; i++){
  //  printf("%d %d %d\n",gtime.tm_year+1900,gtime.tm_mon+1,gtime.tm_mday);
  //  t = t + 3600*dt;
  //  gmtime_r(&t,&gtime);
  //}
  /** Calculate RRMSE **/
  int nvars = 7;
  double rmse[nvars],bs_mean[nvars],rs_mean[nvars];
  double min[nvars],max[nvars];
  int j,pos;
  /** Initialize all the variables to zero **/
  for (j = 0; j < nvars; j++){
    rmse[j] = 0.0;
    bs_mean[j] = 0.0;
    rs_mean[j] = 0.0;
  }
  /** Sum up all the squared errors **/
  t = mktime(&gtime_original);
  for (i = 0; i < n; i++){
    t = t + 3600*dt;
    gmtime_r(&t,&gtime);
    if (gtime.tm_mon == current_month){
      for (j = 0; j < nvars ; j++){          
        pos = i*nvars + j;
        rmse[j] = rmse[j] + pow(rs_output[pos]-bs_output[pos],2);
        bs_mean[j] = bs_mean[j] + bs_output[pos];
        rs_mean[j] = rs_mean[j] + rs_output[pos]; 
//        printf("%d %d %d\n",gtime.tm_mon,t,current_month);
//        printf("%f %f %f\n",rmse[j],bs_mean[j],rs_mean[j]);
      }
    }
  }
  /** Calculate the maximum and minimum values **/
  for (j = 0; j < nvars; j++){
    min[j] = 1000000.0;
    max[j] = -1000000.0;
  }
  t = mktime(&gtime_original);
  for (i = 0; i < n; i++){
    t = t + 3600*dt;
    gmtime_r(&t,&gtime);
    if (gtime.tm_mon == current_month){
      for (j = 0; j < nvars; j++){          
        pos = i*nvars + j;
        if (min[j] > bs_output[pos])min[j] = bs_output[pos];
        if (max[j] < bs_output[pos])max[j] = bs_output[pos];
      }
    }
  }
  /** Calculate the relative root mean squared error**/
  for (j = 0; j < nvars ; j++){
    bs_mean[j] = bs_mean[j]/n;
    rs_mean[j] = rs_mean[j]/n;
    rmse[j] = sqrt(rmse[j]/n);
    if (max[j]-min[j] == 0){rrmse[j] = rrmse[j];}
    else{rrmse[j] = rrmse[j] + rmse[j]/(max[j]-min[j]);}
  }
}
