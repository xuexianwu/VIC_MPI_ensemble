typedef struct {
  FILE *tair;     /* Air temperature (C) */
  FILE *prec;    /* Precipitation (mm) */
  FILE *wind;     /* Wind Speed (m/s) */
  FILE *shum;       /* Specific Humidity */
  FILE *pres;      /* Pressure (kPa) */
  FILE *lwdown;      /* Downward Longwave Radiation (W/m2) */
  FILE *swdown;         /* Downward Shortwave Radiation (W/m2) */
} forcing_filep_struct;

typedef struct {
  char tair[MAXSTRING];      /* Air temperature (C) */
  char prec[MAXSTRING];      /* Precipitation (mm) */
  char wind[MAXSTRING];      /* Wind Speed (m/s) */
  char shum[MAXSTRING];        /* Specific Humidity */
  char pres[MAXSTRING];      /* Pressure (kPa) */
  char lwdown[MAXSTRING];    /* Downward Longwave Radiation (W/m2) */
  char swdown[MAXSTRING];    /* Downward Shortwave Radiation (W/m2) */
} forcing_name_struct;

/**typedef struct{
  int nx;
  int ny;
  double res;
  double undef;
  double minlat;
  double minlon;
  int nt;
  int year;
  int month;
  int day;
  int hour;
  int nt_prec;
  int nt_tair;
  int nt_wind;
  int nt_shum;
  int nt_pres;
  int nt_lwdown;
  int nt_swdown;
  int nt_netcdf;
  int dt;
} grads_file_struct;**/

typedef struct{
  double *tair;
  double *prec;
  double *wind;
  double *shum;
  double *pres;
  double *lwdown;
  double *swdown;
} forcing_cell_struct;
