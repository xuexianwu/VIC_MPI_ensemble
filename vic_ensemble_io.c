#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
#include <vic_ensemble_io.h>
#include <netcdf.h>

void read_forcing();
void downscale_data();

void open_forcing_files(forcing_filep_struct *forcing_filep, forcing_name_struct *forcing_name){
  /** Open up the forcing files **/
  forcing_filep->tair = fopen(forcing_name->tair,"r"); //Air temperature
  forcing_filep->prec = fopen(forcing_name->prec,"r"); //Precipitation
  forcing_filep->wind = fopen(forcing_name->wind,"r"); //Wind speed
  forcing_filep->shum = fopen(forcing_name->shum,"r"); //Specific Humidity
  forcing_filep->pres = fopen(forcing_name->pres,"r"); //Pressure
  forcing_filep->lwdown = fopen(forcing_name->lwdown,"r"); //Downward Longwave Radiation
  forcing_filep->swdown = fopen(forcing_name->swdown,"r"); //Downward Shortwave Radiation
  //printf("Opened all the forcing files\n");
}

int open_netcdf_forcing_files(char netcdf_forcing_filename[MAXSTRING]){
  int ncid;
  /** Open netcdf forcing file **/
  nc_open(netcdf_forcing_filename,NC_NOWRITE,&ncid);
  return ncid;
}

void close_forcing_files(forcing_filep_struct *forcing_filep){
  /** Close the forcing files **/
  fclose(forcing_filep->tair); //Air temperature
  fclose(forcing_filep->prec); //Precipitation
  fclose(forcing_filep->wind); //Wind speed
  fclose(forcing_filep->shum); //Specific Humidity
  fclose(forcing_filep->pres); //Pressure
  fclose(forcing_filep->lwdown); //Downward Longwave Radiation
  fclose(forcing_filep->swdown); //Downward Shortwave Radiation
  //printf("Closed all the forcing files\n");
}

void extract_cell_netcdf(int ncid, grads_file_struct *grads_file, forcing_cell_struct *forcing_cell, int cell_id){

 int ncells = 15836;
 int nt = grads_file->nt_netcdf;//93504;
 int varid;
 int cell_ids[ncells];
 int icell;
 size_t count[2],start[2];
 count[0] = nt;
 count[1] = 1;
 start[0] = 0;
 float data_prec[nt];
 float data_pres[nt];
 float data_wind[nt];
 float data_swdown[nt];
 float data_lwdown[nt];
 float data_tair[nt];
 float data_shum[nt];

 //Read in the cell id array
 nc_inq_varid(ncid,"id",&varid);
 nc_get_var_int(ncid,varid,&cell_ids[0]);

 //Figure out the cells placement in the ids array
 for (icell = 0; icell < ncells; icell++){
  if (cell_ids[icell] == cell_id){break;}
 }
 start[1] = cell_id;

 //Extract the data per variable
  /** Precipitation **/
  nc_inq_varid(ncid,"prec",&varid);
  nc_get_vara_float(ncid,varid,start,count,&data_prec[0]);
  for (int i = 0; i < nt; i++){
   printf("%f\n",data_prec[i]);
  }
  //read_forcing(grads_file,forcing_filep->prec,&data_prec,i,j,grads_file->nt_prec);
  //downscale_data(grads_file->nt,grads_file->nt_prec,data_prec,&data,1);
  //for (t = 0; t < grads_file->nt; t++){forcing_cell[0].prec[t] = data[t];}
}

void extract_cell(forcing_filep_struct *forcing_filep, grads_file_struct *grads_file,
                  double lat, double lon, forcing_cell_struct *forcing_cell){
  /** Function to extract the data for a given variable given the file dimensions and a lat/lon **/
  /** Determine the i,j on the grid **/
  int i = (int)round((lat - grads_file->minlat)/grads_file->res) + 1;
  int j = (int)round((lon - grads_file->minlon)/grads_file->res) + 1;

  /** Iterate through each time step, extracting the desired data **/
  int t;
  double data_prec[grads_file->nt_prec];
  double data_pres[grads_file->nt_pres];
  double data_wind[grads_file->nt_wind];
  double data_swdown[grads_file->nt_swdown];
  double data_lwdown[grads_file->nt_lwdown];
  double data_tair[grads_file->nt_tair];
  double data_shum[grads_file->nt_shum];
  double data[grads_file->nt];
  /** Extract data for variable **/
  /** Precipitation **/
  //printf("Reading in the Precipitation\n");
  //printf("%d\n",grads_file->nt_prec); 
  read_forcing(grads_file,forcing_filep->prec,&data_prec,i,j,grads_file->nt_prec); 
  //for (t = 0; t < grads_file->nt_prec; t++){if(t<100){printf("%d %d %f\n",i,j,data_prec[t]);}}
  downscale_data(grads_file->nt,grads_file->nt_prec,data_prec,&data,1);
  for (t = 0; t < grads_file->nt; t++){forcing_cell[0].prec[t] = data[t];}

  /** Pressure **/
  //printf("Reading in the Pressure\n");
  read_forcing(grads_file,forcing_filep->pres,&data_pres,i,j,grads_file->nt_pres);
  downscale_data(grads_file->nt,grads_file->nt_pres,data_pres,&data,0);
  for (t = 0; t < grads_file->nt; t++){forcing_cell[0].pres[t] = data[t];}

  /** Wind Speed **/
 // printf("Reading in the Wind Speed\n");
  read_forcing(grads_file,forcing_filep->wind,&data_wind,i,j,grads_file->nt_wind);
  downscale_data(grads_file->nt,grads_file->nt_wind,data_wind,&data,0);
  for (t = 0; t < grads_file->nt; t++){forcing_cell[0].wind[t] = data[t];}

  /** Downward Shortwave Radiation **/
  //printf("Reading in the Shortwave Radiation\n");
  read_forcing(grads_file,forcing_filep->swdown,&data_swdown,i,j,grads_file->nt_swdown); 
  downscale_data(grads_file->nt,grads_file->nt_swdown,data_swdown,&data,0);
  for (t = 0; t < grads_file->nt; t++){forcing_cell[0].swdown[t] = data[t];}

  /** Downward Longwave Radiation **/
  //printf("Reading in the Longwave Radiation\n");
  read_forcing(grads_file,forcing_filep->lwdown,&data_lwdown,i,j,grads_file->nt_lwdown);
  downscale_data(grads_file->nt,grads_file->nt_lwdown,data_lwdown,&data,0);
  for (t = 0; t < grads_file->nt; t++){forcing_cell[0].lwdown[t] = data[t];}

  /** Air temperature **/
  //printf("Reading in the Air Temperature\n");
  read_forcing(grads_file,forcing_filep->tair,&data_tair,i,j,grads_file->nt_tair);
  downscale_data(grads_file->nt,grads_file->nt_tair,data_tair,&data,0);
  for (t = 0; t < grads_file->nt; t++){forcing_cell[0].tair[t] = data[t];}

  /** Specific Humidity **/
  //printf("Reading in the Specific Humidity\n");
  read_forcing(grads_file,forcing_filep->shum,&data_shum,i,j,grads_file->nt_shum);
  downscale_data(grads_file->nt,grads_file->nt_shum,data_shum,&data,0);
  for (t = 0; t < grads_file->nt; t++){forcing_cell[0].shum[t] = data[t];}

  /** Convert data to vic format **/
  for (t = 0; t < grads_file->nt; t++){
    /** tair **/
    forcing_cell[0].tair[t] = forcing_cell[0].tair[t] - 273.15; //Kelvin to Celsius
    forcing_cell[0].pres[t] = forcing_cell[0].pres[t]/1000; //Pa to kPa
    //printf("%f %f %f %f %f %f %f\n",forcing_cell[0].prec[t],forcing_cell[0].pres[t],
    //  forcing_cell[0].wind[t],forcing_cell[0].swdown[t],forcing_cell[0].lwdown[t],
    //  forcing_cell[0].tair[t],forcing_cell[0].shum[t]);
  }
}

/*
void read_forcing(grads_file_struct *grads_file,FILE *filep,double *data,int i, int j){
  int pos,ipos,t;
  float var;
  ipos = (i-1)*grads_file->nx + j-1;
  // Seek to the initial position //
  fseek(filep,sizeof(var)*ipos,SEEK_SET);
  for (t = 0; t < grads_file->nt; t++){
    fread(&var,sizeof(var),1,filep);
    data[t]  = (double)var;
    //Update position
    pos = ipos + grads_file->nx*grads_file->ny*t;
    fseek(filep,sizeof(var)*(grads_file->nx*grads_file->ny-1),SEEK_CUR);
  }
}
*/
void read_forcing(grads_file_struct *grads_file,FILE *filep,double *data,int i, int j, int nt){
  unsigned long long ipos;
  int t;
  float var;
  ipos = (unsigned long long)grads_file->nx*(unsigned long long)nt*((unsigned long long)i-1);
  ipos = (unsigned long long)nt*((unsigned long long)j-1) + ipos;
  //printf("%llu\n",ipos);
  // Seek to the initial position of the variable//
  fseek(filep,sizeof(var)*ipos,SEEK_SET);
  for (t = 0; t < nt; t++){
  //for (t=0; t < 10; t++){
    fread(&var,sizeof(var),1,filep);
    data[t]  = (double)var;
    //printf("%lf\n",data[t]);
    //Update position
    //pos = ipos + 1;
    //fseek(filep,sizeof(var)*(grads_file->nx*grads_file->ny-1),SEEK_CUR);
  }
}

/** Downscale in time **/
void downscale_data(int nt, int nt_orig ,double *data_us, double *data_ds, int flag_average){
  int t,tsum,torig,dt;
  dt = nt/nt_orig;
  tsum = 0;
  torig = 0;
  for (t=0;t<nt;t++){
    if (tsum == dt){tsum = 0;torig = torig + 1;}
    tsum = tsum + 1;
    if (flag_average == 0){
      data_ds[t] = data_us[torig];
      }
    else{
      data_ds[t] = data_us[torig]/(double)dt;
      }
  }
} 

void allocate_forcing(forcing_cell_struct **forcing_cell, grads_file_struct *grads_file, int ncells){
  int i;
  *forcing_cell = (forcing_cell_struct *) calloc(ncells, sizeof(forcing_cell_struct));
  for (i = 0; i < ncells; i++) {
    (*forcing_cell)[i].prec = (double *) calloc(grads_file->nt, sizeof(double));
    (*forcing_cell)[i].tair = (double *) calloc(grads_file->nt, sizeof(double));
    (*forcing_cell)[i].shum = (double *) calloc(grads_file->nt, sizeof(double));
    (*forcing_cell)[i].wind = (double *) calloc(grads_file->nt, sizeof(double));
    (*forcing_cell)[i].lwdown = (double *) calloc(grads_file->nt, sizeof(double));
    (*forcing_cell)[i].swdown = (double *) calloc(grads_file->nt, sizeof(double));
    (*forcing_cell)[i].pres = (double *) calloc(grads_file->nt, sizeof(double));
  }
}

void deallocate_forcing(forcing_cell_struct **forcing_cell, grads_file_struct *grads_file, int ncells){
  int i;
  //Deallocate the memory used for storing the forcing data
  for (i = 0; i < ncells; i++) {
    free((*forcing_cell)[i].tair);
  }
  free(*forcing_cell);
}
