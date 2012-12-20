#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
#include <vic_ensemble_io.h>

void read_forcing();

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

void extract_cell(forcing_filep_struct *forcing_filep, grads_file_struct *grads_file,
                  double lat, double lon, forcing_cell_struct *forcing_cell){
  /** Function to extract the data for a given variable given the file dimensions and a lat/lon **/
  /** Determine the i,j on the grid **/
  int i = (int)round((lat - grads_file->minlat)/grads_file->res) + 1;
  int j = (int)round((lon - grads_file->minlon)/grads_file->res) + 1;

  /** Iterate through each time step, extracting the desired data **/
  int t;
  double data[grads_file->nt];
  /** Extract data for variable **/

  /** Precipitation **/
  //printf("Reading in the Precipitation\n");
  read_forcing(grads_file,forcing_filep->prec,&data,i,j); 
  for (t = 0; t < grads_file->nt; t++){forcing_cell[0].prec[t] = data[t];}

  /** Pressure **/
  //printf("Reading in the Pressure\n");
  read_forcing(grads_file,forcing_filep->pres,&data,i,j);
  for (t = 0; t < grads_file->nt; t++){forcing_cell[0].pres[t] = data[t];}

  /** Wind Speed **/
 // printf("Reading in the Wind Speed\n");
  read_forcing(grads_file,forcing_filep->wind,&data,i,j);
  for (t = 0; t < grads_file->nt; t++){forcing_cell[0].wind[t] = data[t];}

  /** Downward Shortwave Radiation **/
  //printf("Reading in the Shortwave Radiation\n");
  read_forcing(grads_file,forcing_filep->swdown,&data,i,j); 
  for (t = 0; t < grads_file->nt; t++){forcing_cell[0].swdown[t] = data[t];}

  /** Downward Longwave Radiation **/
  //printf("Reading in the Longwave Radiation\n");
  read_forcing(grads_file,forcing_filep->lwdown,&data,i,j);
  for (t = 0; t < grads_file->nt; t++){forcing_cell[0].lwdown[t] = data[t];}

  /** Air temperature **/
  //printf("Reading in the Air Temperature\n");
  read_forcing(grads_file,forcing_filep->tair,&data,i,j);
  for (t = 0; t < grads_file->nt; t++){forcing_cell[0].tair[t] = data[t];}

  /** Specific Humidity **/
  //printf("Reading in the Specific Humidity\n");
  read_forcing(grads_file,forcing_filep->shum,&data,i,j);
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
void read_forcing(grads_file_struct *grads_file,FILE *filep,double *data,int i, int j){
  unsigned long ipos;
  int t;
  float var;
  ipos = grads_file->nt*(j-1) + grads_file->nx*grads_file->nt*(i-1);
  // Seek to the initial position of the variable//
  fseek(filep,sizeof(var)*ipos,SEEK_SET);
  for (t = 0; t < grads_file->nt; t++){
  //for (t=0; t < 10; t++){
    fread(&var,sizeof(var),1,filep);
    data[t]  = (double)var;
    //printf("%lf\n",data[t]);
    //Update position
    //pos = ipos + 1;
    //fseek(filep,sizeof(var)*(grads_file->nx*grads_file->ny-1),SEEK_CUR);
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
