#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
#include <global.h>
#include <vic_ensemble_io.h>
#include <mpi.h>

int vicNl_cell();
void resample();
void Comparision_Statistics();
void open_forcing_files();
void close_forcing_files();
void extract_cell();
void allocate_forcing();
void deallocate_forcing();
void initialize_atmos_BLUEWATERS(atmos_data_struct *, dmy_struct *,// FILE **,
                        forcing_cell_struct *forcing_cell,
                        soil_con_struct *);

/** Main Program **/

int main(int argc, char **argv)
{
  /* Initialize MPI */
  MPI_Init(&argc, &argv);
  int np,myid;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  /*Declare variables*/
  extern global_param_struct global_param;
  extern veg_lib_struct *veg_lib;
  int nvars = 7 ;
  int i,j;
  filenames_struct names;
  filep_struct filep;
  soil_con_struct soil_con; 
  veg_con_struct *veg_con;
  atmos_data_struct *atmos; 
  dmy_struct *dmy; 
  char MODEL_DONE;
  char RUN_MODEL;
  int cell_cnt;
  int Nveg_type;
  //Resampling variables
  double dt_rs;
  double dt_bs = 1; //hours
  int ipos; //initial position
  //Provide the forcing file dimensions
  grads_file_struct grads_file;
  grads_file.res = 1.0;
  grads_file.nx = 360;
  grads_file.ny = 180;
  grads_file.minlat = -89.50;
  grads_file.minlon = 0.50;
  grads_file.undef = -999.0;
  grads_file.nt = 32144;
  //Forcing variables
  forcing_filep_struct forcing_filep;
  forcing_name_struct forcing_name;
  forcing_cell_struct *forcing_cell;
  //Assign the forcing names
  strcpy(forcing_name.tair,"/scratch/02179/chaneyna/1.0deg/3hourly/tas/tas_2000-2010_cell.bin");
  strcpy(forcing_name.prec,"/scratch/02179/chaneyna/1.0deg/3hourly/prec/prec_2000-2010_cell.bin");
  strcpy(forcing_name.wind,"/scratch/02179/chaneyna/1.0deg/3hourly/wind/wind_2000-2010_cell.bin");
  strcpy(forcing_name.shum,"/scratch/02179/chaneyna/1.0deg/3hourly/shum/shum_2000-2010_cell.bin");
  strcpy(forcing_name.pres,"/scratch/02179/chaneyna/1.0deg/3hourly/pres/pres_2000-2010_cell.bin");
  strcpy(forcing_name.lwdown,"/scratch/02179/chaneyna/1.0deg/3hourly/dlwrf/dlwrf_2000-2010_cell.bin");
  strcpy(forcing_name.swdown,"/scratch/02179/chaneyna/1.0deg/3hourly/dswrf/dswrf_2000-2010_cell.bin");
  //Initialize global structure
  //strcpy(names.global,"/home/ice/nchaney/PROJECTS/BLUE_WATERS/TEST/DATA/global_param_4.1.2_MPI");
  //filep.globalparam = open_file(names.global,"r");
  //global_param = get_global_param(&names);//, filep.globalparam);
  //fclose(filep.globalparam);
  initialize_global();          
  //return;
  /** Obtain the global information for VIC **/
  global_param = get_global_param(&names); //Fill up the global parameter file

  /** Open up the forcing files **/
  open_forcing_files(&forcing_filep,&forcing_name);

  //MPI_Init(&argc, &argv);
  //int np,myid;
  //MPI_Status status;
  //MPI_Comm_size(MPI_COMM_WORLD,&np);
  //MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  /** Initialize the structures to hold the atmospheric data (before passing to VIC)**/
  int ncells = 1; //The number of cells to read in at once
  allocate_forcing(&forcing_cell,&grads_file,ncells);

  //Assign vegetation/soils names
  //strcpy(names.soil,"/share/home/02179/chaneyna/BLUE_WATERS/TEST/DATA/global_soils_3hourly_calib_smoothed.txt");
  strcpy(names.soil,"/share/home/02179/chaneyna/BLUE_WATERS/TEST/DATA/global_soils_default.txt");
  strcpy(names.veglib,"/share/home/02179/chaneyna/BLUE_WATERS/TEST/DATA/veglib.dat");
  strcpy(names.veg,"/share/home/02179/chaneyna/BLUE_WATERS/TEST/DATA/global_lai.txt");
  global_param.nrecs = grads_file.nt; //HARD CODED
  double *bs_output = (double *) malloc(sizeof(double)*global_param.nrecs*nvars);
  double *rs_output = (double *) malloc(sizeof(double)*global_param.nrecs*nvars);

  //Allocate memory for the atmospheric data
  alloc_atmos(global_param.nrecs, &atmos);

  /** Make Date Data Structure **/
  dmy = make_dmy(&global_param);

  /** Check and Open Files **/
  check_files(&filep, &names);

  /** Read Vegetation Library File **/
  veg_lib = read_veglib(filep.veglib,&Nveg_type);

  /** Set up all variables for the iteration through all cells **/
  int icell;
  int soil_ncells = 15836;
  MODEL_DONE=FALSE;
  RUN_MODEL=FALSE;
  cell_cnt = 0;
  int nrs=5;
  FILE *fp_metrics;
  //FILE *fp_output;
  double *rrmse = (double *) malloc(sizeof(double)*nvars);
  char str[100];
  sprintf (str,"/share/home/02179/chaneyna/BLUE_WATERS/TEST/DATA/Output/metrics_output_%d.txt", myid);
  fp_metrics = fopen(str,"w");
  /** Iterate for all cells in the soil file **/
  icell = myid;
  int linen = -1;
  while (icell < soil_ncells){
    for (i = 0; i < nvars; i++){
      rrmse[i] = 0.0;
    }
    /** Read the soil data **/
    //rewind(filep.soilparam);
    while (linen < icell){
      soil_con = read_soilparam(filep.soilparam, names.soil_dir, &cell_cnt, &RUN_MODEL, &MODEL_DONE);
      linen = linen + 1;
    }
    printf("%d %f %f\n",soil_con.gridcel,soil_con.lat,soil_con.lng);
    /** Read the vegetation parameters **/
    veg_con = read_vegparam(filep.vegparam,soil_con.gridcel,Nveg_type);
    calc_root_fractions(veg_con, &soil_con);
    /** Read in the forcing data for a single cell **/
    extract_cell(&forcing_filep,&grads_file,soil_con.lat,soil_con.lng,forcing_cell);
    /** Copy the forcing data to the VIC structure **/
    initialize_atmos_BLUEWATERS(atmos, dmy, forcing_cell,&soil_con);
    /** Run the VIC-ensemble framework **/
    /** Run a simulation to get all the global parameters the same THIS NEEDS TO CHANGE **/
    vicNl_cell(bs_output,soil_con,veg_con,dmy,atmos);
    //printf("Running the baseline simulation\n");
    vicNl_cell(bs_output,soil_con,veg_con,dmy,atmos);
    /** Resample the precipitation and run the new simulations **/
    fprintf(fp_metrics,"%d %f %f\n",soil_con.gridcel,soil_con.lat,soil_con.lng);
    for (i = 0; i < nrs; i++){
      /** Resample the precipitation **/
      dt_rs = (i+1)*dt_bs;
      //printf("Sampling Period %f\n",dt_rs);
      for (ipos = 0; ipos < (int)dt_rs/dt_bs; ipos++){
       // printf("Initial Position: %d\n",ipos);
        resample(atmos,global_param.nrecs,dt_rs,dt_bs,ipos);
        /** Run the model **/
        vicNl_cell(rs_output,soil_con,veg_con,dmy,atmos);
        /** Compare the output to the baseline **/
        Comparision_Statistics(bs_output,rs_output,global_param.nrecs,rrmse);
      }
      /** Compute the average rrmse for this sampling frequency **/
      fprintf(fp_metrics,"%f ",dt_rs);
      for (j = 0; j < nvars; j++){
        rrmse[j] = rrmse[j]/(dt_rs/dt_bs);
        fprintf(fp_metrics,"%f ",rrmse[j]);
      }
      fprintf(fp_metrics,"\n");
    }
    //Next cell
    icell = icell + np;
  }
  fclose(fp_metrics);
  //Free used memory//
  free_atmos(global_param.nrecs, &atmos);
  free_dmy(&dmy);
  free_veglib(&veg_lib);
  free_vegcon(&veg_con);
  //Close all the files
  fclose(filep.soilparam);
  fclose(filep.vegparam);
  fclose(filep.veglib);
  /** Close up the forcing files **/
  close_forcing_files(&forcing_filep);
  /** Deallocate the global forcing structure **/
  deallocate_forcing(&forcing_cell,&grads_file,ncells);
  /** Free memory **/
  free(rrmse);
  /** Close up MPI section **/
  MPI_Finalize();
  return 0;
}
