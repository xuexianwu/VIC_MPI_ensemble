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
  int nvars = 7; //Number of forcing files
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
  double dt_bs = 3; //hours
  int nrs;//5; Number of resampling intervals
  int soil_ncells;
  int ipos; //initial position
  grads_file_struct grads_file;

  //Open the global parameter file
  char global_filename[MAXSTRING];
  strcpy(global_filename,argv[1]);
  forcing_name_struct forcing_name;
  FILE *global_fp;
  global_fp=fopen(global_filename,"r");
  char tmp_s[MAXSTRING];
  float tmp_f;
  int tmp_i;
  char metrics_filename[MAXSTRING],metrics_root[MAXSTRING];
  grads_file.year = 2000;
  grads_file.month = 1;
  grads_file.day = 1;
  grads_file.hour = 0;

  //Grads file info
  fscanf(global_fp,"%s %f",tmp_s,&tmp_f); grads_file.res = tmp_f;
  //printf("%s %f\n",tmp_s,grads_file.res);
  fscanf(global_fp,"%s %d",tmp_s,&tmp_i); grads_file.nx = tmp_i;
  //printf("%s %d\n",tmp_s,grads_file.nx);
  fscanf(global_fp,"%s %d",tmp_s,&tmp_i); grads_file.ny = tmp_i;
  //printf("%s %d\n",tmp_s,grads_file.ny);
  fscanf(global_fp,"%s %f",tmp_s,&tmp_f); grads_file.minlat = tmp_f;
  //printf("%s %f\n",tmp_s,grads_file.minlat);
  fscanf(global_fp,"%s %f",tmp_s,&tmp_f); grads_file.minlon = tmp_f;
  //printf("%s %f\n",tmp_s,grads_file.minlon);
  fscanf(global_fp,"%s %f",tmp_s,&tmp_f); grads_file.undef = tmp_f;
  //printf("%s %f\n",tmp_s,grads_file.undef);
  fscanf(global_fp,"%s %d",tmp_s,&tmp_i); grads_file.nt = tmp_i;
  //printf("%s %d\n",tmp_s,grads_file.nt);
  //Meteorological input
  fscanf(global_fp,"%s %s",tmp_s,&forcing_name.tair);
  //printf("%s %s\n",tmp_s,forcing_name.tair);
  fscanf(global_fp,"%s %s",tmp_s,&forcing_name.prec);
  //printf("%s %s\n",tmp_s,forcing_name.prec);
  fscanf(global_fp,"%s %s",tmp_s,&forcing_name.wind);
  //printf("%s %s\n",tmp_s,forcing_name.wind);
  fscanf(global_fp,"%s %s",tmp_s,&forcing_name.shum);
  //printf("%s %s\n",tmp_s,forcing_name.shum);
  fscanf(global_fp,"%s %s",tmp_s,&forcing_name.pres);
  //printf("%s %s\n",tmp_s,forcing_name.pres);
  fscanf(global_fp,"%s %s",tmp_s,&forcing_name.lwdown);
  //printf("%s %s\n",tmp_s,forcing_name.lwdown);
  fscanf(global_fp,"%s %s",tmp_s,&forcing_name.swdown);
  //printf("%s %s\n",tmp_s,forcing_name.swdown);
  //Land Info
  fscanf(global_fp,"%s %s",tmp_s,&names.soil);
  //printf("%s %s\n",tmp_s,names.soil);
  fscanf(global_fp,"%s %s",tmp_s,&names.veglib);
  //printf("%s %s\n",tmp_s,names.veglib);
  fscanf(global_fp,"%s %s",tmp_s,&names.veg);
  //printf("%s %s\n",tmp_s,names.veg);
  //Output Metrics
  fscanf(global_fp,"%s %s",tmp_s,&metrics_root);
  //printf("%s %s\n",tmp_s,metrics_root);
  sprintf (metrics_filename,"%s/metrics_output_%d.txt",metrics_root,myid);
  //Number of cells
  fscanf(global_fp,"%s %d",tmp_s,&soil_ncells);
  //Number of resampling intervals
  fscanf(global_fp,"%s %d",tmp_s,&nrs);

  //Forcing variables
  forcing_filep_struct forcing_filep;
  forcing_cell_struct *forcing_cell;
  //Assign the forcing names
  //Close the global parameter file
  fclose(global_fp);

  //Initialize global structure;
  initialize_global();          
  //return;
  /** Obtain the global information for VIC **/
  global_param = get_global_param(&names); //Fill up the global parameter file

  /** Open up the forcing files **/
  open_forcing_files(&forcing_filep,&forcing_name);

  /** Initialize the structures to hold the atmospheric data (before passing to VIC)**/
  int ncells = 1; //The number of cells to read in at once
  allocate_forcing(&forcing_cell,&grads_file,ncells);

  //Assign vegetation/soils names
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
  //int soil_ncells = 2;//15836;
  MODEL_DONE=FALSE;
  RUN_MODEL=FALSE;
  cell_cnt = 0;
  FILE *fp_metrics;
  //FILE *fp_output;
  double *rrmse = (double *) malloc(sizeof(double)*nvars);
  double rrmse_month[12][nvars];

  fp_metrics = fopen(metrics_filename,"w");
  /** Iterate for all cells in the soil file **/
  icell = myid;
  int linen = -1;
  int current_month;
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
        for (current_month = 0; current_month < 12; current_month++){
          /** Compare the output to the baseline **/
          Comparision_Statistics(bs_output,rs_output,global_param.nrecs,rrmse,grads_file.year,
			       grads_file.month,grads_file.day,grads_file.hour,dt_bs,current_month);
          for (j = 0; j < nvars; j++){
            rrmse_month[current_month][j] = rrmse[j];
	  }
	}
      }
      /** Compute the average rrmse for this sampling frequency **/
      for (current_month = 0; current_month < 12; current_month++){
        fprintf(fp_metrics,"%f %d ",dt_rs,current_month+1);
        for (j = 0; j < nvars; j++){
          rrmse_month[current_month][j] = rrmse_month[current_month][j]/(dt_rs/dt_bs);
          fprintf(fp_metrics,"%f ",rrmse_month[current_month][j]);
        }
        fprintf(fp_metrics,"\n");
      }
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
