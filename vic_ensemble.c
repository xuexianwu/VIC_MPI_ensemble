#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
#include <global.h>
#include <vic_ensemble_io.h>
#include <mpi.h>
#include "borg/borg.h"
#include <time.h>

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
float obs[12];
                        
                        
// Calibration wrapper function to be used with Borg
// Receives parameters and writes objective values to an array
void vic_calibration_wrapper(double* vars, double* objs);

// Make these variables global so the wrapper function can use them
// (Will this wreck anything else?)
soil_con_struct soil_con; 
veg_con_struct *veg_con;
atmos_data_struct *atmos; 
dmy_struct *dmy; 
double *observed_data;
double *simulated_data;
forcing_cell_struct *forcing_cell;

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
  
  char MODEL_DONE;
  char RUN_MODEL;
  int cell_cnt;
  int Nveg_type;
  //Resampling variables
  double dt_rs;
  double dt_bs = 1; //hours
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
  char metrics_root[MAXSTRING];
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
  fscanf(global_fp,"%s %s %d",tmp_s,&forcing_name.tair,&tmp_i); grads_file.nt_tair = tmp_i;
  //printf("%s %s\n",tmp_s,forcing_name.tair);
  fscanf(global_fp,"%s %s %d",tmp_s,&forcing_name.prec,&tmp_i); grads_file.nt_prec = tmp_i;
  //printf("%s %s\n",tmp_s,forcing_name.prec);
  fscanf(global_fp,"%s %s %d",tmp_s,&forcing_name.wind,&tmp_i); grads_file.nt_wind = tmp_i;
  //printf("%s %s\n",tmp_s,forcing_name.wind);
  fscanf(global_fp,"%s %s %d",tmp_s,&forcing_name.shum,&tmp_i); grads_file.nt_shum = tmp_i;
  //printf("%s %s\n",tmp_s,forcing_name.shum);
  fscanf(global_fp,"%s %s %d",tmp_s,&forcing_name.pres,&tmp_i); grads_file.nt_pres = tmp_i;
  //printf("%s %s\n",tmp_s,forcing_name.pres);
  fscanf(global_fp,"%s %s %d",tmp_s,&forcing_name.lwdown,&tmp_i); grads_file.nt_lwdown = tmp_i;
  //printf("%s %s\n",tmp_s,forcing_name.lwdown);
  fscanf(global_fp,"%s %s %d",tmp_s,&forcing_name.swdown,&tmp_i); grads_file.nt_swdown = tmp_i;
  //printf("%d\n",grads_file.nt_swdown);
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
  //sprintf (metrics_filename,"%s/metrics_output_%d.txt",metrics_root,myid);
  //Number of cells
  fscanf(global_fp,"%s %d",tmp_s,&soil_ncells);
  //Number of resampling intervals
  fscanf(global_fp,"%s %d",tmp_s,&nrs);

  //Forcing variables
  forcing_filep_struct forcing_filep;
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

  /** Open up the observations files **/
  FILE *obs_fp;
  obs_fp = fopen("/scratch/sciteam/nchaney/data/PSU_AERO_PU/1.0deg/OBS/GRDC_Monthly_Climatology.txt","r");

  /** Initialize the structures to hold the atmospheric data (before passing to VIC)**/
  int ncells = 1; //The number of cells to read in at once
  allocate_forcing(&forcing_cell,&grads_file,ncells);

  //Assign vegetation/soils names
  global_param.nrecs = grads_file.nt; //HARD CODED
  observed_data = (double *) malloc(sizeof(double)*global_param.nrecs*nvars);
  simulated_data = (double *) malloc(sizeof(double)*global_param.nrecs*nvars);

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
  //FILE *fp_metrics;
  //FILE *fp_output;
  double *rrmse = (double *) malloc(sizeof(double)*nvars);
  double rrmse_month[12][nvars];

  //fp_metrics = fopen(metrics_filename,"w");
  /** Iterate for all cells in the soil file **/
  icell = myid;
  int linen = -1;
  int current_month;
  
  while (icell < soil_ncells){
  
    for (i = 0; i < nvars; i++){
      rrmse[i] = 0.0;
      for (current_month = 0; current_month<12; current_month++){
        rrmse_month[current_month][i] = 0.0;
      }
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
    
	// Copy the forcing data to the VIC structure
    initialize_atmos_BLUEWATERS(atmos, dmy, forcing_cell,&soil_con);
    
    // Extract monthly runoff observations from file
    rewind(obs_fp);
    fgetc(obs_fp);
    fgetc(obs_fp);//Read the first line
    float lat,lon,tmp,linem; 
    int flag_cell = 1;
    linem = 0;
    while (linem <  14548){
      fscanf(obs_fp,"%f,%f",&lat,&lon);
      fgetc(obs_fp);
      if (lat == soil_con.lat && lon == soil_con.lng){
        //printf("%f %f\n",lat,lon);
        for (int i = 0; i < 12; i++){
         fscanf(obs_fp,"%f",&tmp);
         fgetc(obs_fp);
         obs[i] = tmp;
         //printf("%f \n",obs[i]);
         }
        //printf("\n");
        flag_cell = 0;
        break;
      }
      linem = linem + 1;
    } 
    if (flag_cell == 1){
      //We don't have runoff observations for this cell so there is no need to calibrate...
      printf("There are no observations for this cell\n");
      icell = icell + np;
      continue; 
    }
	
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Changes for hypercube version start here
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // For each cell ... (this is inside the cell loop)
    // (1) Read in parameter values from hypercube file
    // (2) Evaluate the model for this parameter set
	
	// Pass in size of hypercube sample on the command line
	int num_hypercube = atoi(argv[2]);
	
	// Open the file of hypercube samples for reading into an array
	FILE *hcube_fp;
    hcube_fp = fopen("/u/sciteam/jdh33/projects/VIC/vic_hypercube_100.txt", "r");
	double *hcube_params = (double *) malloc(sizeof(double)*4); // 4 parameters
	
	// Open the output file for writing objective(s) after each evaluation
	FILE *hcube_output_fp;
	char hcube_output_filename[MAXSTRING];
	sprintf(hcube_output_filename, "%s/hcube_lat_%f_long_%f_cell_%d.txt", metrics_root, lat, lon, icell);
	hcube_output_fp = fopen(hcube_output_filename, "w");
	
	double *hcube_obj = (double *) malloc(sizeof(double)*1); // 1 objective
	
	for(i = 0; i < num_hypercube; i++) {
	
		fscanf(hcube_fp,"%f %f %f %f", &hcube_params[0], &hcube_params[1], &hcube_params[2], &hcube_params[3]);
		fgetc(hcube_fp); // skip EOL character
		
		// test values
		printf("%f %f %f %f\n", hcube_params[0], hcube_params[1], hcube_params[2], hcube_params[3]);
		
		// Run the model with these parameters. record objective(s).
		vic_calibration_wrapper(hcube_params, hcube_obj);
		
		fprintf(hcube_output_fp, "%f\n", hcube_obj[0]);
		
	}
	
    fclose(hcube_fp);
	fclose(hcube_output_fp);
	free(hcube_params);
    free(hcube_obj);
    
    //Next cell
    icell = icell + np;
  
  }
  
  
  //fclose(fp_metrics);
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


// Calibration wrapper function to be used with Borg
// Receives parameters and writes objective values to an array
void vic_calibration_wrapper(double* vars, double* objs) {
    
    float sim[12],count[12],qbase,qsurf;
    time_t t;
    int n = global_param.nrecs;
    int dt = 1;
    int nvars = 7;
    struct tm gtime;
    struct tm gtime_original;
    gtime.tm_sec = 0;
    gtime.tm_min = 0;
    gtime.tm_hour = 0;
    gtime.tm_mday = 0;
    gtime.tm_mon = 0;
    gtime.tm_year = 2000 - 1900;
    gtime_original = gtime;
    // Initialize all the simulated data
    for (int i = 0; i < n*nvars; i++){
      simulated_data[i] = 0.0;
    }
    // Initialize all the monthly count and simulated data (REASONS WHY I HATE C)
    for (int i = 0; i< 12; i++){
      count[i] = 0.0;
      sim[i] = 0.0;
    }
    //if(m==1)output_array[out_pos] = out_data[OUT_RUNOFF].data[0];//Surface Runoff [mm]      
    //if(m==2)output_array[out_pos] = out_data[OUT_BASEFLOW].data[0];//Baseflow [mm]

    // Set parameter values in soil struct
    soil_con.b_infilt = vars[0];
    soil_con.Ds = vars[1];
    soil_con.Dsmax = vars[2];
    soil_con.Ws = vars[3];
    
    // Run the model
    vicNl_cell(simulated_data, soil_con, veg_con, dmy, atmos);

    // Find the monthly averages of simulated runoff (baseflow + surface runoff)

    int tcount = 0;
    int mon = 0;
    t = mktime(&gtime_original);
    for (int i = 0; i < n; i++){
      t = t + 3600*dt;
      gmtime_r(&t,&gtime);
      qsurf = simulated_data[i*nvars+1];
      qbase = simulated_data[i*nvars+2];
      sim[gtime.tm_mon] = sim[gtime.tm_mon] + qsurf + qbase;
      if (mon != gtime.tm_mon){
        count[gtime.tm_mon] = count[gtime.tm_mon] + 1;
        mon = gtime.tm_mon;
      }
    }
    
    // Calculate objective values by comparing simulated_data to observed_data
    // You just need to compare the arrays obs and sim (12 values per array)
    // Fill out the array objs[] with these values
    
    // For now, make up an objective function
    objs[0] = vars[0]*vars[1] + vars[2]*vars[3];
    
}
