#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
#include <global.h>
#include <vic_ensemble_io.h>
#include <mpi.h>
//#include "borg/borg.h"
#include <time.h>
#include <netcdf.h>
#include <math.h>

int vicNl_cell();
void resample();
void Comparision_Statistics();
void open_forcing_files();
int open_netcdf_forcing_files();
void close_forcing_files();
void extract_cell();
void extract_cell_netcdf();
void allocate_forcing();
void deallocate_forcing();
void initialize_atmos_BLUEWATERS(atmos_data_struct *, dmy_struct *,// FILE **,
                        forcing_cell_struct *forcing_cell,
                        soil_con_struct *);
float obs[12];
                        
                        
// Calibration wrapper function to be used with Borg
// Receives parameters and writes objective values to an array
void vic_calibration_wrapper(double* vars, double* objs,grads_file_struct *);

// Make these variables global so the wrapper function can use them
// (Will this wreck anything else?)
soil_con_struct soil_con; 
veg_con_struct *veg_con;
atmos_data_struct *atmos; 
dmy_struct *dmy; 
double *observed_data;
double *simulated_data;
forcing_cell_struct *forcing_cell;
float *rmin;
int Nveg_type;

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
  int i,j,p;
  filenames_struct names;
  filep_struct filep;
  
  char MODEL_DONE;
  char RUN_MODEL;
  int cell_cnt;
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
  char netcdf_forcing_filename[MAXSTRING];
  //grads_file.year = 2000;
  //grads_file.month = 1;
  //grads_file.day = 1;
  grads_file.hour = 0;
  //grads_file.dt = 3;

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
  fscanf(global_fp,"%s %d",tmp_s,&tmp_i); grads_file.dt = tmp_i;
  fscanf(global_fp,"%s %d",tmp_s,&tmp_i); grads_file.year = tmp_i;
  fscanf(global_fp,"%s %d",tmp_s,&tmp_i); grads_file.month = tmp_i;
  fscanf(global_fp,"%s %d",tmp_s,&tmp_i); grads_file.day = tmp_i;
  //printf("%s %d\n",tmp_s,grads_file.nt);
  //Land Info
  fscanf(global_fp,"%s %s",tmp_s,&names.soil);
  //printf("%s %s\n",tmp_s,names.soil);
  fscanf(global_fp,"%s %s",tmp_s,&names.veglib);
  //printf("%s %s\n",tmp_s,names.veglib);
  fscanf(global_fp,"%s %s",tmp_s,&names.veg);
  //printf("%s %s\n",tmp_s,names.veg);
  //snowbands file
  fscanf(global_fp,"%s %s",tmp_s,&names.snowband);
  //netcdf forcing file
  fscanf(global_fp,"%s %s %d",tmp_s,&netcdf_forcing_filename,&tmp_i); grads_file.nt_netcdf = tmp_i;
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
  global_param = get_global_param(&names,&grads_file); //Fill up the global parameter file

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
  
  /** Copy the rmin **/
  rmin = (float *) malloc(sizeof(float)*Nveg_type);
  for (i=0;i<Nveg_type;i++){
   rmin[i] = veg_lib[i].rmin;
  }

  /** Set up all variables for the iteration through all cells **/
  int icell;
  //soil_ncells = 2//15836; TEMPORARY FIX
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
  
  while (icell < soil_ncells) {
  
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
    /** Read Elevation Band Data if Used **/
    read_snowband(filep.snowband, &soil_con);
    /** Read the vegetation parameters **/
    veg_con = read_vegparam(filep.vegparam,soil_con.gridcel,Nveg_type);
    calc_root_fractions(veg_con, &soil_con);
    /** Read in the forcing data for a single cell **/
    //extract_cell(&forcing_filep,&grads_file,soil_con.lat,soil_con.lng,forcing_cell);
    /** Get data for the cell **/
    int forcing_ncid;
    forcing_ncid = open_netcdf_forcing_files(&netcdf_forcing_filename);
    extract_cell_netcdf(forcing_ncid,&grads_file,forcing_cell,soil_con.gridcel);
    nc_close(forcing_ncid);

    //extract_cell_netcdf(forcing_file,&grads_file,soil_celln);
    
	  // Copy the forcing data to the VIC structure
    initialize_atmos_BLUEWATERS(atmos, dmy, forcing_cell,&soil_con);
    
    // Extract monthly runoff observations from file
    //rewind(obs_fp);
    //fgetc(obs_fp);
    //fgetc(obs_fp);//Read the first line
    float lat,lon,tmp,linem; 
    int flag_cell = 1;
    linem = 0;
	
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Changes for hypercube version start here
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // For each cell ... (this is inside the cell loop)
    // (1) Read in parameter values from hypercube file
    // (2) Evaluate the model for this parameter set
	
  	// Pass in size of hypercube sample on the command line
    // and the number of the hypercube file (0-9)
  	int num_hypercube = atoi(argv[2]);
  	int filenum = atoi(argv[3]);
  	// Open the file of hypercube samples for reading into an array
  	FILE *hcube_fp;
    char hcube_input_filename[MAXSTRING];
    sprintf(hcube_input_filename, "/u/sciteam/jdh33/projects/VIC/vic_hypercube_9_250_%d.txt", filenum);
    hcube_fp = fopen(hcube_input_filename, "r");
  	double *hcube_params = (double *) malloc(sizeof(double)*6); // 6 parameters
  	
  	// Open the output file for writing objective(s) after each evaluation
  	// The filename will contain 6 digits after the decimals by default
  	FILE *hcube_output_fp;
  	char hcube_output_filename[MAXSTRING];
  	char hcube_output_nc4_filename[MAXSTRING];
  	char hcube_output_nc3_filename[MAXSTRING];
  	sprintf(hcube_output_filename, "%s/file_%d/txt/hcube_lat_%f_long_%f.txt", metrics_root, filenum, soil_con.lat, soil_con.lng);
  	sprintf(hcube_output_nc4_filename, "%s/file_%d/nc/hcube_lat_%f_long_%f.nc4", metrics_root, filenum, soil_con.lat, soil_con.lng);
  	sprintf(hcube_output_nc3_filename, "%s/file_%d/nc/hcube_lat_%f_long_%f.nc3", metrics_root, filenum, soil_con.lat, soil_con.lng);
  	hcube_output_fp = fopen(hcube_output_filename, "w");
  	
  	double *hcube_obj = (double *) malloc(sizeof(double)*12); // 1 objective (12 to print monthly sim values)
  	
  	// test: print the observations to the output file first row
  	/*for(i = 0; i < 12; i++) {
  		fprintf(hcube_output_fp, "%f", obs[i]);
  		if(i < 11) fprintf(hcube_output_fp, " ");
  		else fprintf(hcube_output_fp, "\n");
  	}*/
        // Initialize netcdf file
        int status;
        int ncid;
        status = nc_create(hcube_output_nc4_filename,NC_CLOBBER|NC_NETCDF4, &ncid);
        //Define the netcdf dimensions
        int timeid;
        int ens_id;
        status = nc_def_dim(ncid,"time",global_param.nrecs,&timeid);
        status = nc_def_dim(ncid,"ensemble",num_hypercube,&ens_id);
        //Define the variables
        int prec_id,qsurf_id,qbase_id,evap_id,sm1_id,sm2_id,sm3_id;
        int var_dimids[2];
        int prec_dimids[1];
        var_dimids[0] = ens_id;
        var_dimids[1] = timeid;
 	prec_dimids[0] = timeid;
        status = nc_def_var(ncid,"prec",NC_FLOAT,2,var_dimids,&prec_id);
        status = nc_def_var(ncid,"qbase",NC_FLOAT,2,var_dimids,&qbase_id);
        status = nc_def_var(ncid,"qsurf",NC_FLOAT,2,var_dimids,&qsurf_id);
        status = nc_def_var(ncid,"evap",NC_FLOAT,2,var_dimids,&evap_id);
        status = nc_def_var(ncid,"sm1",NC_FLOAT,2,var_dimids,&sm1_id);
        status = nc_def_var(ncid,"sm2",NC_FLOAT,2,var_dimids,&sm2_id);
        status = nc_def_var(ncid,"sm3",NC_FLOAT,2,var_dimids,&sm3_id);
        //Set up compression
        nc_def_var_deflate(ncid,evap_id,1,1,3);
        nc_def_var_deflate(ncid,qbase_id,1,1,3);
        nc_def_var_deflate(ncid,qsurf_id,1,1,3);
        nc_def_var_deflate(ncid,prec_id,1,1,3);
        nc_def_var_deflate(ncid,sm1_id,1,1,3);
        nc_def_var_deflate(ncid,sm2_id,1,1,3);
        nc_def_var_deflate(ncid,sm3_id,1,1,3);
        //Close the file
        status = nc_close(ncid);
	//Define output variables
        int itime;
        float qsurf[global_param.nrecs];
        float qbase[global_param.nrecs];
        float evap[global_param.nrecs];
        float prec[global_param.nrecs];
        float sm1[global_param.nrecs];
        float sm2[global_param.nrecs];
        float sm3[global_param.nrecs];
        size_t index[2],count[2],start[2];
	
  	for(i = 0; i < num_hypercube; i++) {
  	  
      // Read in 9 parameters from sample file
      for(p = 0; p < 9; p++) {
  		  fscanf(hcube_fp,"%lf", &hcube_params[p]);
  		}
      fgetc(hcube_fp); // skip EOL character
  		
  		// print run details to stdout to keep track of what's happening
  		// printf("Cell %d, Sim %d: %f %f %f %f\n", icell, i, hcube_params[0], hcube_params[1], hcube_params[2], hcube_params[3]);
  		
  		// Run the model with these parameters. record objective(s).
  		vic_calibration_wrapper(hcube_params, hcube_obj, &grads_file);
  		
  		for(j = 0; j < 12; j++) {
  		  fprintf(hcube_output_fp, "%f", hcube_obj[j]);
  		  if(j < 11) fprintf(hcube_output_fp, " ");
  		  else fprintf(hcube_output_fp, "\n");
  		}

		//Place output in netcdf file
                for (itime = 0; itime < global_param.nrecs; itime++){
                 evap[itime] = simulated_data[itime*7+0];
                 qsurf[itime] = simulated_data[itime*7+1];
                 qbase[itime] = simulated_data[itime*7+2];
                 prec[itime] = simulated_data[itime*7+3];
                 sm1[itime] = simulated_data[itime*7+4];
                 sm2[itime] = simulated_data[itime*7+5];
                 sm3[itime] = simulated_data[itime*7+6];
                 //index[0] = i;
		 //index[1] = itime;
                 //status = nc_put_var1_float(ncid,prec_id,index,&prec[itime]);
                 //status = nc_put_var1_float(ncid,evap_id,index,&evap[itime]);
                 //status = nc_put_var1_float(ncid,sm1_id,index,&sm1[itime]);
                 //status = nc_put_var1_float(ncid,sm2_id,index,&sm2[itime]);
                 //status = nc_put_var1_float(ncid,sm3_id,index,&sm3[itime]);
                 //status = nc_put_var1_float(ncid,qbase_id,index,&qbase[itime]);
                 //status = nc_put_var1_float(ncid,qsurf_id,index,&qsurf[itime]);
                }
 		count[0] = 1;
                count[1] = global_param.nrecs;
                start[0] = i;
                start[1] = 0;

                //Open file
                nc_open(hcube_output_nc4_filename,NC_WRITE,&ncid);
                nc_inq_varid(ncid,"prec",&prec_id);
		nc_put_vara_float(ncid,prec_id,start,count,&prec[0]);
                nc_inq_varid(ncid,"evap",&evap_id);
		status = nc_put_vara_float(ncid,evap_id,start,count,&evap[0]);
                nc_inq_varid(ncid,"sm1",&sm1_id);
		status = nc_put_vara_float(ncid,sm1_id,start,count,&sm1[0]);
                nc_inq_varid(ncid,"sm2",&sm2_id);
		status = nc_put_vara_float(ncid,sm2_id,start,count,&sm2[0]);
                nc_inq_varid(ncid,"sm3",&sm3_id);
		status = nc_put_vara_float(ncid,sm3_id,start,count,&sm3[0]);
                nc_inq_varid(ncid,"qbase",&qbase_id);
		status = nc_put_vara_float(ncid,qbase_id,start,count,&qbase[0]);
                nc_inq_varid(ncid,"qsurf",&qsurf_id);
		status = nc_put_vara_float(ncid,qsurf_id,start,count,&qsurf[0]);
       
                // Close the netcdf file
                status = nc_close(ncid);
  		
  		// buffer flush after each evaluation
  		fflush(hcube_output_fp);
  		
  	}

      
       // Convert from nc4 to nc3 (This is embarassing... but it works!)
       char system_call_string[MAXSTRING];
       sprintf(system_call_string,"/opt/cray/netcdf/4.1.3/cray/73/bin/nccopy -k 1 %s %s", hcube_output_nc4_filename,hcube_output_nc3_filename);
       system(system_call_string);
       // Remove netcdf 4 file
       sprintf(system_call_string,"rm %s",hcube_output_nc4_filename);
       system(system_call_string);
       sprintf(system_call_string, "chmod 770 %s", hcube_output_nc3_filename);
       system(system_call_string);
	
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
  fclose(filep.snowband);
  fclose(filep.soilparam);
  fclose(filep.vegparam);
  fclose(filep.veglib);
  /** Close up the forcing files **/
  //nc_close(forcing_ncid);
  //close_forcing_files(&forcing_filep);
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
void vic_calibration_wrapper(double* vars, double* objs, grads_file_struct *grads_file) {
    
  float sim[12],count[12],qbase,qsurf;
  time_t t;
  int n = global_param.nrecs;
	int i;
  int dt = grads_file->dt;
  int nvars = 7;
  struct tm gtime;
  struct tm gtime_original;
  gtime.tm_sec = 0;
  gtime.tm_min = 0;
  gtime.tm_hour = 0;
  gtime.tm_mday = 0;
  gtime.tm_mon = 0;
  gtime.tm_year = grads_file->year - 1900;//grads_file.year - 1900;
  gtime_original = gtime;
  // Initialize all the simulated data
  for (i = 0; i < n*nvars; i++){
    simulated_data[i] = 0.0;
  }
  // Initialize all the monthly count and simulated data (REASONS WHY I HATE C)
  for (i = 0; i< 12; i++){
    count[i] = 0.0;
    sim[i] = 0.0;
  }
  //if(m==1)output_array[out_pos] = out_data[OUT_RUNOFF].data[0];//Surface Runoff [mm]      
  //if(m==2)output_array[out_pos] = out_data[OUT_BASEFLOW].data[0];//Baseflow [mm]

  // Set parameter values in soil struct
  // Note: to run default parameters, comment these out.
  soil_con.b_infilt = pow(10,vars[0]);
  soil_con.Ds = pow(10,vars[1]);
  soil_con.Dsmax = pow(10,vars[2]);
  soil_con.Ws = vars[3];
  soil_con.depth[1] = vars[4];
  soil_con.depth[2] = vars[5];

  //Here we are adding the new parameters 
  //[vars[6] should be between 0.1 and 10.0]
  for (i = 0; i < Nveg_type; i++){
    veg_lib[i].rmin = pow(10,vars[6])*rmin[i];
  }
  for (i = 0; i < 3; i++){
    //expt - between 1.0 and 30.0
    soil_con.expt[i] = vars[7];
    //ksat - between 100 and 10,000 [mm/day]
    soil_con.Ksat[i] = pow(10,vars[8]);
  }
    
  // Run the model
  vicNl_cell(simulated_data, soil_con, veg_con, dmy, atmos);

  // Find the monthly averages of simulated runoff (baseflow + surface runoff)
  int tcount = 0;
  int mon = 0;
  t = mktime(&gtime_original);
  for (i = 0; i < n; i++){
    t = t + 3600*dt;
    gmtime_r(&t,&gtime);
    if (gtime.tm_year <= (grads_file->year - 1900)){continue;}
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
  
  // For now, save each sim monthly value (average) as an objective
	for(i = 0; i < 12; i++) {
		objs[i] = sim[i]/count[i];
	}
    
}
