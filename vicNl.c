#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
#include <vic_ensemble_io.h>
//#include <global.h>

static char vcid[] = "$Id: vicNl.c,v 5.14.2.20 2012/02/05 00:15:44 vicadmin Exp $";

int vicNl_cell(double *output_array,soil_con_struct soil_con, 
               veg_con_struct *veg_con,dmy_struct *dmy,atmos_data_struct *atmos,struct netcdf_output_struct netcdf_output)//,global_param_struct global_param)
/**********************************************************************
	vicNl.c		Dag Lohmann		January 1996

  This program controls file I/O and variable initialization as well as
  being the primary driver for the model.

  For details about variables, input files and subroutines check:
	http://ce.washington.edu/~hydro/Lettenmaier/Models/VIC/VIC_home.html

  UNITS: unless otherwise marked:
         all water balance components are in mm
	 all energy balance components are in mks
	 depths, and lengths are in m

**********************************************************************/
{

  extern veg_lib_struct *veg_lib;
  extern option_struct options;
  extern Error_struct Error;
  extern global_param_struct global_param;

  /** Variable Declarations **/

  char                     NEWCELL;
  char                     LASTREC;
  char                     MODEL_DONE;
  char                     RUN_MODEL;
  char                    *init_STILL_STORM;
  char                     ErrStr[MAXSTRING];
  int                      rec, i, j;
  int                      veg;
  int                      dist;
  int                      band;
  int                      Ndist;
  int                      Nveg_type;
  int                      cellnum;
  int                      index;
  int                     *init_DRY_TIME;
  int                      Ncells;
  int                      cell_cnt;
  int                      startrec;
  int                      ErrorFlag;
  float                    mu;
  double                   storage;
  double                   veg_fract;
  double                   band_fract;
  double                   Clake;
  //dmy_struct              *dmy;
  //atmos_data_struct       *atmos;
  //veg_con_struct          *veg_con;
  //soil_con_struct          soil_con;
  dist_prcp_struct         prcp; /* stores information about distributed 
				    precipitation */
  //filenames_struct         filenames;
  filep_struct             filep;
  lake_con_struct          lake_con;
  out_data_file_struct     *out_data_files;
  out_data_struct          *out_data;
  save_data_struct         save_data;
  //function call implementation
  int 			   out_pos;
  int 			   m;
  
//#if VERBOSE
//  display_current_settings(DISP_VERSION,(filenames_struct*)NULL,(global_param_struct*)NULL);
//#endif

  /** Set up output data structures **/
  out_data = create_output_list();
  out_data_files = set_output_defaults(out_data);
  //fclose(filep.globalparam);
  //filep.globalparam = open_file(filenames.global,"r");
  //parse_output_info(&filenames, filep.globalparam, &out_data_files, out_data);
  //return;

  /** Check and Open Files **/
  //check_files(&filep, &filenames);

#if !OUTPUT_FORCE

  /** Read Vegetation Library File **/
  //veg_lib = read_veglib(filep.veglib,&Nveg_type);

#endif // !OUTPUT_FORCE

  /** Initialize Parameters **/
  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;
  cellnum = -1;

  /** Make Date Data Structure **/
  //dmy      = make_dmy(&global_param);

  /** allocate memory for the atmos_data_struct **/
  //alloc_atmos(global_param.nrecs, &atmos);

  /** Initial state **/
  startrec = 0;
//#if !OUTPUT_FORCE
//  if ( options.INIT_STATE ) 
//    filep.init_state = check_state_file(filenames.init_state, dmy, 
//					 &global_param, options.Nlayer, 
//					 options.Nnode, &startrec);

  /** open state file if model state is to be saved **/
//  if ( options.SAVE_STATE && strcmp( filenames.statefile, "NONE" ) != 0 )
//    filep.statefile = open_state_file(&global_param, filenames, options.Nlayer,
//                                         options.Nnode);
 // else 
 filep.statefile = NULL;

//#endif // !OUTPUT_FORCE

  /************************************
    Run Model for all Active Grid Cells
    ************************************/
  MODEL_DONE = FALSE;
  RUN_MODEL = TRUE;
  cell_cnt=0;
 // while(!MODEL_DONE) {
    //printf("%s\n%s\n%d\n%d\n%d\n",filep.soilparam,filenames.soil_dir,cell_cnt,RUN_MODEL,MODEL_DONE);
    //soil_con = read_soilparam(filep.soilparam, filenames.soil_dir, &cell_cnt, &RUN_MODEL, &MODEL_DONE);

    if(RUN_MODEL) {

#if QUICK_FS
      /** Allocate Unfrozen Water Content Table **/
      if(options.FROZEN_SOIL) {
	for(i=0;i<MAX_LAYERS;i++) {
	  soil_con.ufwc_table_layer[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));
	  for(j=0;j<QUICK_FS_TEMPS+1;j++) 
	    soil_con.ufwc_table_layer[i][j] = (double *)malloc(2*sizeof(double));
	}
	for(i=0;i<MAX_NODES;i++) {
	  soil_con.ufwc_table_node[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));

	  for(j=0;j<QUICK_FS_TEMPS+1;j++) 
	    soil_con.ufwc_table_node[i][j] = (double *)malloc(2*sizeof(double));
	}
      }
#endif /* QUICK_FS */

      NEWCELL=TRUE;
      cellnum++;

#if !OUTPUT_FORCE

      /** Read Grid Cell Vegetation Parameters **/
      //veg_con = read_vegparam(filep.vegparam, soil_con.gridcel,
      //                       Nveg_type);
      //calc_root_fractions(veg_con, &soil_con);

      //if ( options.LAKES ) 
//	lake_con = read_lakeparam(filep.lakeparam, soil_con, veg_con);

#endif // !OUTPUT_FORCE

      /** Build Gridded Filenames, and Open **/
      //make_in_and_outfiles(&filep, &filenames, &soil_con);//, out_data_files);

      //if (options.PRT_HEADER) {
        /** Write output file headers **/
      //  write_header(out_data_files, out_data, dmy, global_param);
      //}

#if !OUTPUT_FORCE

      /** Read Elevation Band Data if Used **/
      //read_snowband(filep.snowband, &soil_con);

      /** Make Precipitation Distribution Control Structure **/
      prcp     = make_dist_prcp(veg_con[0].vegetat_type_num);

#endif // !OUTPUT_FORCE

      /**************************************************
         Initialize Meteological Forcing Values That
         Have not Been Specifically Set
       **************************************************/

#if VERBOSE
      fprintf(stderr,"Initializing Forcing Data\n");
#endif /* VERBOSE */

//      initialize_atmos(atmos, dmy, filep.forcing,
//#if OUTPUT_FORCE
//		       &soil_con, out_data_files, out_data); 
//#else /* OUTPUT_FORCE */
//                       &soil_con); 
//#endif /* OUTPUT_FORCE */

#if !OUTPUT_FORCE

      /**************************************************
        Initialize Energy Balance and Snow Variables 
      **************************************************/

#if VERBOSE
      fprintf(stderr,"Model State Initialization\n");
#endif /* VERBOSE */
      ErrorFlag = initialize_model_state(&prcp, dmy[0], &global_param, filep, 
			     soil_con.gridcel, veg_con[0].vegetat_type_num,
			     options.Nnode, Ndist, 
			     atmos[0].air_temp[NR],
			     &soil_con, veg_con, lake_con,
			     &init_STILL_STORM, &init_DRY_TIME);
      if ( ErrorFlag == ERROR ) {
	if ( options.CONTINUEONERROR == TRUE ) {
	  // Handle grid cell solution error
	  fprintf(stderr, "ERROR: Grid cell %i failed in record %i so the simulation has not finished.  An incomplete output file has been generated, check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
	  //break;
	} else {
	  // Else exit program on cell solution error as in previous versions
	  sprintf(ErrStr, "ERROR: Grid cell %i failed in record %i so the simulation has ended. Check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
	  vicerror(ErrStr);
	}
      }
      
#if VERBOSE
      fprintf(stderr,"Running Model\n");
#endif /* VERBOSE */

      /** Update Error Handling Structure **/
      Error.filep = filep;
      Error.out_data_files = out_data_files;

      /** Initialize the storage terms in the water and energy balances **/
      /** Sending a negative record number (-global_param.nrecs) to dist_prec() will accomplish this **/
      ErrorFlag = dist_prec(&atmos[0], &prcp, &soil_con, veg_con,
		  &lake_con, dmy, &global_param, &filep, out_data_files,
		  out_data, &save_data, -global_param.nrecs, cellnum,
                  NEWCELL, LASTREC, init_STILL_STORM, init_DRY_TIME);

      /** Pass initialized values to output array**/
      rec = 0;
      output_array[0] = out_data[OUT_EVAP].data[0];//Evaporation [mm]
      output_array[1] = out_data[OUT_RUNOFF].data[0];//Surface Runoff [mm]
      output_array[2] = out_data[OUT_BASEFLOW].data[0];//Baseflow [mm]
      output_array[3] = out_data[OUT_PREC].data[0];//Precipitation [mm]
      output_array[4] = out_data[OUT_SOIL_MOIST].data[0];//Soil moisture 1 [mm]
      output_array[5] = out_data[OUT_SOIL_MOIST].data[1];//Soil moisture 2 [mm]
      output_array[6] = out_data[OUT_SOIL_MOIST].data[2];//Soil moisture 3 [mm]

      /******************************************
	Run Model in Grid Cell for all Time Steps
	******************************************/

      for ( rec = startrec ; rec < global_param.nrecs; rec++ ) {

        if ( rec == global_param.nrecs - 1 ) LASTREC = TRUE;
        else LASTREC = FALSE;
	//printf("%f\n",out_data[103].data[0]);
        ErrorFlag = dist_prec(&atmos[rec], &prcp, &soil_con, veg_con,
		  &lake_con, dmy, &global_param, &filep,
		  out_data_files, out_data, &save_data, rec, cellnum,
                  NEWCELL, LASTREC, init_STILL_STORM, init_DRY_TIME);

        /** Place the data that we desire into the output array**/
        netcdf_output.SNOW_DEPTH[rec] = out_data[OUT_SNOW_DEPTH].data[0];
        netcdf_output.SWE[rec] = out_data[OUT_SWE].data[0];
        netcdf_output.EVAP_BARE[rec] = out_data[OUT_EVAP_BARE].data[0];
        netcdf_output.EVAP_CANOP[rec] = out_data[OUT_EVAP_CANOP].data[0];
        netcdf_output.SNOW_MELT[rec] = out_data[OUT_SNOW_MELT].data[0];
        netcdf_output.TRANSP_VEG[rec] = out_data[OUT_TRANSP_VEG].data[0];
        netcdf_output.SURF_TEMP[rec] = out_data[OUT_SURF_TEMP].data[0];
        netcdf_output.GRND_FLUX[rec] = out_data[OUT_GRND_FLUX].data[0];
        netcdf_output.LATENT[rec] = out_data[OUT_LATENT].data[0];
        netcdf_output.NET_LONG[rec] = out_data[OUT_NET_LONG].data[0];
        netcdf_output.NET_SHORT[rec] = out_data[OUT_NET_SHORT].data[0];
        netcdf_output.R_NET[rec] = out_data[OUT_R_NET].data[0];
        netcdf_output.SENSIBLE[rec] = out_data[OUT_SENSIBLE].data[0];
        netcdf_output.AERO_COND[rec] = out_data[OUT_AERO_COND].data[0];
        netcdf_output.SURF_COND[rec] = out_data[OUT_SURF_COND].data[0];
        netcdf_output.EVAP[rec] = out_data[OUT_EVAP].data[0];
        netcdf_output.RUNOFF[rec] = out_data[OUT_RUNOFF].data[0];
        netcdf_output.BASEFLOW[rec] = out_data[OUT_BASEFLOW].data[0];
        netcdf_output.SOIL_MOIST1[rec] = out_data[OUT_SOIL_MOIST].data[0];
        netcdf_output.SOIL_MOIST2[rec] = out_data[OUT_SOIL_MOIST].data[1];
        netcdf_output.SOIL_MOIST3[rec] = out_data[OUT_SOIL_MOIST].data[2];

        for (m=0;m<7;m++){
          out_pos = (rec)*7 + m;
          if(m==0)output_array[out_pos] = out_data[OUT_EVAP].data[0];//Evaporation [mm]
          if(m==1)output_array[out_pos] = out_data[OUT_RUNOFF].data[0];//Surface Runoff [mm]      
          if(m==2)output_array[out_pos] = out_data[OUT_BASEFLOW].data[0];//Baseflow [mm]
          if(m==3)output_array[out_pos] = out_data[OUT_PREC].data[0];//Precipitation [mm]
          if(m==4)output_array[out_pos] = out_data[OUT_SOIL_MOIST].data[0];//Soil moisture 1 [mm]
          if(m==5)output_array[out_pos] = out_data[OUT_SOIL_MOIST].data[1];//Soil moisture 2 [mm] 
          if(m==6)output_array[out_pos] = out_data[OUT_SOIL_MOIST].data[2];//Soil moisture 3 [mm]
	}

        if ( ErrorFlag == ERROR ) {
          if ( options.CONTINUEONERROR == TRUE ) {
            // Handle grid cell solution error
            fprintf(stderr, "ERROR: Grid cell %i failed in record %i so the simulation has not finished.  An incomplete output file has been generated, check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
            break;
          } else {
	    // Else exit program on cell solution error as in previous versions
            sprintf(ErrStr, "ERROR: Grid cell %i failed in record %i so the simulation has ended. Check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
            vicerror(ErrStr);
	  }
        }

        NEWCELL=FALSE;
	for ( veg = 0; veg <= veg_con[0].vegetat_type_num; veg++ )
	  init_DRY_TIME[veg] = -999;

      }	/* End Rec Loop */

#endif /* !OUTPUT_FORCE */

      //close_files(&filep,out_data_files,&filenames); 

#if !OUTPUT_FORCE
#if QUICK_FS
      if(options.FROZEN_SOIL) {
	for(i=0;i<MAX_LAYERS;i++) {
	  for(j=0;j<6;j++) 
	    free((char *)soil_con.ufwc_table_layer[i][j]);
	  free((char *)soil_con.ufwc_table_layer[i]);
	}
	for(i=0;i<MAX_NODES;i++) {
	  for(j=0;j<6;j++) 
	    free((char *)soil_con.ufwc_table_node[i][j]);
	  free((char *)soil_con.ufwc_table_node[i]);
	}
      }
#endif
 /* QUICK_FS */
      free_dist_prcp(&prcp,veg_con[0].vegetat_type_num);
      //free_vegcon(&veg_con);
      //free((char *)soil_con.AreaFract);
      //free((char *)soil_con.BandElev);
      //free((char *)soil_con.Tfactor);
      //free((char *)soil_con.Pfactor);
      //free((char *)soil_con.AboveTreeLine);
      free((char*)init_STILL_STORM);
      free((char*)init_DRY_TIME);
#endif /* !OUTPUT_FORCE */
    }	/* End Run Model Condition */
 // } 	/* End Grid Loop */

  /** cleanup **/
  free_out_data_files(&out_data_files);
  free_out_data(&out_data);

  return EXIT_SUCCESS;
}	/* End Main Program */
