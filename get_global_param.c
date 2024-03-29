#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>
#include <vic_ensemble_io.h>
 
static char vcid[] = "$Id: get_global_param.c,v 5.22.2.39 2012/02/06 23:54:00 vicadmin Exp $";

/********************************************************************/
/*			GLOBAL VARIABLES                            */
/********************************************************************/
int NR;		      /* array index for atmos struct that indicates
			 the model step avarage or sum */
int NF;		      /* array index loop counter limit for atmos
			 struct that indicates the SNOW_STEP values */
 
global_param_struct get_global_param(filenames_struct *names,grads_file_struct *grads_file)
/**********************************************************************
  get_global_param	Keith Cherkauer	            March 1998

  This routine reads the VIC model global control file, getting
  values for global parameters, model options, and debugging controls.

  NOTE: any additions or removals of parameters in this file must also
  be made in display_current_settings.c.

  Modifications:
  7-19-96 Modified to read time step		        KAC
  4-5-98  Modified to read model options and debugging
          controls from a single file                   KAC
  01-20-00 modified to work with new radiation estimation routines,
           new simplified frozen soil moisture, and new new open
           format forcing file rad routines.              KAC
  02-27-01 added reads for lake model parameters          KAC
  04-21-03 added parameters for blowing snow algorithm, printing
           lake variables during debugging and reading Bart's 
           new Arno parameters.                           KAC
  11-May-04 Modified to display compile-time and run-time options
	    if VERBOSE is set to TRUE.						TJB
  13-Oct-04 Added validation for GRND_FLUX option.              		TJB
  01-Nov-04 Added validation for Nnodes with QUICK_FLUX option, as
	    part of fix for QUICK_FLUX state file compatibility.		TJB
  2005-03-08 Added EQUAL_AREA option.						TJB
  2005-03-24 Added ALMA_OUTPUT option.						TJB
  2005-04-07 Fixed state file warning check.					TJB
  2005-Apr-13 Added logic for OUTPUT_FORCE option.				TJB
  2005-Apr-23 Changed ARNO_PARAMS to NIJSSEN2001_BASEFLOW			TJB
  2005-11-29 SAVE_STATE is set in global_param (not at compile time)		GCT
  2005-12-06 Moved setting of statename from open_state_file to here.		GCT
  2005-12-07 Added checks for range of STATEMONTH and STATEDAY			GCT
  2005-12-07 Allow user to use NO_FLUX in addition to NOFLUX for NOFLUX in
             global.param.file							GCT
  2006-09-13 Replaced NIJSSEN2001_BASEFLOW with BASEFLOW option.		TJB/GCT
  2006-Sep-23 Implemented flexible output configuration; removed the
              OPTIMIZE and LDAS_OUTPUT options; implemented aggregation of
	      output variables.							TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct;
	      This included moving global->statename to filenames->statefile;
	      also added f_path_pfx to store forcing file path and prefix.	TJB
  2006-Nov-07 Removed LAKE_MODEL option.					TJB
  2007-Jan-03 Added ALMA_INPUT option.						TJB
  2007-Jan-15 Added PRT_HEADER option.						TJB
  2007-Apr-03 Added CONTINUEONERROR option.					GCT
  2007-Apr-24 Added EXP_TRANS option.						JCA
  2007-Apr-24 Added IMPLICIT option.						JCA
  2007-Apr-23 Added initialization of global parameters.			TJB
  2007-Apr-23 Added check for FULL_ENERGY if lake model is run.			TJB
  2007-Aug-08 Added EXCESS_ICE option.						JCA
  2007-Sep-14 Added initialization of names->soil_dir.				TJB
  2007-Oct-10 Added validation of dt, start date, end date, and nrecs.		TJB
  2007-Oct-31 Added validation of input/output files.				TJB
  2008-Jan-25 Removed setting of SNOW_STEP = global.dt for
	      OUTPUT_FORCE == TRUE.						TJB
  2008-Jan-28 Added check that end date falls AFTER start date.			TJB
  2008-Mar-12 Relocated code validating IMPLICIT and EXCESS_ICE options.	TJB
  2008-Apr-21 Added SNOW_ALBEDO option.						KAC via TJB
  2008-Apr-21 Added SNOW_DENSITY option.					TJB
  2009-Jan-12 Added COMPUTE_TREELINE and JULY_TAVG_SUPPLIED options.		TJB
  2009-Jan-16 Modified aero_resist_used and Ra_used to become arrays of
	      two elements (surface and overstory); added
	      options.AERO_RESIST_CANSNOW.					TJB
  2009-May-17 Added AR_406_LS to options.AERO_RESIST_CANSNOW.			TJB
  2009-May-17 Added options.MIN_LIQ.						TJB
  2009-May-18 Added options.PLAPSE.						TJB
  2009-May-20 Added options.GRND_FLUX_TYPE.					TJB
  2009-May-22 Added TFALLBACK value to options.CONTINUEONERROR.			TJB
  2009-Jun-15 Changed order of options to match global parameter file.		TJB
  2009-Aug-29 Now handles commented lines that are indented.			TJB
  2009-Sep-19 Moved TFALLBACK to its own separate option.			TJB
  2009-Nov-15 Added prohibition of QUICK_FLUX=TRUE when IMPLICIT=TRUE or
	      EXP_TRANS=TRUE.							TJB
  2009-Dec-11 Removed min_liq and options.MIN_LIQ.				TJB
  2010-Apr-28 Replaced GLOBAL_LAI with VEGPARAM_LAI and LAI_SRC.		TJB
  2011-May-31 Removed options.GRND_FLUX.					TJB
  2011-Jun-03 Added options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.					TJB
  2011-Jul-05 Now QUICK_FLUX=FALSE is prohibited in water balance mode
	      (when both FULL_ENERGY and FROZEN_SOIL are FALSE).		TJB
  2011-Nov-04 Added options to access new forcing estimation features.		TJB
  2012-Jan-16 Removed LINK_DEBUG code						BN
  2012-Jan-28 Removed AR_COMBO and GF_FULL.					TJB
  2012-Jan-28 IMPLICIT is now set to TRUE by default if FROZEN_SOIL is
	      TRUE.								TJB
  2012-Jan-28 Changed default values of MIN_WIND_SPEED, MIN_RAIN_TEMP,
	      and MAX_SNOW_TEMP to reflect the most commonly-used values.	TJB
**********************************************************************/
{
  extern option_struct    options;
  extern param_set_struct param_set;
  extern int              NF, NR;

  char cmdstr[MAXSTRING];
  char optstr[MAXSTRING];
  char flgstr[MAXSTRING];
  char ErrStr[MAXSTRING];
  int  file_num;
  int  field;
  int  i;
  int  tmpstartdate;
  int  tmpenddate;
  int  lastvalidday;
  int  lastday[] = {
            31, /* JANUARY */
            28, /* FEBRUARY */
            31, /* MARCH */
            30, /* APRIL */
            31, /* MAY */
            30, /* JUNE */
            31, /* JULY */
            31, /* AUGUST */
            30, /* SEPTEMBER */
            31, /* OCTOBER */
            30, /* NOVEMBER */
            31, /* DECEMBER */
        } ;
  global_param_struct global;

  /** Initialize global parameters (that aren't part of the options struct) **/
  global.dt            = MISSING;
  global.nrecs         = MISSING;
  global.startyear     = MISSING;
  global.startmonth    = MISSING;
  global.startday      = MISSING;
  global.starthour     = MISSING;
  global.endyear       = MISSING;
  global.endmonth      = MISSING;
  global.endday        = MISSING;
  global.resolution    = MISSING;
  global.dt = grads_file->dt;//1; //hours
  global.nrecs = -99999; //Number of time steps
  global.startyear = grads_file->year; //Start year
  global.startmonth = grads_file->month;//Start month
  global.startday = grads_file->day;//Start day
  global.starthour = 0;//Start hour
  global.endyear = 2010;//End year
  global.endmonth = 12;//End month
  global.endday = 31;//End day
  global.resolution = 1.0;//Resolution
  global.MAX_SNOW_TEMP = 0.5;
  global.MIN_RAIN_TEMP = -0.5;
  global.measure_h     = 2.0;
  global.wind_h        = 10.0;
  for(i = 0; i < 2; i++) {
    global.forceyear[i]  = MISSING;
    global.forcemonth[i] = 1;
    global.forceday[i]   = 1;
    global.forcehour[i]  = 0;
    global.forceskip[i]  = 0;
    strcpy(names->f_path_pfx[i],"MISSING");
  }
  file_num             = 0;
  global.skipyear      = 0;
  strcpy(names->init_state,   "MISSING");
  global.stateyear     = MISSING;
  global.statemonth    = MISSING;
  global.stateday      = MISSING;
  strcpy(names->statefile,    "MISSING");
  //strcpy(names->soil,         "MISSING");
  strcpy(names->soil_dir,     "MISSING");
  //strcpy(names->veg,          "MISSING");
  //strcpy(names->veglib,       "MISSING");
  //strcpy(names->snowband,     "MISSING");
  strcpy(names->lakeparam,    "MISSING");
  strcpy(names->result_dir,   "MISSING");
  global.out_dt        = MISSING;


      /*************************************
       Define Global Parameters
      *************************************/
      options.Nlayer = 3;
      options.Nnode = 20;//20;
      global.dt = grads_file->dt;//1; //hours
      options.SNOW_STEP = grads_file->dt;//1;
      options.FULL_ENERGY=TRUE;      
      options.FROZEN_SOIL=FALSE;//TRUE;
      options.QUICK_FLUX=FALSE;
      options.QUICK_SOLVE=TRUE; 
      options.NOFLUX=TRUE;
      options.IMPLICIT=FALSE;
      options.EXP_TRANS=FALSE;
      options.DIST_PRCP = TRUE;
      options.PREC_EXPT = 0.6;
      options.CORRPREC=FALSE;
      options.MIN_WIND_SPEED = 0.1;
      global.MIN_RAIN_TEMP = -0.5;
      global.MAX_SNOW_TEMP = 0.5;
      options.CONTINUEONERROR=FALSE;
      options.COMPUTE_TREELINE=FALSE;
      options.EQUAL_AREA=FALSE;
      global.resolution = 1.0;
      options.PLAPSE = TRUE;

      /*************************************
       Define forcing files
      *************************************/
      param_set.FORCE_FORMAT[0] = ASCII;
      param_set.FORCE_ENDIAN[0] = LITTLE;
      param_set.N_TYPES[0] = 7;
      param_set.FORCE_DT[0] = grads_file->dt;
      global.forceyear[0] = grads_file->year;//2000;
      global.forcemonth[0] = grads_file->month;//1;
      global.forceday[0] = grads_file->day;//1;
      global.forcehour[0] = 0;
      options.GRID_DECIMAL = 2;
      global.wind_h = 2.0;
      global.measure_h = 2.0;
      strcpy(names->f_path_pfx[0],"/home/ice/nchaney/PROJECTS/BLUE_WATERS/TEST/DATA/Input/forcings_");
      //strcpy(names->soil,"/home/ice/nchaney/PROJECTS/BLUE_WATERS/TEST/DATA/Soil.txt");
      //strcpy(names->veglib,"/home/ice/nchaney/PROJECTS/BLUE_WATERS/TEST/DATA/veglib.dat");
      //strcpy(names->veg,"/home/ice/nchaney/PROJECTS/BLUE_WATERS/TEST/DATA/LAI.txt");
      //options.SNOW_BAND = 10;
      options.VEGPARAM_LAI=TRUE;
      options.LAI_SRC = LAI_FROM_VEGPARAM;
      options.ROOT_ZONES = 2;//3;
      strcpy(names->result_dir,"/home/ice/nchaney/PROJECTS/BLUE_WATERS/TEST/DATA/Output"); 
      global.out_dt = 0;
      global.skipyear = 0;
      options.COMPRESS = FALSE;
      options.BINARY_OUTPUT = FALSE;
      field = 0;
      strcpy(cmdstr,"FORCE_TYPE PREC");
      get_force_type(cmdstr,file_num,&field);
      strcpy(cmdstr,"FORCE_TYPE AIR_TEMP");
      get_force_type(cmdstr,file_num,&field);
      strcpy(cmdstr,"FORCE_TYPE WIND");
      get_force_type(cmdstr,file_num,&field);
      //strcpy(cmdstr,"FORCE_TYPE REL_HUMID");
      strcpy(cmdstr,"FORCE_TYPE QAIR");
      get_force_type(cmdstr,file_num,&field);
      strcpy(cmdstr,"FORCE_TYPE SHORTWAVE");
      get_force_type(cmdstr,file_num,&field);
      strcpy(cmdstr,"FORCE_TYPE LONGWAVE");
      get_force_type(cmdstr,file_num,&field);
      strcpy(cmdstr,"FORCE_TYPE PRESSURE");
      get_force_type(cmdstr,file_num,&field);

  /******************************************
    Check for undefined required parameters
  ******************************************/

  // Validate model time step
  if (global.dt == MISSING)
    nrerror("Model time step has not been defined.  Make sure that the global file defines TIME_STEP.");
  else if (global.dt < 1) {
    sprintf(ErrStr,"The specified model time step (%d) < 1 hour.  Make sure that the global file defines a positive number of hours for TIME_STEP.",global.dt);
    nrerror(ErrStr);
  }

  // Validate the output step
  if (global.out_dt == 0 || global.out_dt == MISSING) {
    global.out_dt = global.dt;
  }
  else if (global.out_dt < global.dt || global.out_dt > 24 || (float)global.out_dt/(float)global.dt != (float)(global.out_dt/global.dt)){
    nrerror("Invalid output step specified.  Output step must be an integer multiple of the model time step; >= model time step and <= 24");
  }

  // Validate SNOW_STEP and set NR and NF
  if (global.dt < 24 && global.dt != options.SNOW_STEP)
    nrerror("If the model step is smaller than daily, the snow model should run\nat the same time step as the rest of the model.");
  if (global.dt % options.SNOW_STEP != 0 || options.SNOW_STEP > global.dt)
    nrerror("SNOW_STEP should be <= TIME_STEP and divide TIME_STEP evenly ");
  NF = global.dt/options.SNOW_STEP;
  if (NF == 1)
    NR = 0;
  else
    NR = NF;

  // Validate simulation start date
  if (global.startyear == MISSING)
    nrerror("Simulation start year has not been defined.  Make sure that the global file defines STARTYEAR.");
  else if (global.startyear < 0) {
    sprintf(ErrStr,"The specified simulation start year (%d) < 0.  Make sure that the global file defines a positive integer for STARTYEAR.",global.startyear);
    nrerror(ErrStr);
  }
  if (global.startmonth == MISSING)
    nrerror("Simulation start month has not been defined.  Make sure that the global file defines STARTMONTH.");
  else if (global.startmonth < 0) {
    sprintf(ErrStr,"The specified simulation start month (%d) < 0.  Make sure that the global file defines a positive integer for STARTMONTH.",global.startmonth);
    nrerror(ErrStr);
  }
  if (global.startday == MISSING)
    nrerror("Simulation start day has not been defined.  Make sure that the global file defines STARTDAY.");
  else if (global.startday < 0) {
    sprintf(ErrStr,"The specified simulation start day (%d) < 0.  Make sure that the global file defines a positive integer for STARTDAY.",global.startday);
    nrerror(ErrStr);
  }
  if (global.starthour == MISSING) {
    if (global.dt == 24)
      global.starthour = 0;
    else
      nrerror("Simulation start hour has not been defined, yet model time step is less than 24 hours.  Make sure that the global file defines STARTHOUR.");
  }
  else if (global.starthour < 0) {
    sprintf(ErrStr,"The specified simulation start hour (%d) < 0.  Make sure that the global file defines a positive integer for STARTHOUR.",global.starthour);
    nrerror(ErrStr);
  }

  // Validate simulation end date and/or number of timesteps
  if (global.nrecs == MISSING && global.endyear == MISSING && global.endmonth == MISSING && global.endday == MISSING)
    nrerror("The model global file MUST define EITHER the number of records to simulate (NRECS), or the year (ENDYEAR), month (ENDMONTH), and day (ENDDAY) of the last full simulation day");
  else if (global.nrecs == MISSING) {
    if (global.endyear == MISSING)
      nrerror("Simulation end year has not been defined.  Make sure that the global file defines ENDYEAR.");
    else if (global.endyear < 0) {
      sprintf(ErrStr,"The specified simulation end year (%d) < 0.  Make sure that the global file defines a positive integer for ENDYEAR.",global.endyear);
      nrerror(ErrStr);
    }
    if (global.endmonth == MISSING)
      nrerror("Simulation end month has not been defined.  Make sure that the global file defines ENDMONTH.");
    else if (global.endmonth < 0) {
      sprintf(ErrStr,"The specified simulation end month (%d) < 0.  Make sure that the global file defines a positive integer for ENDMONTH.",global.endmonth);
      nrerror(ErrStr);
    }
    if (global.endday == MISSING)
      nrerror("Simulation end day has not been defined.  Make sure that the global file defines ENDDAY.");
    else if (global.endday < 0) {
      sprintf(ErrStr,"The specified simulation end day (%d) < 0.  Make sure that the global file defines a positive integer for ENDDAY.",global.endday);
      nrerror(ErrStr);
    }
    tmpstartdate = global.startyear*10000 + global.startmonth*100 + global.startday;
    tmpenddate = global.endyear*10000 + global.endmonth*100 + global.endday;
    if (tmpenddate < tmpstartdate) {
      sprintf(ErrStr,"The specified simulation end date (%04d-%02d-%02d) is EARLIER than the specified start date (%04d-%02d-%02d).",global.endyear,global.endmonth,global.endday,global.startyear,global.startmonth,global.startday);
      nrerror(ErrStr);
    }
  }
  else if (global.nrecs < 1) {
    sprintf(ErrStr,"The specified duration of simulation (%d) < 1 time step.  Make sure that the global file defines a positive integer for NRECS.",global.nrecs);
    nrerror(ErrStr);
  }

  // Validate forcing files and variables
  if ( strcmp ( names->f_path_pfx[0], "MISSING" ) == 0 )
    nrerror("No forcing file has been defined.  Make sure that the global file defines FORCING1.");
  for(i=0;i<2;i++) {
    if ( i == 0 || (i == 1 && param_set.N_TYPES[i] != MISSING) ) {
      if (param_set.N_TYPES[i] == MISSING) {
        sprintf(ErrStr,"Need to specify the number forcing variables types in forcing file %d.", i);
        nrerror(ErrStr);
      }
      if (param_set.FORCE_FORMAT[i] == MISSING) {
        sprintf(ErrStr,"Need to specify the INPUT_FORMAT (ASCII or BINARY) for forcing file %d.",i);
        nrerror(ErrStr);
      }
      if (param_set.FORCE_INDEX[i][param_set.N_TYPES[i]-1] == MISSING) {
        sprintf(ErrStr,"Did not define enough forcing variables in forcing file %d.",i);
        nrerror(ErrStr);
      }
      if(param_set.FORCE_DT[i] == MISSING ) {
        sprintf(ErrStr,"Must define time steps (FORCE_DT <dt>) in control file for focing file %d.",file_num);
        nrerror(ErrStr);
      }
    }
  }
  if(param_set.N_TYPES[1] != MISSING && global.forceyear[1] == MISSING) {
    global.forceyear[1] = global.forceyear[0];
    global.forcemonth[1] = global.forcemonth[0];
    global.forceday[1] = global.forceday[0];
    global.forcehour[1] = global.forcehour[0];
    global.forceskip[1] = 0;
  }

  // Validate result directory
  if ( strcmp ( names->result_dir, "MISSING" ) == 0 )
    nrerror("No results directory has been defined.  Make sure that the global file defines the result directory on the line that begins with \"RESULT_DIR\".");

  // Validate soil parameter file information
  //if ( strcmp ( names->soil, "MISSING" ) == 0 )
  //  nrerror("No soil parameter file has been defined.  Make sure that the global file defines the soil parameter file on the line that begins with \"SOIL\".");
  if (options.ARC_SOIL && strcmp ( names->soil_dir, "MISSING" ) == 0)
    nrerror("\"ARC_SOIL\" was specified as TRUE, but no soil parameter directory (\"SOIL_DIR\") has been defined.  Make sure that the global file defines the soil parameter directory on the line that begins with \"SOIL_DIR\".");

  /*******************************************************************************
    Validate parameters required for normal simulations but NOT for OUTPUT_FORCE
  *******************************************************************************/

#if !OUTPUT_FORCE

  // Validate veg parameter information
  //if ( strcmp ( names->veg, "MISSING" ) == 0 )
  //  nrerror("No vegetation parameter file has been defined.  Make sure that the global file defines the vegetation parameter file on the line that begins with \"VEGPARAM\".");
  //if ( strcmp ( names->veglib, "MISSING" ) == 0 )
  //  nrerror("No vegetation library file has been defined.  Make sure that the global file defines the vegetation library file on the line that begins with \"VEGLIB\".");
  if(options.ROOT_ZONES<0)
    nrerror("ROOT_ZONES must be defined to a positive integer greater than 0, in the global control file.");
  if (options.LAI_SRC == LAI_FROM_VEGPARAM && !options.VEGPARAM_LAI) {
      sprintf(ErrStr, "\"LAI_SRC\" was specified as \"LAI_FROM_VEGPARAM\", but \"VEGPARAM_LAI\" was set to \"FALSE\" in the global parameter file.  If you want VIC to read LAI values from the vegparam file, you MUST make sure the veg param file contains 1 line of 12 monthly LAI values for EACH veg tile in EACH grid cell, and you MUST specify \"VEGPARAM_LAI\" as \"TRUE\" in the global parameter file.  Alternatively, if you want VIC to read LAI values from the veg library file, set \"LAI_SRC\" ro \"LAI_FROM_VEGLIB\" in the global parameter file.  In either case, the setting of \"VEGPARAM_LAI\" must be consistent with the contents of the veg param file (i.e. whether or not it contains LAI values).");
      nrerror(ErrStr);
  }

  // Validate the elevation band file information
  if(options.SNOW_BAND > 1) {
    if ( strcmp ( names->snowband, "MISSING" ) == 0 ) {
      sprintf(ErrStr, "\"SNOW_BAND\" was specified with %d elevation bands, but no elevation band file has been defined.  Make sure that the global file defines the elevation band file on the line that begins with \"SNOW_BAND\" (after the number of bands).", options.SNOW_BAND);
      nrerror(ErrStr);
    }
    if(options.SNOW_BAND > MAX_BANDS) {
      sprintf(ErrStr,"Global file wants more snow bands (%d) than are defined by MAX_BANDS (%d).  Edit user_def.h and recompile.",options.SNOW_BAND,MAX_BANDS);
      nrerror(ErrStr);
    }
  }
  else if (options.SNOW_BAND <= 0) {
    sprintf(ErrStr,"Invalid number of elevation bands specified in global file (%d).  Number of bands must be >= 1.",options.SNOW_BAND);
    nrerror(ErrStr);
  }

  // Validate the input state file information
  if( options.INIT_STATE ) {
    if ( strcmp ( names->init_state, "MISSING" ) == 0 )
      nrerror("\"INIT_STATE\" was specified, but no input state file has been defined.  Make sure that the global file defines the inputstate file on the line that begins with \"INIT_STATE\".");
  }

  // Validate the output state file information
  if( options.SAVE_STATE ) {
    if ( strcmp ( names->statefile, "MISSING" ) == 0)
      nrerror("\"SAVE_STATE\" was specified, but no output state file has been defined.  Make sure that the global file defines the output state file on the line that begins with \"SAVE_STATE\".");
    if ( global.stateyear == MISSING || global.statemonth == MISSING || global.stateday == MISSING )  {
      sprintf(ErrStr,"Incomplete specification of the date to save state for state file (%s).\nSpecified date (yyyy-mm-dd): %04d-%02d-%02d\nMake sure STATEYEAR, STATEMONTH, and STATEDAY are set correctly in your global parameter file.\n", names->statefile, global.stateyear, global.statemonth, global.stateday);
      nrerror(ErrStr);
    }
    // Check for month, day in range
    lastvalidday = lastday[global.statemonth - 1];
    if ( global.statemonth == 2 ) {
      if ( (global.stateyear % 4) == 0 && ( (global.stateyear % 100) != 0 || (global.stateyear % 400) == 0 ) ){
        lastvalidday = 29;
      }
    }
    if ( global.stateday > lastvalidday || global.statemonth > 12 || global.statemonth < 1 || global.stateday > 31 || global.stateday < 1 ){
      sprintf(ErrStr,"Unusual specification of the date to save state for state file (%s).\nSpecified date (yyyy-mm-dd): %04d-%02d-%02d\nMake sure STATEYEAR, STATEMONTH, and STATEDAY are set correctly in your global parameter file.\n", names->statefile, global.stateyear, global.statemonth, global.stateday);
      nrerror(ErrStr);
    }
  }
  // Set the statename here to be able to compare with INIT_STATE name
  if( options.SAVE_STATE ) {
    sprintf(names->statefile,"%s_%04i%02i%02i", names->statefile,
          global.stateyear, global.statemonth, global.stateday);
  }
  if( options.INIT_STATE && options.SAVE_STATE && (strcmp( names->init_state, names->statefile ) == 0))  {
      sprintf(ErrStr,"The save state file (%s) has the same name as the initialize state file (%s).  The initialize state file will be destroyed when the save state file is opened.", names->statefile, names->init_state);
      nrerror(ErrStr);
  }

  // Validate soil parameter/simulation mode combinations
  if(options.QUICK_FLUX) {
    if(options.Nnode != 3) {
      fprintf(stderr,"WARNING: To run the model QUICK_FLUX=TRUE, you must define exactly 3 soil thermal nodes.  Currently Nnodes is set to %d.  Setting Nnodes to 3.\n",options.Nnode);
      options.Nnode = 3;
    }
    if(options.IMPLICIT || options.EXP_TRANS) {
      sprintf(ErrStr,"To run the model with QUICK_FLUX=TRUE, you cannot have IMPLICIT=TRUE or EXP_TRANS=TRUE.");
      nrerror(ErrStr);
    }
  }
  else {
    if(!options.FULL_ENERGY && !options.FROZEN_SOIL) {
      sprintf(ErrStr,"To run the model in water balance mode (both FULL_ENERGY and FROZEN_SOIL are FALSE) you MUST set QUICK_FLUX to TRUE (or leave QUICK_FLUX out of your global parameter file).");
      nrerror(ErrStr);
    }
  }
  if((options.FULL_ENERGY || options.FROZEN_SOIL) && options.Nlayer<3) {
    sprintf(ErrStr,"You must define at least 3 soil moisture layers to run the model in FULL_ENERGY or FROZEN_SOIL modes.  Currently Nlayers is set to  %d.",options.Nlayer);
    nrerror(ErrStr);
  }
  if((!options.FULL_ENERGY && !options.FROZEN_SOIL) && options.Nlayer<1) {
    sprintf(ErrStr,"You must define at least 1 soil moisture layer to run the model.  Currently Nlayers is set to  %d.",options.Nlayer);
    nrerror(ErrStr);
  }
  if(options.IMPLICIT)  {
    if ( QUICK_FS ) 
      fprintf(stderr,"WARNING: IMPLICIT and QUICK_FS are both TRUE.\n\tThe QUICK_FS option is ignored when IMPLICIT=TRUE\n");
  }
  if( EXCESS_ICE ) {
    if ( !options.FULL_ENERGY )
      nrerror("set FULL_ENERGY = TRUE to run EXCESS_ICE option.");
    if ( !options.FROZEN_SOIL )
      nrerror("set FROZEN_SOIL = TRUE to run EXCESS_ICE option.");
    if ( options.QUICK_SOLVE ) {
      fprintf(stderr,"WARNING: QUICK_SOLVE and EXCESS_ICE are both TRUE.\n\tThis is an incompatible combination.  Setting QUICK_SOLVE to FALSE.\n");
      options.QUICK_SOLVE=FALSE;  
    }    
    if ( QUICK_FS ) 
      nrerror("QUICK_FS = TRUE and EXCESS_ICE = TRUE are incompatible options.");
  }
  if(options.Nlayer > MAX_LAYERS) {
    sprintf(ErrStr,"Global file wants more soil moisture layers (%d) than are defined by MAX_LAYERS (%d).  Edit user_def.h and recompile.",options.Nlayer,MAX_LAYERS);
    nrerror(ErrStr);
  }
  if(options.Nnode > MAX_NODES) {
    sprintf(ErrStr,"Global file wants more soil thermal nodes (%d) than are defined by MAX_NODES (%d).  Edit user_def.h and recompile.",options.Nnode,MAX_NODES);
    nrerror(ErrStr);
  }

  // Validate lake parameter information
  if (options.LAKES) {
    if (!options.FULL_ENERGY) {
      sprintf(ErrStr, "FULL_ENERGY must be TRUE if the lake model is to be run.");
      nrerror(ErrStr);
    }
    if ( strcmp ( names->lakeparam, "MISSING" ) == 0 )
      nrerror("\"LAKES\" was specified, but no lake parameter file has been defined.  Make sure that the global file defines the lake parameter file on the line that begins with \"LAKES\".");
    if (global.resolution == 0) {
      sprintf(ErrStr, "The model grid cell resolution (RESOLUTION) must be defined in the global control file when the lake model is active.");
      nrerror(ErrStr);
    }
    if (global.resolution > 360 && !options.EQUAL_AREA) {
      sprintf(ErrStr, "For EQUAL_AREA=FALSE, the model grid cell resolution (RESOLUTION) must be set to the number of lat or lon degrees per grid cell.  This cannot exceed 360.");
      nrerror(ErrStr);
    }
    if (options.COMPUTE_TREELINE) {
      sprintf(ErrStr, "LAKES = TRUE and COMPUTE_TREELINE = TRUE are incompatible options.");
      nrerror(ErrStr);
    }
  }

  /*********************************
    Output major options to stderr
  *********************************/
#if VERBOSE
  display_current_settings(DISP_ALL,names,&global);
#else
  display_current_settings(DISP_VERSION,names,&global);
#endif

#if VERBOSE
  fprintf(stderr,"Time Step = %d hour(s)\n",global.dt);
  fprintf(stderr,"Simulation start date = %02i/%02i/%04i\n",
	  global.startday, global.startmonth, global.startyear);
  if ( global.nrecs > 0 )
    fprintf(stderr,"Number of Records = %d\n\n",global.nrecs);
  else 
    fprintf(stderr,"Simulation end date = %02i/%02i/%04i\n\n",
	    global.endday, global.endmonth, global.endyear);
  fprintf(stderr,"Full Energy...................(%d)\n",options.FULL_ENERGY);
  fprintf(stderr,"Use Distributed Precipitation.(%d)\n",options.DIST_PRCP);
  if(options.DIST_PRCP)
    fprintf(stderr,"..Using Precipitation Exponent of %f\n",options.PREC_EXPT);
  fprintf(stderr,"Ground heat flux will be estimated ");
  if ( options.QUICK_FLUX ) 
    fprintf(stderr,"using Liang, Wood and Lettenmaier (1999).\n");
  else 
    fprintf(stderr,"using Cherkauer and Lettenmaier (1999).\n");
  fprintf(stderr,"Use Frozen Soil Model.........(%d)\n",options.FROZEN_SOIL);
  if( options.IMPLICIT ) 
    fprintf(stderr,".... Using the implicit solution for the soil heat equation.\n");
  else
    fprintf(stderr,".... Using the explicit solution for the soil heat equation.\n");
  if( options.EXP_TRANS )
    fprintf(stderr,".... Thermal nodes are exponentially distributed with depth.\n");
  else
    fprintf(stderr,".... Thermal nodes are linearly distributed with depth (except top two nodes).\n");
  if( EXCESS_ICE )
    fprintf(stderr,".... Excess ground ice is being considered.\n\t\tTherefore, ground ice (as a volumetric fraction) must be initialized for each\n\t\t   soil layer in the soil file.\n\t\tCAUTION: When excess ice melts, subsidence occurs.\n\t\t  Therefore, soil layer depths, damping depth, thermal node depths,\n\t\t     bulk densities, porosities, and other properties are now dynamic!\n\t\t  EXERCISE EXTREME CAUTION IN INTERPRETING MODEL OUTPUT.\n\t\t  It is recommended to add OUT_SOIL_DEPTH to your list of output variables.\n");
  if ( QUICK_FS ){
    fprintf(stderr,".... Using linearized UFWC curve with %d temperatures.\n", QUICK_FS_TEMPS);
  }
  fprintf(stderr,"Run Snow Model Using a Time Step of %d hours\n", 
	  options.SNOW_STEP);
  fprintf(stderr,"Compress Output Files.........(%d)\n",options.COMPRESS);
  fprintf(stderr,"Correct Precipitation.........(%d)\n",options.CORRPREC);
  fprintf(stderr,"\n");
  fprintf(stderr,"Using %d Snow Bands\n",options.SNOW_BAND);
  fprintf(stderr,"Using %d Root Zones\n",options.ROOT_ZONES);
  if ( options.SAVE_STATE )
    fprintf(stderr,"Model state will be saved on = %02i/%02i/%04i\n\n",
	    global.stateday, global.statemonth, global.stateyear);
  if ( options.BINARY_OUTPUT ) 
    fprintf(stderr,"Model output is in standard BINARY format.\n");
  else 
    fprintf(stderr,"Model output is in standard ASCII format.\n");

#endif // VERBOSE

#endif // !OUTPUT_FORCE

  return global;

}
