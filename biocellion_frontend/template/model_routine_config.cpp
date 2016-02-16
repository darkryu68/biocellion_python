/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */
#include "biocellion.h"
#include "model_routine.h"

/* MODEL START */

#include "model_define.h"


/* MODEL END */

S32  INI_N_CELLS;

using namespace std;

void ModelRoutine::updateIfGridSpacing( REAL& ifGridSpacing ) {
	/* MODEL START */

	ifGridSpacing = IF_GRID_SPACING;

	/* MODEL END */

	return;
}

void ModelRoutine::updateOptModelRoutineCallInfo( OptModelRoutineCallInfo& callInfo ) {
	/* MODEL START */

	callInfo.numUpdateIfGridVarPreStateAndGridStepRounds = 1;
	callInfo.numUpdateIfGridVarPostStateAndGridStepRounds = 0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateDomainBdryType( domain_bdry_type_e a_domainBdryType[DIMENSION] ) {
	/* MODEL START */

	a_domainBdryType[0] = DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL;
	a_domainBdryType[1] = DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL;
	a_domainBdryType[2] = DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL;

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferBdryType( pde_buffer_bdry_type_e& pdeBufferBdryType ) {
	/* MODEL START */

	pdeBufferBdryType = PDE_BUFFER_BDRY_TYPE_HARD_WALL;

	/* MODEL END */

	return;
}

void ModelRoutine::updateTimeStepInfo( TimeStepInfo& timeStepInfo ) {
	/* MODEL START */

	timeStepInfo.durationBaselineTimeStep = BASELINE_TIME_STEP_DURATION;
	timeStepInfo.numStateAndGridTimeStepsPerBaseline = NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE;

	/* MODEL END */

	return;
}

void ModelRoutine::updateSyncMethod( sync_method_e& extraMechIntrctSyncMethod, sync_method_e& updateIfGridVarSyncMethod/* dummy if both callUpdateIfGridVarPreStateAndGridStep and callUpdateIfGridVarPostStateAndGridStep are set to false in ModelRoutine::updateOptModelRoutineCallInfo */ ) {
	/* MODEL START */

	extraMechIntrctSyncMethod = SYNC_METHOD_OVERWRITE;
	updateIfGridVarSyncMethod = SYNC_METHOD_OVERWRITE;

	/* MODEL END */

	return;
}


void ModelRoutine::updateJunctionEndInfo( Vector<JunctionEndInfo>& v_junctionEndInfo ) {
	/* MODEL START */

	v_junctionEndInfo.resize( NUM_JUNCTION_END_TYPES );
	for( S32 i = 0 ; i < NUM_JUNCTION_END_TYPES ; i++ ) {
		v_junctionEndInfo[i].numModelReals = 0;/* we do not associate any variable to a junction end type */
		v_junctionEndInfo[i].numModelInts = 0;
	} 

	/* MODEL END */

	return;
}


void ModelRoutine::updateRNGInfo( Vector<RNGInfo>& v_rngInfo ) {
	/* MODEL START */

	CHECK( NUM_MODEL_RNGS == 2 );

	v_rngInfo.resize( 2 );

	RNGInfo rngInfo;

	rngInfo.type = RNG_TYPE_UNIFORM;
	rngInfo.param0 = 0.0;
	rngInfo.param1 = 1.0;
	rngInfo.param2 = 0.0;/* dummy */
	v_rngInfo[MODEL_RNG_UNIFORM] = rngInfo;

        rngInfo.type = RNG_TYPE_UNIFORM;
        rngInfo.param0 = 0.9;
        rngInfo.param1 = 1.1;
        rngInfo.param2 = 0.0;/* dummy */
        v_rngInfo[MODEL_RNG_UNIFORM_10PERCENT] = rngInfo ;

	/* MODEL END */

	return;
}

void ModelRoutine::updateFileOutputInfo( FileOutputInfo& fileOutputInfo ) {
	/* MODEL START */

	fileOutputInfo.particleOutput = true;
	fileOutputInfo.particleNumExtraOutputVars = 3;
	fileOutputInfo.v_gridPhiOutput.assign( NUM_DIFFUSIBLE_ELEMS, true );
	fileOutputInfo.v_gridPhiOutputDivideByKappa.assign( NUM_DIFFUSIBLE_ELEMS, false );

	/* MODEL END */

	return;
}

void ModelRoutine::updateGlobalData( Vector<U8>& v_globalData ) {
   /* MODEL START */

   Vector<string> v_modelParam;
   S8* p_modelParams;
   S8* p_param;

   UBInitData* p_ubInitData;
   IniCellData* p_IniCellData ;

   FILE* p_file;
   S8* p_buf;
   size_t bufSize;
   ssize_t numCharsRead;

   S32 IniNumberCells = 0;
   S32 NumberVoxels = 0;

   VIdx ifRegionSize;
   ifRegionSize[0] = Info::getDomainSize(0) - AGAR_HEIGHT ;
   ifRegionSize[1] = Info::getDomainSize(1) ;
   ifRegionSize[2] = Info::getDomainSize(2) ;

   NumberVoxels = ifRegionSize[0] * ifRegionSize[1] * ifRegionSize[2] ;

   Vector<VIdx> v_habitableVIdxOffset;

   p_modelParams = new S8[Info::getModelParam().size() + 1 ]; // + 1 for null termination
   strcpy( p_modelParams, Info::getModelParam().c_str() );

   p_param = strtok( p_modelParams, " \t\n" );
   while( p_param != NULL ) {
      v_modelParam.push_back( string( p_param ) );
      p_param = strtok( NULL, " \t\n" );
   }

   delete[] p_modelParams;

   if( v_modelParam.size() != 1 ) {
           ERROR( "invalid model parameters: invalid number of parameters." );
   }
 
   // read Boxes with yeast cells 
   p_file = fopen( v_modelParam[0].c_str(), "r" );
   if( p_file == NULL ) {
       ERROR( "invalid model parameters: invalid blocked_grid_points file path." );
   }

   p_buf = NULL;
   bufSize = 0;
   S32 count = 0;
   if (( numCharsRead = getline( &p_buf, &bufSize, p_file ) ) != -1 ) {
      S8* p_token;
      p_token = strtok( p_buf, " \t\n" );
      IniNumberCells = (S32) strtol( p_token, NULL, 10 );
   }
   else {
      ERROR( "Unable to read number of cells" );
   }
   INI_N_CELLS = IniNumberCells ;

   // init global data 
   v_globalData.resize( NumberVoxels*sizeof(UBInitData) + IniNumberCells*sizeof(IniCellData));
   p_ubInitData = ( UBInitData* )&( v_globalData[ 0  ] );
   p_IniCellData = ( IniCellData * )&( v_globalData[ NumberVoxels*sizeof(UBInitData) ] );
   // initialize the  the data (p_ubInitData),  X_SIZE * Z_SIZE  uvInitDatas
   for( idx_t i = 0 ; i <  ifRegionSize[0] ; i++ ) {
      for( idx_t j = 0 ; j < ifRegionSize[1]  ; j++ ) {
         for( idx_t k = 0 ; k < ifRegionSize[2]  ; k++ ) {
            UBInitData& ubInitData = p_ubInitData[VIdx::getIdx3DTo1D( i, j, k, ifRegionSize )];
            ubInitData.numCells = 0;
            ubInitData.x_offset = 0.0 ;
            ubInitData.y_offset = 0.0 ;
            ubInitData.z_offset = 0.0 ;
            ubInitData.a_type = 0 ;
            ubInitData.biomass = 0.0 ;
            ubInitData.inert = 0.0 ;
            ubInitData.IdxIniCellData  = -1 ;
         }
      }
   }
   for ( S32 i=0; i < IniNumberCells; i++){
      IniCellData& cellInitData = p_IniCellData[i];
      cellInitData.x_offset = 0.0 ;
      cellInitData.y_offset = 0.0 ;
      cellInitData.z_offset = 0.0 ;
      cellInitData.a_type = 0 ;
      cellInitData.biomass = 0.0 ;
      cellInitData.inert = 0.0 ;
      cellInitData.IdxIniCellData  = -1 ;
   }

   while( ( numCharsRead = getline( &p_buf, &bufSize, p_file ) ) != -1 ) {

      VIdx vIdxOffset;
      VReal posOffset;


      S8* p_token;
      p_token = strtok( p_buf, " \t\n" );

      for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {
         if( p_token == NULL ) {
              ERROR( "invalid blocked_grid_points file." );
         }
         REAL Position  = ( REAL )strtod( p_token, NULL );
         vIdxOffset[dim] = ( idx_t ) ( Position/IF_GRID_SPACING) ;
         posOffset[dim] = Position - (( (REAL)vIdxOffset[dim]*IF_GRID_SPACING ) + 0.5 * IF_GRID_SPACING ) ;

         if( posOffset[dim] < IF_GRID_SPACING * -0.5 )
             posOffset[dim] = IF_GRID_SPACING  * -0.5;
         else if( posOffset[dim] > IF_GRID_SPACING  * 0.5 )
             posOffset[dim] = IF_GRID_SPACING * 0.5;

         p_token = strtok( NULL, " \t\n" );

         if( ( vIdxOffset[dim] < 0 )   ) {
             ERROR( "invalid blocked_grid_points file." );
         }
      }
     
      UBInitData& ubInitData = p_ubInitData[VIdx::getIdx3DTo1D( vIdxOffset, ifRegionSize )];
      IniCellData& cellInitData = p_IniCellData[count];

      if( p_token == NULL ) {
          ERROR( "invalid blocked_grid_points file." );
      }
      else {
          // ubInitData.radius  = strtol( p_token, NULL, 10 )  
          cellInitData.biomass  = ( REAL )strtod( p_token, NULL );
      }

      p_token = strtok( NULL, " \t\n" );
      if( p_token == NULL ) {
          ERROR( "invalid blocked_grid_points file." );
      }
      else {
          cellInitData.inert  =( REAL)strtod( p_token, NULL ) ;
      } 

      p_token = strtok( NULL, " \t\n" );
      if( p_token == NULL ) {
          ERROR( "invalid blocked_grid_points file." );
      }
      else {
          cellInitData.a_type  = strtol( p_token, NULL, 10 ) ;
      }

      p_token = strtok( NULL, " \t\n" );
      if( p_token != NULL ) {
          ERROR( "invalid blocked_grid_points file." );
      }

      cellInitData.x_offset = (REAL) posOffset[0] ;
      cellInitData.y_offset = (REAL) posOffset[1] ;
      cellInitData.z_offset = (REAL)  posOffset[2] ;

      if ( ubInitData.numCells == 0 ) {
          ubInitData.IdxIniCellData  = count ; //  p_IniCellData[count];
      }
      else{
         S32 nextIdx = ubInitData.IdxIniCellData;
         if ( nextIdx == -1 )
             ERROR( "invalid blocked_grid_points file." );

         for ( S32 i = 1; i < ubInitData.numCells; i++) {
              if ( nextIdx == -1 )
                 ERROR( "invalid blocked_grid_points file." );

              IniCellData& cellData = p_IniCellData[ nextIdx  ];
              nextIdx = cellData.IdxIniCellData  ;          
         }

         IniCellData& cellData = p_IniCellData[ nextIdx  ];
         if ( cellData.IdxIniCellData != -1 )
              ERROR( "invalid blocked_grid_points file." );
         cellData.IdxIniCellData = count;
      }
      ubInitData.numCells++ ;
      count++ ;
   }

   if( p_buf != NULL ) {
      free( p_buf );
   }
   fclose( p_file );

   /* MODEL END */

   return;
}

void ModelRoutine::init( void ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

void ModelRoutine::term( void ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

void ModelRoutine::setPDEBuffer( const VIdx& startVIdx, const VIdx& regionSize, BOOL& isPDEBuffer ) {
	/* MODEL START */

	if( (startVIdx[0] + regionSize[0]) >= AGAR_HEIGHT ) {
		isPDEBuffer = false;
	}
	else {
		isPDEBuffer = true;
	}

	/* MODEL END */

	return;
}

void ModelRoutine::setHabitable( const VIdx& vIdx, BOOL& isHabitable ) {
	/* MODEL START */

	if( vIdx[0] >=  AGAR_HEIGHT ) {
		isHabitable = true;
	}
	else {/* currently, cells cannot penetrate agar */
		isHabitable = false;
	}

	/* MODEL END */

	return;
}

