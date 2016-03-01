/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

/* MODEL START */

#include "model_define.h"

/* MODEL END */

using namespace std;

#if HAS_SPAGENT
void ModelRoutine::addSpAgents( const BOOL init, const VIdx& startVIdx, const VIdx& regionSize, const IfGridBoxData<BOOL>& ifGridHabitableBoxData, Vector<VIdx>& v_spAgentVIdx, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentOffset ) {
   /* MODEL START */

   if( init == true ) {

       const Vector<U8>& v_globalData = Info::getGlobalDataRef(); // data passed to all mpi proce
       VIdx ifRegionSize;  // Help computing Index for data in global data
       const UBInitData* p_ubInitData; // pointer to the data
       const IniCellData* p_IniCellData; 

       ifRegionSize[0] = Info::getDomainSize(0) - AGAR_HEIGHT ;  // index wil 
       ifRegionSize[1] = Info::getDomainSize(1) ;
       ifRegionSize[2] = Info::getDomainSize(2) ;

       CHECK( v_globalData.size() == ( ifRegionSize[0]*ifRegionSize[1]*ifRegionSize[2]*sizeof( UBInitData )) + INI_N_CELLS*sizeof(IniCellData));

       //  Reading the address where data is located.
       //p_ubInitData = ( const UBInitData* )&( v_globalData[3*sizeof(S32) + 3*sizeof(REAL)] );
       p_ubInitData = ( const UBInitData* )&( v_globalData[0] );
       p_IniCellData = ( const IniCellData* )&( v_globalData[ ifRegionSize[0]*ifRegionSize[1]*ifRegionSize[2]*sizeof(UBInitData)] );

       // place cells according to input file (image)  
       for ( S32 i = startVIdx[0]; i < startVIdx[0] + regionSize[0] ; i++ ) {
         for ( S32 j = startVIdx[1]; j < startVIdx[1] + regionSize[1] ; j++ ) {
           for ( S32 k = startVIdx[2]; k < startVIdx[2] + regionSize[2] ; k++ ) {
               if (  i >= AGAR_HEIGHT  ) {
                   VIdx posVIdx;     // Voxel in the biocellion
                   VIdx vIdxOffset;  // Voxel in the input file
 
                   posVIdx[0] = i;
                   posVIdx[1] = j;
                   posVIdx[2] = k;

                   vIdxOffset = posVIdx;
                   vIdxOffset[0] -=  AGAR_HEIGHT   ;

                   const UBInitData& ubInitData = p_ubInitData[VIdx::getIdx3DTo1D( vIdxOffset, ifRegionSize )];
                   //cout << i << " "<< j << " " << k << " " <<  ubInitData.numCells << endl ; 
                   if ( ubInitData.numCells > 0 ) {
                      S32 IdxCell = ubInitData.IdxIniCellData ;
                      for ( S32 i = 0; i < ubInitData.numCells ; i++) { 
           
                         const IniCellData& cellInitData = p_IniCellData[IdxCell]; 
                      
                         SpAgentState state;
                         VReal posOffset;
                         CHECK( MAX_CELL_RADIUS >= MIN_CELL_RADIUS );

                         posOffset[0] = cellInitData.x_offset ;
                         posOffset[1] = cellInitData.y_offset ;
                         posOffset[2] = cellInitData.z_offset ;

                         state.setType( cellInitData.a_type );
                         REAL biomass  = cellInitData.biomass ; 
                         REAL inert =  cellInitData.inert ;
                         REAL total_biomass = biomass + inert;  
                  
		         REAL volume = total_biomass / A_DENSITY_BIOMASS[ubInitData.a_type]; 
                         REAL radius = radius_from_volume( volume ); 
   
                         if ( radius >  MAX_CELL_RADIUS  )
                             radius = MAX_CELL_RADIUS;
                         else if ( radius < MIN_CELL_RADIUS )
                             radius = MIN_CELL_RADIUS;
                       
                         state.setRadius( radius );
                         state.setModelReal( CELL_MODEL_REAL_BIOMAS, biomass );
                         state.setModelReal( CELL_MODEL_REAL_INERT, inert );
                         state.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, 1.0 );
                         state.setModelReal( CELL_MODEL_REAL_SECRETION_PCT, 1.0 );
                         state.setModelInt( CELL_MODEL_INT_BOND_B, 0 ) ;
   
                         Vector<REAL> v_odeNetVal; 
                         v_odeNetVal = state.getODEValArray( 0 );
                         if ( v_odeNetVal.size() > 0 ) {
                            for( S32 i = 0 ; i < ( S32 )v_odeNetVal.size() ; i++ ) {
                                v_odeNetVal[i] = 0.0;
                            }
                            if ( A_BIOMASS_ODE_INDEX[ cellInitData.a_type ] != -1 ) {
                               S32 ode_idx = A_BIOMASS_ODE_INDEX[cellInitData.a_type ] ;
                               v_odeNetVal[ode_idx] = biomass;
                            }  
                         }

                         CHECK( ifGridHabitableBoxData.get( posVIdx ) == true );
                         v_spAgentVIdx.push_back( posVIdx );
                         v_spAgentState.push_back( state );
                         v_spAgentOffset.push_back( posOffset ); 
                         //cout << "  Agent Addedd " << posOffset[0] <<" "<< posOffset[1] <<" "<< posOffset[2] <<  endl;

                         IdxCell = cellInitData.IdxIniCellData ; 
                      }

                   }
               }
           }
         }
       }

   }
   /* MODEL END */
   return;
}


void ModelRoutine::updateSpAgentState( const VIdx& vIdx, const AgentJunctionInfo& junctionInfo, const VReal& vOffset, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, SpAgentState& state/* INOUT */ ) {
    /* MODEL START */

    REAL aaa_ratio[3][3][3];

    Util::computeSphereUBVolOvlpRatio( SPHERE_UB_VOL_OVLP_RATIO_MAX_LEVEL, vOffset, state.getRadius(), aaa_ratio );

    REAL uptakePct = state.getModelReal( CELL_MODEL_REAL_UPTAKE_PCT );
    REAL secretionPct = state.getModelReal( CELL_MODEL_REAL_SECRETION_PCT );



    S32 type = state.getType() ; // id of cell type 
    //S32 pdeIdx  = A_AGENT_GROWTH_SOURCE[type] ;//index of the sol needed to update biomass

    /// this part can be encapsulated in another function
    if( uptakePct > 0.0 ) {/* live */
       REAL a_avg_diff[ NUM_DIFFUSIBLE_ELEMS  ]  ;
       for ( S32 i = 0; i < NUM_DIFFUSIBLE_ELEMS; i++ )
          a_avg_diff[ i  ] = 0.0;    

       if ( NUM_DIFFUSIBLE_ELEMS > 0 ) {

          for( S32 i = -1 ; i <= 1 ; i++ ) {
             for( S32 j = -1 ; j <= 1 ; j++ ) {
                for( S32 k = -1 ; k <= 1 ; k++ ) {
                   if( aaa_ratio[i + 1][j+1][k + 1] > 0.0 ) {
                      CHECK( ( idx_t )( vIdx[0] + i ) < Info::getDomainSize( 0 ) );
                      for (S32 dIndx =0 ; dIndx < NUM_DIFFUSIBLE_ELEMS; dIndx++ ) { 
                           a_avg_diff[ dIndx  ]  += v_gridPhiNbrBox[ dIndx  ].getVal( i, j, k ) * aaa_ratio[i + 1][j+1][k + 1];
                      }        
                   }
                }
             }
          }
       }
         //cout << "VIdx:"<<vIdx[0] << "," << vIdx<<[1]<<","<<vIDx[0]

         //cout<<"phi " <<  v_gridPhiNbrBox[0].getVal(0,0,0) << " phi -z:" <<v_gridPhiNbrBox[0].getVal(0,0,1)<< " phi +z:" <<v_gridPhiNbrBox[0].getVal(0,1,+1) << " phi -x-y:" << v_gridPhiNbrBox[0].getVal(-1,-1,0) <<" phi -x:" << v_gridPhiNbrBox[0].getVal(-1,0,0) << " phi -x+y:" << v_gridPhiNbrBox[0].getVal(-1,1,0)   <<endl;
         
       //REAL dt = BASELINE_TIME_STEP_DURATION/NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE;
       REAL Biomas = state.getModelReal(CELL_MODEL_REAL_BIOMAS );
       REAL Inert = state.getModelReal(CELL_MODEL_REAL_INERT ); 

       if ( A_BIOMASS_ODE_INDEX[ type ] != -1 ) {
          S32 ode_idx = A_BIOMASS_ODE_INDEX[ type ] ; 
          Biomas = state.getODEVal( 0 ,  ode_idx  );  
       }
       
       REAL cellVol = ( Biomas + Inert ) / A_DENSITY_BIOMASS[ type ] ;

       if ( cellVol > MAX_CELL_VOL ) {
          cellVol = MAX_CELL_VOL ;
          Biomas = cellVol * A_DENSITY_BIOMASS[ type ]  - Inert ;
       } 
       CHECK( cellVol >= 0.0 );
       REAL newRadius = radius_from_volume( cellVol );     
         //cout<<"biomas:"<<Biomas<<" rad:"<<newRadius<<  " uscale:"<<uScale<<" dt: " << dt <<    endl;
          
       state.setRadius( newRadius );
       state.setModelReal( CELL_MODEL_REAL_BIOMAS, Biomas );

       if( uptakePct < 1.0 ) {
           CHECK( UPTAKE_PCT_INC_RATIO >= 1.0 );
           uptakePct *= UPTAKE_PCT_INC_RATIO;
           if( uptakePct > 1.0 ) {
              uptakePct = 1.0;
           }
           state.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, uptakePct );
       }

       if( secretionPct < 1.0 ) {
          CHECK( SECRETION_PCT_CHANGE_RATIO >= 1.0 );
          secretionPct *= SECRETION_PCT_CHANGE_RATIO;
          if( secretionPct > 1.0 ) {
             secretionPct = 1.0;
          }
          state.setModelReal( CELL_MODEL_REAL_SECRETION_PCT, secretionPct );
       }

         // check if bnd should be gnerated 
 
       if ( A_AGENT_BOND_BOUNDARY_S[type] > 0.0 ) {
          REAL R0 =A_AGENT_SHOVING_SCALE[type]*state.getRadius();
          REAL x=((REAL)vIdx[0] +0.5+vOffset[0])*IF_GRID_SPACING;
          REAL dist_b = FABS(x - AGAR_HEIGHT*IF_GRID_SPACING);
          if (state.getModelInt(CELL_MODEL_INT_BOND_B)==0){
             if (dist_b < R0*A_AGENT_BOND_BOUNDARY_CREATE[type]){
                state.setModelInt(CELL_MODEL_INT_BOND_B,1); 
             } 
          } 
          else{
             if (dist_b> R0*A_AGENT_BOND_BOUNDARY_DESTROY[type]){
               state.setModelInt(CELL_MODEL_INT_BOND_B,0); 
             }   
          }
       }
    }
    else {/* dead */
           // Still nothing to do here
    }

    /* MODEL END */
    return;
}

void ModelRoutine::spAgentSecretionBySpAgent( const VIdx& vIdx, const AgentJunctionInfo& junctionInfo, const VReal& vOffset, const AgentMechIntrctData& mechIntrctData, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, SpAgentState& state/* INOUT */, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentDisp ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

void ModelRoutine::updateSpAgentBirthDeath( const VIdx& vIdx, const SpAgent& spAgent, const AgentMechIntrctData& mechIntrctData, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, BOOL& divide, BOOL& disappear ) {
    /* MODEL START */

    divide = false;
    disappear = false;

    //if ( (Info::getCurBaselineTimeStep() % 20)  == 0 ){ 

       REAL rnd_num = Util::getModelRand(MODEL_RNG_UNIFORM_10PERCENT); //0.9-1.1
       REAL testrad = DIVISION_RADIUS * rnd_num;
     

       if( ( spAgent.state.getModelReal( CELL_MODEL_REAL_UPTAKE_PCT ) > 0.0 )/* live cells */   ) { 
           if( spAgent.state.getRadius() >= testrad ) { 
              divide = true;
           }
           else if ( spAgent.state.getRadius() <= MIN_CELL_RADIUS ) {
              disappear = true;
           }
       }
    //}

    /* MODEL END */

    return;
}

void ModelRoutine::adjustSpAgent( const VIdx& vIdx, const AgentJunctionInfo& junctionInfo, const VReal& vOffset, const AgentMechIntrctData& mechIntrctData, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, SpAgentState& state/* INOUT */, VReal& disp ) {
	/* MODEL START */

    disp = mechIntrctData.force;

    //cout<<"dist: " <<SQRT(disp[0]*disp[0]+disp[1]*disp[1]+disp[2]*disp[2])<<endl;

    // Force due to Bond with agar
    if (state.getModelInt(CELL_MODEL_INT_BOND_B) == 1) {
        REAL x = ((REAL)vIdx[0] + 0.5 + vOffset[0] )*IF_GRID_SPACING;
        REAL dist_b = FABS(x -  AGAR_HEIGHT * IF_GRID_SPACING);
        S32 type = state.getType(); 
        REAL xij = A_AGENT_SHOVING_SCALE[type]*state.getRadius() - dist_b;
#if REAL_IS_FLOAT
        disp[0] +=  xij * tanhf( FABS(xij) * A_AGENT_BOND_BOUNDARY_S[type] );
#else
        disp[0] +=  xij * tanh( FABS(xij) * A_AGENT_BOND_BOUNDARY_S[type] );
#endif

    }   
    
 
    for( S32 dim = 0 ; dim < SYSTEM_DIMENSION; dim++ ) {/* limit the maximum displacement within a single time step */
        if( disp[dim] > MAX_CELL_RADIUS ) {
            disp[dim] = MAX_CELL_RADIUS ;
        }
	else if( disp[dim] < ( MAX_CELL_RADIUS * -1.0 ) ) {
             disp[dim] = MAX_CELL_RADIUS * -1.0 ;
        }
    }

    //cout << "Radius " << state.getRadius() << endl;
    //disp[2] = 0.0; // because cells move only in x-y (2D)
    /* MODEL END */
    return;
}

void ModelRoutine::divideSpAgent( const VIdx& vIdx, const AgentJunctionInfo& junctionInfo, const VReal& vOffset, const AgentMechIntrctData& mechIntrctData, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, SpAgentState& motherState/* INOUT */, VReal& motherDisp, SpAgentState& daughterState, VReal& daughterDisp, Vector<BOOL>& v_junctionDivide, BOOL& motherDaughterLinked, JunctionEnd& motherEnd, JunctionEnd& daughterEnd ) {
	/* MODEL START */

    CHECK( ( motherState.getType() == AGENT_TYPE_MyGrowingYeast ) && ( motherState.getModelReal( CELL_MODEL_REAL_UPTAKE_PCT ) > 0.0 )/* live */  );

    REAL radius;
    VReal dir;
    REAL scale;
    S32 type_id = motherState.getType(); 
    REAL cellVol;
    REAL biomas, mother_biomas, dougther_biomas ;
    REAL inert, mother_inert, dougther_inert;
    REAL rnd_num1 = Util::getModelRand(MODEL_RNG_UNIFORM_10PERCENT);//0.9-1.1


    motherDisp = mechIntrctData.force;

    biomas = motherState.getModelReal( CELL_MODEL_REAL_BIOMAS );
    inert = motherState.getModelReal( CELL_MODEL_REAL_INERT );
    mother_biomas = 0.5 * biomas * rnd_num1 ;
    mother_inert = 0.5 * inert * rnd_num1 ;

    cellVol = ( mother_biomas + mother_inert)  / A_DENSITY_BIOMASS[ type_id ]; 
    if ( cellVol < MIN_CELL_VOL){
       cellVol = MIN_CELL_VOL;
       mother_biomas =  cellVol * A_DENSITY_BIOMASS[ type_id ]  -  mother_inert  ; 
    }

    radius = radius_from_volume( cellVol ); 
    //  set the model variable 
    motherState.setRadius( radius );
    motherState.setModelReal( CELL_MODEL_REAL_BIOMAS, mother_biomas  );
    motherState.setModelReal( CELL_MODEL_REAL_INERT, mother_inert  );
    motherState.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, 1.0 );
    motherState.setModelReal( CELL_MODEL_REAL_SECRETION_PCT, 1.0 );

    dougther_biomas = biomas - mother_biomas;
    dougther_inert = inert - mother_inert;

    cellVol = ( dougther_biomas + dougther_inert) / A_DENSITY_BIOMASS[ type_id ];
    if ( cellVol < MIN_CELL_VOL){
       cellVol = MIN_CELL_VOL;
       dougther_biomas =  cellVol * A_DENSITY_BIOMASS[ type_id ];
    }

    radius = radius_from_volume( cellVol );    
    // set the model variable
    daughterState.setType( type_id );
    daughterState.setRadius( radius );
    daughterState.setModelReal( CELL_MODEL_REAL_BIOMAS, dougther_biomas );
    daughterState.setModelReal( CELL_MODEL_REAL_INERT, dougther_inert );
    daughterState.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, 1.0 );
    daughterState.setModelReal( CELL_MODEL_REAL_SECRETION_PCT, 1.0 );

    // change the new biomass on the  ODEs
    if ( A_BIOMASS_ODE_INDEX[ type_id ] != -1 ) {
       S32 ode_idx = A_BIOMASS_ODE_INDEX[type_id];
       motherState.setODEVal(0, ode_idx, mother_biomas);
       daughterState.setODEVal(0, ode_idx, dougther_biomas);
    }

    //dir[0] = 0.0;
    //dir[1] = 0.0;

    /* divide in a random direction */
    dir = VReal::ZERO;
    do {
        scale = 0.0;
        for( S32 dim = 0; dim < SYSTEM_DIMENSION; dim++ ) {
            dir[dim] = Util::getModelRand( MODEL_RNG_UNIFORM ) - 0.5;
            scale += dir[dim] * dir[dim];
        }
        scale = SQRT( scale );
    } while( scale < REAL_MIN );

    for( S32 dim = 0 ; dim < SYSTEM_DIMENSION ; dim++ ) {
        dir[dim] /= scale;
    }

    //////////////
    //dir[2]=0; // Needs checking for 2D or 3D simulation 
    //cout << "Displacement:" <<  dir[0] << "," << dir[1] << "," << dir[2] << endl;
  
    radius = 0.5* A_AGENT_SHOVING_SCALE[ type_id] *(daughterState.getRadius()+motherState.getRadius());
    motherDisp += dir * radius;
    daughterDisp -= dir * radius;


    for( S32 dim = 0 ; dim < SYSTEM_DIMENSION ; dim++ ) {/* limit the maximum displacement  */
        if( motherDisp[dim] > MAX_CELL_RADIUS ) 
            motherDisp[dim] = MAX_CELL_RADIUS ;
        else if( motherDisp[dim] < ( MAX_CELL_RADIUS * -1.0 ) ) 
            motherDisp[dim] = MAX_CELL_RADIUS * -1.0;
        
        if( daughterDisp[dim] > MAX_CELL_RADIUS  ) 
            daughterDisp[dim] = MAX_CELL_RADIUS ;
        else if( daughterDisp[dim] < ( MAX_CELL_RADIUS  * -1.0 ) ) 
            daughterDisp[dim] = MAX_CELL_RADIUS * -1.0;
    }

    //daughterDisp[2] = 0; // Needs checking for 2D or 3D simulations
    CHECK( junctionInfo.getNumJunctions() == 0 );
    v_junctionDivide.clear();

    motherDaughterLinked = false;

    /* MODEL END */
    return;
}
#endif

