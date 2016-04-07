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
                         S32 atype ;
                         //CHECK( MAX_CELL_RADIUS >= MIN_CELL_RADIUS );

                         posOffset[0] = cellInitData.x_offset ;
                         posOffset[1] = cellInitData.y_offset ;
                         posOffset[2] = cellInitData.z_offset ;
                         atype  = cellInitData.a_type;

                         state.setType( atype );
                         REAL biomass  = cellInitData.biomass ; 
                         REAL inert =  cellInitData.inert ;
                         REAL total_biomass = biomass + inert;  
                  
		         REAL volume = total_biomass / A_DENSITY_BIOMASS[atype]; 
                         REAL radius = radius_from_volume( volume ); 
   
                         if ( radius >  A_MAX_CELL_RADIUS[atype]  )
                             radius = A_MAX_CELL_RADIUS[atype];
                         else if ( radius < A_MIN_CELL_RADIUS[atype] )
                             radius = A_MIN_CELL_RADIUS[atype];
                       
                         state.setRadius( radius );
                         state.setModelReal( CELL_MODEL_REAL_BIOMAS, biomass );
                         state.setModelReal( CELL_MODEL_REAL_INERT, inert );
                         state.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, 1.0 );
                         state.setModelReal( CELL_MODEL_REAL_SECRETION_PCT, 1.0 );
                         state.setModelInt( CELL_MODEL_INT_BOND_B, 0 ) ;
   
                         if (  A_NUM_ODE_NET_VAR[atype] > 0 ) {
                            for( S32 i = 0 ; i < A_NUM_ODE_NET_VAR[atype]   ; i++ ) {
                               state.setODEVal(0,i,0.0);
                            }
                            if ( A_BIOMASS_ODE_INDEX[ atype ] != -1 ) {
                               S32 ode_idx = A_BIOMASS_ODE_INDEX[atype ] ;
                               state.setODEVal(0,ode_idx,biomass);
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


void ModelRoutine::updateSpAgentState( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */ ) {
    /* MODEL START */

    REAL aaa_ratio[3][3][3];

    Util::computeSphereUBVolOvlpRatio( SPHERE_UB_VOL_OVLP_RATIO_MAX_LEVEL, vOffset, state.getRadius(), aaa_ratio );

    REAL uptakePct = state.getModelReal( CELL_MODEL_REAL_UPTAKE_PCT );
    REAL secretionPct = state.getModelReal( CELL_MODEL_REAL_SECRETION_PCT );
    S32 type = state.getType() ; // id of cell type 

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
                           a_avg_diff[ dIndx  ]  += nbrUBEnv.getPhi(i,j,k,dIndx)*aaa_ratio[i+1][j+1][k+1];
                      }        
                   }
                }
             }
          }
          for ( S32 i = 1; i <= NUM_DIFFUSIBLE_ELEMS ; i++ )
             state.setModelReal( CELL_NUM_MODEL_REALS - i, a_avg_diff[NUM_DIFFUSIBLE_ELEMS-i]) ;
       }
       
     
       REAL Biomas = state.getModelReal(CELL_MODEL_REAL_BIOMAS);
       REAL Inert = state.getModelReal(CELL_MODEL_REAL_INERT); 

       if ( A_BIOMASS_ODE_INDEX[ type ] != -1 ) {
          S32 ode_idx = A_BIOMASS_ODE_INDEX[ type ] ; 
          Biomas = state.getODEVal( 0, ode_idx);  
       }
       
       REAL cellVol = (Biomas + Inert)/A_DENSITY_BIOMASS[type];

       if ( cellVol > A_MAX_CELL_VOL[type] ) {
          cellVol = A_MAX_CELL_VOL[type] ;
          Biomas = cellVol * A_DENSITY_BIOMASS[type]  - Inert ;
       } 
       CHECK( cellVol >= 0.0 );
       REAL newRadius = radius_from_volume( cellVol );     
          
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
          REAL x=((REAL)vIdx[0] +0.5)*IF_GRID_SPACING + vOffset[0];
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

void ModelRoutine::spAgentSecretionBySpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentDisp ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

void ModelRoutine::updateSpAgentBirthDeath( const VIdx& vIdx, const SpAgent& spAgent, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, BOOL& divide, BOOL& disappear ) {
    /* MODEL START */

    divide = false;
    disappear = false;

    S32 type = spAgent.state.getType() ;

    REAL rnd_num = Util::getModelRand(MODEL_RNG_UNIFORM_10PERCENT); //0.9-1.1
    REAL testrad = A_DIVISION_RADIUS[type] * rnd_num;
     
    if ((spAgent.state.getModelReal(CELL_MODEL_REAL_UPTAKE_PCT)>0.0)/* live cells */   ) { 
        if( spAgent.state.getRadius() >= testrad ) { 
            divide = true;
        }
        else if ( spAgent.state.getRadius() <= A_MIN_CELL_RADIUS[type] ) {
            disappear = true;
        }
    }
    
    // Remove cells when the touch the Borders
    if ( A_AGENT_BORDER_DISAPPEAR[0] ){
        if ((vIdx[0] == AGAR_HEIGHT) || (vIdx[0] == Info::getDomainSize(0)-1))
            disappear = true ;
    }
    if ( A_AGENT_BORDER_DISAPPEAR[1] ){
        if ((vIdx[1] == 0) || (vIdx[1] == Info::getDomainSize(1)-1))
            disappear = true ;
    }
    if ( A_AGENT_BORDER_DISAPPEAR[2] ){
        if ((vIdx[2] == 0) || (vIdx[2] == Info::getDomainSize(2)-1))
            disappear = true ;
    }

    /* MODEL END */

    return;
}

void ModelRoutine::adjustSpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */, VReal& disp ) {
    /* MODEL START */

    disp[0] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_X );
    disp[1] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_Y );
    disp[2] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_Z );
    S32 type = state.getType(); 

    // Force due to Bond with agar
    if (state.getModelInt(CELL_MODEL_INT_BOND_B) == 1) {
        REAL x=((REAL)vIdx[0]+0.5)*IF_GRID_SPACING + vOffset[0];
        REAL dist_b = FABS(x -  AGAR_HEIGHT * IF_GRID_SPACING);
        REAL xij = A_AGENT_SHOVING_SCALE[type]*state.getRadius() - dist_b;
#if REAL_IS_FLOAT
        disp[0] +=  xij * tanhf( FABS(xij) * A_AGENT_BOND_BOUNDARY_S[type] );
#else
        disp[0] +=  xij * tanh( FABS(xij) * A_AGENT_BOND_BOUNDARY_S[type] );
#endif

    }   
 
    // Random movement (Brownian)
    if ( A_DIFFUSION_COEFF_CELLS[type] > 0.0 ){
        REAL F_prw = SQRT( 2*A_DIFFUSION_COEFF_CELLS[type] * BASELINE_TIME_STEP_DURATION );
        for( S32 dim = 0 ; dim < SYSTEM_DIMENSION ; dim++ ) 
           disp[dim]+= F_prw* Util::getModelRand(MODEL_RNG_GAUSSIAN);
    }

    // limiting the displacement 
    for( S32 dim = 0 ; dim < SYSTEM_DIMENSION; dim++ ) {/* limit the maximum displacement within a single time step */
        if( disp[dim] > A_MAX_CELL_RADIUS[type] ) {
            disp[dim] = A_MAX_CELL_RADIUS[type] ;
        }
	else if( disp[dim] < ( A_MAX_CELL_RADIUS[type] * -1.0 ) ) {
             disp[dim] = A_MAX_CELL_RADIUS[type] * -1.0 ;
        }
    }

    for( S32 epr = 0; epr < NUM_E_PERTURBATIONS; epr++) {
      ExtConditions econd = A_E_PERTURBATIONS[epr];
      if ( Info::getCurBaselineTimeStep() == econd.TimeStep){
        if ( (econd.AgentType==-1) || ( econd.AgentType==type) ){
           REAL xx = ( REAL(vIdx[0]) + 0.5 )* IF_GRID_SPACING + vOffset[0] ;
           REAL yy = ( REAL(vIdx[1]) + 0.5 )* IF_GRID_SPACING + vOffset[1] ;
           REAL zz = ( REAL(vIdx[2]) + 0.5 )* IF_GRID_SPACING + vOffset[2] ;
           xx -= (REAL)(AGAR_HEIGHT*IF_GRID_SPACING);

           BOOL inside = false ;
           if ((xx > econd.xo )  && (xx < econd.xf)) {
             if ((yy > econd.yo) && (yy < econd.yf)) {
               if ( SYSTEM_DIMENSION == 2)
                   inside = true;
               else if ( (zz>econd.zo) && (zz<econd.zf) )
                   inside = true;
             }
           }
           if ( inside )  {
              if ( econd.Var_Index != -1 )
                 state.setModelReal( econd.Var_Index, econd.Var_Value ) ;
              if ( (econd.ODE_Index != -1) && ( econd.AgentType == -1 ) )
                 state.setODEVal(0, econd.ODE_Index, econd.ODE_Value ) ;
           }
        }
      }
    }

    /* MODEL END */
    return;
}

void ModelRoutine::divideSpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& motherState/* INOUT */, VReal& motherDisp, SpAgentState& daughterState, VReal& daughterDisp, Vector<BOOL>& v_junctionDivide, BOOL& motherDaughterLinked, JunctionEnd& motherEnd, JunctionEnd& daughterEnd ) {
	/* MODEL START */

    CHECK( ( motherState.getType() == AGENT_TYPE_MyGrowingYeast ) && ( motherState.getModelReal( CELL_MODEL_REAL_UPTAKE_PCT ) > 0.0 )/* live */  );

    REAL radius;
    VReal dir;
    REAL scale;
    S32 type_id = motherState.getType(); 
    REAL OldVol, MotherVol, DougtherVol ; 
    REAL biomas, mother_biomas, dougther_biomas ;
    REAL inert, mother_inert, dougther_inert;
    REAL rnd_num1 = Util::getModelRand(MODEL_RNG_UNIFORM_10PERCENT);//0.9-1.1

    motherDisp[0] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_X );
    motherDisp[1] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_Y );
    motherDisp[2] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_Z );

    biomas = motherState.getModelReal( CELL_MODEL_REAL_BIOMAS );
    inert = motherState.getModelReal( CELL_MODEL_REAL_INERT );
    OldVol = (biomas + inert)/A_DENSITY_BIOMASS[ type_id ];
 
    mother_biomas = 0.5 * biomas * rnd_num1 ;
    mother_inert = 0.5 * inert * rnd_num1 ;

    MotherVol = ( mother_biomas + mother_inert)  / A_DENSITY_BIOMASS[ type_id ]; 
    if ( MotherVol < A_MIN_CELL_VOL[type_id] ){
       MotherVol = A_MIN_CELL_VOL[type_id];
       mother_biomas = MotherVol*A_DENSITY_BIOMASS[type_id] - mother_inert; 
    }

    radius = radius_from_volume( MotherVol ); 
    //  set the model variable 
    motherState.setRadius( radius );
    motherState.setModelReal( CELL_MODEL_REAL_BIOMAS, mother_biomas  );
    motherState.setModelReal( CELL_MODEL_REAL_INERT, mother_inert  );
    motherState.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, 1.0 );
    motherState.setModelReal( CELL_MODEL_REAL_SECRETION_PCT, 1.0 );
    
    dougther_biomas = biomas - mother_biomas;
    dougther_inert = inert - mother_inert;

    DougtherVol = (dougther_biomas + dougther_inert)/A_DENSITY_BIOMASS[type_id];
    if ( DougtherVol < A_MIN_CELL_VOL[type_id] ){
       DougtherVol = A_MIN_CELL_VOL[type_id] ;
       dougther_biomas = DougtherVol*A_DENSITY_BIOMASS[type_id];
    }

    radius = radius_from_volume( DougtherVol );    
    // set the model variable
    daughterState.setType( type_id );
    daughterState.setRadius( radius );
    daughterState.setModelReal( CELL_MODEL_REAL_BIOMAS, dougther_biomas );
    daughterState.setModelReal( CELL_MODEL_REAL_INERT, dougther_inert );
    daughterState.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, 1.0 );
    daughterState.setModelReal( CELL_MODEL_REAL_SECRETION_PCT, 1.0 );
  
    // Copy values of ODEs from mother to dauther
    if ( A_NUM_ODE_NET_VAR[type_id] > 0 ) {
       for( S32 i = 0; i<A_NUM_ODE_NET_VAR[type_id]; i++ ) {
          REAL mol = motherState.getODEVal(0,i);          
          motherState.setODEVal(0,i, mol*OldVol/MotherVol)  ;
          daughterState.setODEVal(0,i, mol*OldVol/DougtherVol)  ;
       }
    }

    // change the new biomass on the  ODEs
    if ( A_BIOMASS_ODE_INDEX[ type_id ] != -1 ) {
       S32 ode_idx = A_BIOMASS_ODE_INDEX[type_id];
       motherState.setODEVal(0, ode_idx, mother_biomas);
       daughterState.setODEVal(0, ode_idx, dougther_biomas);
    }

    // divide in a random direction
    dir = VReal::ZERO;
    do {
        scale = 0.0;
        for( S32 dim = 0; dim < SYSTEM_DIMENSION; dim++ ) {
            dir[dim] = Util::getModelRand( MODEL_RNG_UNIFORM ) - 0.5;
            scale += dir[dim] * dir[dim];
        }
        scale = SQRT( scale );
    } while( scale > 0.5 );

    for( S32 dim = 0 ; dim < SYSTEM_DIMENSION ; dim++ ) {
        dir[dim] /= scale;
    }

    //////////////
  
    radius = 0.5* A_AGENT_SHOVING_SCALE[ type_id] *(daughterState.getRadius()+motherState.getRadius());
    motherDisp += dir * radius;
    daughterDisp -= dir * radius;


    for( S32 dim = 0 ; dim < SYSTEM_DIMENSION ; dim++ ) {/* limit the maximum displacement  */
        if( motherDisp[dim] > A_MAX_CELL_RADIUS[type_id] ) 
            motherDisp[dim] = A_MAX_CELL_RADIUS[type_id] ;
        else if( motherDisp[dim] < ( A_MAX_CELL_RADIUS[type_id] * -1.0 ) ) 
            motherDisp[dim] = A_MAX_CELL_RADIUS[type_id] * -1.0;
        
        if( daughterDisp[dim] > A_MAX_CELL_RADIUS[type_id]  ) 
            daughterDisp[dim] = A_MAX_CELL_RADIUS[type_id] ;
        else if( daughterDisp[dim] < ( A_MAX_CELL_RADIUS[type_id]  * -1.0 ) ) 
            daughterDisp[dim] = A_MAX_CELL_RADIUS[type_id] * -1.0;
    }

    //daughterDisp[2] = 0; // Needs checking for 2D or 3D simulations
    CHECK( junctionInfo.getNumJunctions() == 0 );
    v_junctionDivide.clear();

    motherDaughterLinked = false;

    /* MODEL END */
    return;
}
#endif

