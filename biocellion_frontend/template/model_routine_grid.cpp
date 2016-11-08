/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

#include "model_define.h"

using namespace std;

void ModelRoutine::initIfGridVar( const VIdx& vIdx, const UBAgentData& ubAgentData, UBEnv& ubEnv ) {
    /* MODEL START */
    for ( S32 sol = 0 ; sol < NUM_DIFFUSIBLE_ELEMS ; sol++ )  {  

       if(  vIdx[0] >= AGAR_HEIGHT ) {
           ubEnv.setPhi( sol,  A_INIT_IF_CONCENTRATION[sol]);
       }
       else {
           ubEnv.setPhi( sol, A_INIT_AGAR_CONCENTRATION[sol]) ;
       }

       if ( (SYSTEM_DIMENSION == 2) && (vIdx[2]>0) )
           ubEnv.setPhi(sol,0.0);
    }

    for( S32 sol = 0 ; sol < NUM_GRID_MODEL_REALS ; sol++ ) {
         ubEnv.setModelReal( sol, 0.0 ); 
    }
    /* MODEL END */
    return;
}

void ModelRoutine::updateIfGridVar( const BOOL pre, const S32 iter, const VIdx& vIdx, const NbrUBAgentData& nbrUBAgentData, NbrUBEnv& nbrUBEnv/*[INOUT]*/) {
    /* MODEL START */
    CHECK( pre == true );
    CHECK( iter == 0 );

    if ( ( Info::getCurBaselineTimeStep() < 1 ) && ( Info::getSimInitType() == SIM_INIT_TYPE_CODE )/* not starting from checkpoint data */ ) {
    /* wait for one baseline time step to reduce cell overlaps */
    }
    else {
        if( vIdx[0]  >= AGAR_HEIGHT  ) {/* the cell growth region (+ a top fraction of the agarose cylinder), we are ignoring cell penetration to the agarose cylinder */
            REAL dt;
            REAL colonyVolRatio = 0.0;
            dt = BASELINE_TIME_STEP_DURATION/NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE;

	    /* iterate over 3 * 3 * 3 boxes */
	    for( S32 i = -1 ; i <= 1 ; i++ ) {
	     for( S32 j = -1 ; j <= 1 ; j++ ) {
	      for( S32 k = -1 ; k <= 1 ; k++ ) {

                  if (((idx_t)(vIdx[0]+i)>=Info::getDomainSize(0)) || ((idx_t)(vIdx[0]+i) < 0)) 
                      continue;

                  if (((idx_t)(vIdx[1]+j)>=Info::getDomainSize(1)) || ((idx_t)(vIdx[1]+j) < 0)) 
                      continue;

                  if (((idx_t)(vIdx[2]+k)>=Info::getDomainSize(2)) || ((idx_t)(vIdx[2]+k) < 0)) 
                      continue;

                                 
                  const UBAgentData& ubAgentData = *( nbrUBAgentData.getConstPtr(i,j,k));
                  VIdx ubVIdxOffset;
                  ubVIdxOffset[0] = i* -1;
                  ubVIdxOffset[1] = j* -1;
                  ubVIdxOffset[2] = k* -1;

                  for(ubAgentIdx_t l=0; l<(ubAgentIdx_t)ubAgentData.v_spAgent.size(); l++ ) {
                      const SpAgent& spAgent = ubAgentData.v_spAgent[l];
                      REAL ratio = Util::computeSphereUBVolOvlpRatio( SPHERE_UB_VOL_OVLP_RATIO_MAX_LEVEL, spAgent.vOffset, spAgent.state.getRadius(), ubVIdxOffset );
                      
                      if( ratio > 0.0 ) {
                          REAL radius;
			  REAL vol;

                          /* compute colony volume */
                          radius = spAgent.state.getRadius();
                          vol = ( 4.0 * MY_PI / 3.0 ) * radius * radius * radius;
                          colonyVolRatio += vol * ratio;
                      }
		  }
              }
             }
	    }
        
            colonyVolRatio /= ( IF_GRID_SPACING * IF_GRID_SPACING * IF_GRID_SPACING ) * UB_FULL_COLONY_VOL_RATIO;
	    if( colonyVolRatio > 1.0 ) {
	        colonyVolRatio = 1.0;
	    }

            nbrUBEnv.setModelReal(0,0,0,GRID_MODEL_REAL_COLONY_VOL_RATIO,colonyVolRatio);
        }
        else {/* the agarose cylinder (excluding the top fraction) */
	    CHECK( vIdx[0]   <  AGAR_HEIGHT  );
	    CHECK( nbrUBEnv.getModelReal(0,0,0,GRID_MODEL_REAL_COLONY_VOL_RATIO ) <= 0.0 );
        }
    }

/* MODEL END */
return;

}

void ModelRoutine::updateIfSubgridKappa( const S32 pdeIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridKappa/* to consider cell volume exclusion in computing diffusion flux, set to 1.0 to ignore the volume occupied by cells */ ) {
	/* MODEL START */

        gridKappa = 1.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridAlpha( const S32 elemIdx, const VIdx& vIdx, const VIdx& subgridVOffset,  const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridAlpha/* decay (-) */ ) {
	/* MODEL START */

        // Decay Rate can be specified here
	gridAlpha = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridBetaInIfRegion( const S32 elemIdx, const S32 dim, const VIdx& vIdx0, const VIdx& subgridVOffset0, const UBAgentData& ubAgentData0, const UBEnv& ubEnv0, const VIdx& vIdx1, const VIdx& subgridVOffset1, const UBAgentData& ubAgentData1, const UBEnv& ubEnv1, REAL& gridBeta ) {
	/* MODEL START */

    REAL gridBeta0;
    REAL gridBeta1;

    if( vIdx0[0] < AGAR_HEIGHT ) {
        gridBeta0 = A_DIFFUSION_COEFF_AGAR[elemIdx];
    }
    else {
        gridBeta0 = A_DIFFUSION_COEFF_COLONY[elemIdx];
        if( A_DIFFUSION_COLONY_ONLY[elemIdx] ) {
            REAL scale = ubEnv0.getModelReal(GRID_MODEL_REAL_COLONY_VOL_RATIO);
            CHECK( scale <= 1.0 );
            gridBeta0 *= scale; 
        }
    }

    if( vIdx1[0] < AGAR_HEIGHT ) {
        gridBeta1 = A_DIFFUSION_COEFF_AGAR[elemIdx];
    }
    else {
        gridBeta1 = A_DIFFUSION_COEFF_COLONY[elemIdx];       
        if( A_DIFFUSION_COLONY_ONLY[elemIdx] ) {
            REAL scale = ubEnv1.getModelReal(GRID_MODEL_REAL_COLONY_VOL_RATIO);
            CHECK( scale <= 1.0 );
            gridBeta1 *= scale;
        }
    }

    if ( SYSTEM_DIMENSION == 2 ) {
        if ( vIdx0[2] > 0 )
            gridBeta0 = 0.0 ;
        if ( vIdx1[2] > 0 )
            gridBeta1 = 0.0 ;
    }

    if( ( gridBeta0 > 0.0 ) && ( gridBeta1 > 0.0 ) ) {
        gridBeta = 1.0 / ( ( 1.0 / gridBeta0 + 1.0 / gridBeta1 ) * 0.5 );/* harmonic mean */
    }
    else {
        gridBeta = 0.0;
    }

    /* MODEL END */
    return;
}

void ModelRoutine::updateIfSubgridBetaPDEBufferBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridBeta ) {
	/* MODEL START */

	CHECK( vIdx[0] > 0 );
	CHECK( vIdx[0]  < AGAR_HEIGHT );
	gridBeta = A_DIFFUSION_COEFF_AGAR[elemIdx];

        

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridBetaDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridBeta ) {
	/* MODEL START */

//	CHECK( dim == 2 );
//	CHECK( vIdx[2] == ( idx_t )( Info::getDomainSize( dim ) - 1 ) );

	gridBeta = 0.0;

	/* MODEL END */

	return;
}


//void ModelRoutine::adjustIfSubgridRHSTimeDependentLinear( const S32 elemIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData,  const UBEnv& ubEnv, REAL& gridRHS/* INOUT, uptake(-) and secretion (+) */ ) {
	/* MODEL START */

//	ERROR( "unimplemented." );

	/* MODEL END */

//	return;
//}

void ModelRoutine::updateIfSubgridRHSTimeDependentSplitting( const S32 pdeIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnvModelVar& ubEnvModelVar, const Vector<double>& v_gridPhi/* [idx] */, Vector<double>& v_gridRHS/* [idx], uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAMRTags( const VIdx& vIdx, const NbrUBAgentData& nbrUBAgentData, const NbrUBEnv& nbrUBEnv, Vector<S32>& v_finestLevel/* [pdeIdx] */ ) {
    /* MODEL START */

    v_finestLevel.assign(NUM_DIFFUSIBLE_ELEMS, 0 );/* coarsest level */
    //for( S32 pdeIdx = 0 ; pdeIdx < NUM_DIFFUSIBLE_ELEMS ; pdeIdx++ ) {
    //        v_finestLevel[pdeIdx] = A_NUM_AMR_LEVELS[pdeIdx] - 1;
    // }

    /* MODEL END */
    return;
}

void ModelRoutine::updateIfGridDirichletBCVal( const S32 elemIdx, const VReal& pos, const S32 dim, const BOOL lowSide, const UBEnvModelVar a_ubEnvModelVar[3], const Vector<REAL> a_gridPhi[3], REAL& bcVal ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridNeumannBCVal( const S32 elemIdx, const VReal& pos, const S32 dim, const BOOL lowSide, const UBEnvModelVar a_ubEnvModelVar[3], const Vector<REAL>  a_gridPhi[3], REAL& bcVal ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::initPDEBufferPhi( const S32 pdeIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, Vector<REAL>& v_gridPhi ) {
	/* MODEL START */

#if ENABLE_CHECK
	CHECK( startVIdx[0]  < AGAR_HEIGHT );
#endif

    //for( S32 pdeIdx = 0 ; pdeIdx < NUM_DIFFUSIBLE_ELEMS ; pdeIdx++ ) {
	v_gridPhi[0] = A_INIT_AGAR_CONCENTRATION[ pdeIdx ];
    //}


	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferKappa( const S32 pdeIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridKappa/* to consider cell volume exclusion in computing diffusion flux, set to 1.0 to ignore the volume occupied by cells */ ) {
	/* MODEL START */

	gridKappa = 1.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferAlpha( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridAlpha/* decay (-) */ ) {
	/* MODEL START */
        // Decay here
	gridAlpha = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferBetaInPDEBufferRegion( const S32 elemIdx, const S32 dim, const VIdx& startVIdx0, const VIdx& startVIdx1, const VIdx& pdeBufferBoxSize, REAL& gridBeta ) {
	/* MODEL START */

#if ENABLE_CHECK
	CHECK( startVIdx0[0] + pdeBufferBoxSize[0]  < AGAR_HEIGHT );
#endif

	gridBeta = A_DIFFUSION_COEFF_AGAR[elemIdx];

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferBetaDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridBeta ) {
	/* MODEL START */

	//CHECK( dim == 2 );
	//CHECK( startVIdx[2] == 0 );

	gridBeta = 0.0;

	/* MODEL END */

	return;
}

//void ModelRoutine::updatePDEBufferAdvectionVelocityInPDEBufferRegion( const S32 elemIdx, const S32 dim, const VIdx& startVIdx0, const VIdx& startVIdx1, const VIdx& pdeBufferBoxSize, REAL& gridAdvectionVelocity ) {
//	/* MODEL START */
//
//	ERROR( "unimplemented." );
//
//	/* MODEL END */
//
//	return;
//}

//void ModelRoutine::updatePDEBufferAdvectionVelocityDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridAdvectionVelocity ) {
//	/* MODEL START */
//
//	ERROR( "unimplemented." );
//
//	/* MODEL END */
//
//	return;
//}

void ModelRoutine::updatePDEBufferRHSLinear( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, const REAL gridPhi, REAL& gridRHS/* uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	gridRHS = 0.0;

	/* MODEL END */

	return;
}

//void ModelRoutine::adjustPDEBufferRHSTimeDependentLinear( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, const REAL gridPhi, REAL& gridRHS/* INOUT, uptake(-) and secretion (+) */ ) {
	/* MODEL START */

//	ERROR( "unimplemented." );

	/* MODEL END */

//	return;
//}

void ModelRoutine::updatePDEBufferRHSTimeDependentSplitting( const S32 pdeIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, const Vector<double>& v_gridPhi/* [idx] */, Vector<double>& v_gridRHS/* [idx], uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferDirichletBCVal( const S32 elemIdx, const VReal& startPos, const VReal& pdeBufferBoxSize, const S32 dim, const BOOL lowSide, REAL& bcVal ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferNeumannBCVal( const S32 elemIdx, const VReal& startPos, const VReal& pdeBufferBoxSize, const S32 dim, const BOOL lowSide, REAL& bcVal ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

