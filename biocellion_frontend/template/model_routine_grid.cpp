/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

#include "model_define.h"

using namespace std;

void ModelRoutine::initIfGridVar( const VIdx& vIdx, const UBAgentData& ubAgentData, Vector<REAL>& v_gridPhi/* [elemIdx] */, Vector<REAL>& v_gridModelReal/* [elemIdx] */, Vector<S32>& v_gridModelInt ) {
    /* MODEL START */


    for ( S32 sol = 0 ; sol < NUM_DIFFUSIBLE_ELEMS ; sol++ )  {  

       if(  vIdx[0] >= AGAR_HEIGHT ) {
          v_gridPhi[ sol ] = A_INIT_IF_CONCENTRATION[ sol ]; ;
       }
       else {
          v_gridPhi[ sol ] = A_INIT_AGAR_CONCENTRATION[ sol ];
       }

       if ( (SYSTEM_DIMENSION == 2) && (vIdx[2]>0) )
           v_gridPhi[ sol ] = 0.0;
    }


    for( S32 elemIdx = 0 ; elemIdx < NUM_GRID_MODEL_REALS ; elemIdx++ ) {
        v_gridModelReal[elemIdx] = 0.0;
    }
    /* MODEL END */
    return;
}

void ModelRoutine::updateIfGridVar( const BOOL pre, const S32 round, const VIdx& vIdx, const NbrBox<const UBAgentData*>& ubAgentDataPtrNbrBox, Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* INOUT, [elemIdx] */, Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* INOUT, [elemIdx] */ ) {
	/* MODEL START */

CHECK( pre == true );
CHECK( round == 0 );


if ( ( Info::getCurBaselineTimeStep() < 1 ) && ( Info::getSimInitType() == SIM_INIT_TYPE_CODE )/* not starting from checkpoint data */ ) {
    /* do nothing - wait for one baseline time step to reduce cell overlaps */
}
else {
    if( vIdx[0]  >= AGAR_HEIGHT  ) {/* the cell growth region (+ a top fraction of the agarose cylinder), we are ignoring cell penetration to the agarose cylinder */

        REAL dt;

        REAL colonyVolRatio = 0.0;
        //REAL a_uPrime[NUM_DIFFUSIBLE_ELEMS];/* rhs = uPrime * uScale * -1.0 + q */
        //REAL a_q[NUM_DIFFUSIBLE_ELEMS];
        //REAL a_uScale[NUM_DIFFUSIBLE_ELEMS];
        //REAL a_rhsToLow[NUM_DIFFUSIBLE_ELEMS];
        //REAL rhs;

        dt = BASELINE_TIME_STEP_DURATION / NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE;

        //for( S32 elemIdx = 0 ; elemIdx < NUM_DIFFUSIBLE_ELEMS ; elemIdx++ ) {
        //    a_uPrime[elemIdx] = 0.0;
	//    a_q[elemIdx] = 0.0;
	//    a_uScale[elemIdx] = 0.0;
	//    a_rhsToLow[elemIdx] = 0.0;
	//}

	/* iterate over 3 * 3 * 3 boxes */
	for( S32 i = -1 ; i <= 1 ; i++ ) {
	   for( S32 j = -1 ; j <= 1 ; j++ ) {
	      for( S32 k = -1 ; k <= 1 ; k++ ) {

                 if ( ((idx_t)(vIdx[0]+i)>=Info::getDomainSize(0)) || ((idx_t)(vIdx[0]+i) < 0) ) 
                     continue;

                 if ( ((idx_t)(vIdx[1]+j)>=Info::getDomainSize(1)) || ((idx_t)(vIdx[1]+j) < 0) ) 
                     continue;

                 if ( ((idx_t)(vIdx[2]+k)>=Info::getDomainSize(2)) || ((idx_t)(vIdx[2]+k) < 0) ) 
                     continue;

              
                 if( (idx_t)(vIdx[0] + i) < Info::getDomainSize(0)) {
                    
                    const UBAgentData& ubAgentData = *(ubAgentDataPtrNbrBox.getVal(i,j,k));
                    VIdx ubVIdxOffset;
                    ubVIdxOffset[0] = i* -1;
                    ubVIdxOffset[1] = j* -1;
                    ubVIdxOffset[2] = k* -1;


                    for(ubAgentIdx_t l=0; l<(ubAgentIdx_t)ubAgentData.v_spAgent.size(); l++ ) {
                        const SpAgent& spAgent = ubAgentData.v_spAgent[l];
                        //agentType_t type = spAgent.state.getType();
                        REAL ratio = Util::computeSphereUBVolOvlpRatio( SPHERE_UB_VOL_OVLP_RATIO_MAX_LEVEL, spAgent.vOffset, spAgent.state.getRadius(), ubVIdxOffset );
                        //S32 sourceIdx = A_AGENT_GROWTH_SOURCE[type];
                    

			if( ratio > 0.0 ) {
			    REAL radius;
			    REAL vol;
			    //REAL uptakePct;
			    //REAL secretionPct;
                            //REAL biomass;

			    /* compute colony volume */
                            radius = spAgent.state.getRadius();
                            vol = ( 4.0 * MY_PI / 3.0 ) * radius * radius * radius;
                            colonyVolRatio += vol * ratio;
                            //cout <<  "vol " << vol << " ratio " << ratio <<endl;
			    /* update uPrime and q */
                            //uptakePct = spAgent.state.getModelReal( CELL_MODEL_REAL_UPTAKE_PCT );
                            //secretionPct = spAgent.state.getModelReal( CELL_MODEL_REAL_SECRETION_PCT );
                            //biomass = spAgent.state.getModelReal( CELL_MODEL_REAL_BIOMAS) ; 
                            // here should go double for loop 
                            //a_uPrime[sourceIdx] += uptakePct*ratio * biomass   ;
                            //cout << "uprime:" <<  a_uPrime[0] <<" uptakePct:"<<uptakePct <<" ratio:"<<ratio<<endl;
			}
		    }
		 }
              }
           }
	}
        

        colonyVolRatio /= ( IF_GRID_SPACING * IF_GRID_SPACING * IF_GRID_SPACING ) * UB_FULL_COLONY_VOL_RATIO;
	if( colonyVolRatio > 1.0 ) {
	    colonyVolRatio = 1.0;
	}

	v_gridModelRealNbrBox[GRID_MODEL_REAL_COLONY_VOL_RATIO].setVal( 0, 0, 0, colonyVolRatio );

	//for( S32 elemIdx = 0 ; elemIdx < NUM_DIFFUSIBLE_ELEMS ; elemIdx++ ) {
	//    a_uPrime[elemIdx] /= IF_GRID_SPACING * IF_GRID_SPACING * IF_GRID_SPACING;
	//    a_q[elemIdx] /= IF_GRID_SPACING * IF_GRID_SPACING * IF_GRID_SPACING;
	//}

	/* compute uScale */
	//for( S32 elemIdx = 0 ; elemIdx < NUM_DIFFUSIBLE_ELEMS ; elemIdx++ ) {
	//    if( a_uPrime[elemIdx] > 0.0 ) {
	//       REAL oldUScale = 0.0;
	//	REAL effPhi;
	//	REAL newUScale;

	//	oldUScale = v_gridModelRealNbrBox[3*elemIdx + 1 ].getVal(0,0,0);

	//	effPhi = v_gridPhiNbrBox[elemIdx].getVal( 0, 0, 0 );
	//	if( vIdx[0]  >= AGAR_HEIGHT ) {
	//	    effPhi /= colonyVolRatio;
	//	}
		//CHECK( effPhi + A_Ks[elemIdx] > 0.0 );
		//newUScale = effPhi / ( effPhi + A_Ks[elemIdx] );
          //      newUScale = BiomasRate(0, 1.0, effPhi, 1.0 ) ;

	//	if( oldUScale > U_SCALE_LIMIT_THRESHOLD ) {
	//	    if( newUScale > oldUScale * U_SCALE_MAX_INC_RATIO ) {
	//	         newUScale = oldUScale * U_SCALE_MAX_INC_RATIO;
	//	    }
	//	}
	//	else {
	//	    if( newUScale > U_SCALE_LIMIT_THRESHOLD * U_SCALE_MAX_INC_RATIO ) {
	//	        newUScale = U_SCALE_LIMIT_THRESHOLD * U_SCALE_MAX_INC_RATIO;
	//	    }
	//	}

	//	a_uScale[elemIdx] = newUScale;
	//    }

	    //if( elemIdx == DIFFUSIBLE_ELEM_GLUCOSE ) {
	//    v_gridModelRealNbrBox[3*elemIdx+1].setVal(0,0,0,a_uScale[elemIdx]);
	    //}
	//}

        //if ( colonyVolRatio > 0.0 ) 
        //    cout <<" u':"<<a_uPrime[0]<<" u_s:"<<a_uScale[0]<<" q:"<<a_q[0]<< endl;

	// update rhs 
	//for( S32 elemIdx = 0 ; elemIdx < NUM_DIFFUSIBLE_ELEMS ; elemIdx++ ) {
	//  rhs = 20.0 *a_uPrime[elemIdx]*a_uScale[elemIdx] * -1.0 + a_q[elemIdx];
	//rhs += glucoseQ;
	//   v_gridModelRealNbrBox[3*elemIdx + 2 ].setVal(0,0,0,rhs);
        //}

	// nutrients floating in the air 
	//if ( ( vIdx[0] >= AGAR_HEIGHT ) && ( colonyVolRatio <= 0.0 ) ) {
	//    for( S32 elemIdx = 0 ; elemIdx < NUM_DIFFUSIBLE_ELEMS ; elemIdx++ ) {
	//        REAL phi = v_gridPhiNbrBox[elemIdx].getVal( 0, 0, 0 );
	//	if( phi > 0.0 ) {
	//	    if( phi > A_FLOATING_NUTRIENTS_MAX_DELTA[elemIdx] ) {
	//	        v_gridPhiNbrBox[elemIdx].setVal( 0, 0, 0, phi - A_FLOATING_NUTRIENTS_MAX_DELTA[elemIdx] );
	//		a_rhsToLow[elemIdx] = A_FLOATING_NUTRIENTS_MAX_DELTA[elemIdx] / dt;
	//	    }
	//	    else {
	//	        v_gridPhiNbrBox[elemIdx].setVal( 0, 0, 0, 0.0 );
	//	        a_rhsToLow[elemIdx] = phi / dt;
	//	    }
	//	}
	//   }
	//}

        //if ( vIdx[0] > AGAR_HEIGHT) { 
	//   for( S32 pdeIdx = 0 ; pdeIdx < NUM_DIFFUSIBLE_ELEMS ; pdeIdx++ ) {
        //      v_gridModelRealNbrBox[3*pdeIdx+3].setVal(-1,0,0,a_rhsToLow[pdeIdx]);
        //   }
        //}

        //if ((vIdx[0]==10) && (vIdx[1]==16) && (vIdx[2]==10000)) {
        //    cout << " VolRatio " << colonyVolRatio << endl;
        //    cout << " uprime " << a_uPrime[1] << endl;
        //    cout << " uscale " << a_uScale[1] << endl;
        //    cout << " q  " << a_q[1] << endl;
        //    cout << " rhs " <<   v_gridModelRealNbrBox[GRID_MODEL_REAL_Glucose_RHS].getVal(0,0,0) << endl;
        //    cout << " phi " <<  v_gridPhiNbrBox[1 ].getVal( 0, 0, 0 )<< endl ;
            
        //}
    }
    else {/* the agarose cylinder (excluding the top fraction) */
	CHECK( vIdx[0]   <  AGAR_HEIGHT  );
	CHECK( v_gridModelRealNbrBox[GRID_MODEL_REAL_COLONY_VOL_RATIO].getVal( 0, 0, 0 ) <= 0.0 );
    }
}

/* MODEL END */
return;

}

void ModelRoutine::updateIfGridKappa( const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi/* [elemIdx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridKappa/* to consider cell volume exclusion in computing diffusion flux, set to 1.0 to ignore the volume occupied by cells */ ) {
	/* MODEL START */

	//REAL z = ( ( REAL )vIdx[0] + 0.5 ) * IF_GRID_SPACING;
        gridKappa = 1.0;

	//if( vIdx[0] < AGAR_HEIGHT ) {
	//	gridKappa = 1.0;
	//}
	//else {
	//	REAL colonyVolRatio = v_gridModelReal[GRID_MODEL_REAL_COLONY_VOL_RATIO];
	//	if( colonyVolRatio <= 0.0 ) {
	//		gridKappa = 1.0;
	//	}
	//	else {
        //                gridKappa = 1.0;
			//gridKappa = colonyVolRatio;
	//	}
	//}

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAlpha( const S32 elemIdx, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi/* [elemIdx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridAlpha/* decay (-) */ ) {
	/* MODEL START */

	gridAlpha = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridBetaInIfRegion( const S32 elemIdx, const S32 dim, const VIdx& vIdx0, const VIdx& vIdx1, const UBAgentData& ubAgentData0, const UBAgentData& ubAgentData1, const Vector<REAL>& v_gridPhi0/* [elemIdx] */, const Vector<REAL>& v_gridPhi1/* [elemIdx] */, const Vector<REAL>& v_gridModelReal0/* [elemIdx] */, const Vector<REAL>& v_gridModelReal1/* [elemIdx] */, const Vector<S32>& v_gridModelInt0/* [elemIdx] */, const Vector<S32>& v_gridModelInt1/* [elemIdx] */, REAL& gridBeta ) {
	/* MODEL START */

    REAL gridBeta0;
    REAL gridBeta1;

    if( vIdx0[0] < AGAR_HEIGHT ) {
        gridBeta0 = A_DIFFUSION_COEFF_AGAR[elemIdx];
    }
    else {
        gridBeta0 = A_DIFFUSION_COEFF_COLONY[elemIdx];
        if( A_DIFFUSION_COLONY_ONLY[elemIdx] ) {
            REAL scale = v_gridModelReal0[GRID_MODEL_REAL_COLONY_VOL_RATIO];
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
            REAL scale = v_gridModelReal1[GRID_MODEL_REAL_COLONY_VOL_RATIO];
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
    
    //if ( v_gridModelReal0[GRID_MODEL_REAL_COLONY_VOL_RATIO] > 0.0 ) 
    //    cout << " BETAA " << gridBeta << endl ;

    /* MODEL END */
    return;
}

void ModelRoutine::updateIfGridBetaPDEBufferBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi/* [elemIdx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridBeta ) {
	/* MODEL START */

	CHECK( vIdx[0] > 0 );
	CHECK( vIdx[0]  < AGAR_HEIGHT );
	gridBeta = A_DIFFUSION_COEFF_AGAR[elemIdx];

        

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridBetaDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi/* [elemIdx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridBeta ) {
	/* MODEL START */

//	CHECK( dim == 2 );
//	CHECK( vIdx[2] == ( idx_t )( Info::getDomainSize( dim ) - 1 ) );

	gridBeta = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAdvectionVelocityInIfRegion( const S32 elemIdx, const S32 dim, const VIdx& vIdx0, const VIdx& vIdx1, const UBAgentData& ubAgentData0, const UBAgentData& ubAgentData1, const Vector<REAL>& v_gridPhi0/* [elemIdx] */, const Vector<REAL>& v_gridPhi1/* [elemIdx] */, const Vector<REAL>& v_gridModelReal0/* [elemIdx] */, const Vector<REAL>& v_gridModelReal1/* [elemIdx] */, const Vector<S32>& v_gridModelInt0/* [elemIdx] */, const Vector<S32>& v_gridModelInt1/* [elemIdx] */, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAdvectionVelocityPDEBufferBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi/* [elemIdx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAdvectionVelocityDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi/* [elemIdx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

//void ModelRoutine::updateIfGridRHSLinear( const S32 elemIdx, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi/* [elemIdx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridRHS/* uptake(-) and secretion (+) */ ) {
	/* MODEL START */

        //gridRHS=rhs+rhs_from_high;//check grid_model_real_e in model_define.h
//	gridRHS = v_gridModelReal[3*elemIdx+2] ;
       
        //if ( v_gridModelReal[GRID_MODEL_REAL_COLONY_VOL_RATIO]  > 0.0 )
        //          cout <<"rhs" <<  gridRHS << endl;        

          

	/* MODEL END */

//	return;
//}

void ModelRoutine::adjustIfGridRHSTimeDependentLinear( const S32 elemIdx, const VIdx& vIdx, const REAL gridPhi, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridRHS/* INOUT, uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridRHSTimeDependentSplitting( const S32 pdeIdx, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<double>& v_gridPhi/* [idx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, Vector<double>& v_gridRHS/* [idx], uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAMRTags( const VIdx& vIdx, const NbrBox<const UBAgentData*>& ubAgentDataPtrNbrBox, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, Vector<S32>& v_finestLevel/* [pdeIdx] */ ) {
    /* MODEL START */

    for( S32 pdeIdx = 0 ; pdeIdx < NUM_DIFFUSIBLE_ELEMS ; pdeIdx++ ) {
        //if( vIdx[0] < AGAR_HEIGHT ) {
        //    v_finestLevel[pdeIdx] = 0;
        //}
        //else {
            v_finestLevel[pdeIdx] = A_NUM_AMR_LEVELS[pdeIdx] - 1;
        //}
    }

    /* MODEL END */
    return;
}

void ModelRoutine::updateIfGridDirichletBCVal( const S32 elemIdx, const VReal& pos, const S32 dim, const BOOL lowSide, const REAL a_gridPhi[3], const Vector<Vector<REAL> >& vv_gridModelReal/* [elemIdx] */, const Vector<Vector<S32> >& vv_gridModelInt/* [elemIdx] */, REAL& bcVal ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridNeumannBCVal( const S32 elemIdx, const VReal& pos, const S32 dim, const BOOL lowSide, const REAL a_gridPhi[3], const Vector<Vector<REAL> >& vv_gridModelReal/* [elemIdx] */, const Vector<Vector<S32> >& vv_gridModelInt/* [elemIdx] */, REAL& bcVal ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::initPDEBufferPhi( const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, Vector<REAL>& v_gridPhi ) {
	/* MODEL START */

#if ENABLE_CHECK
	CHECK( startVIdx[0]  < AGAR_HEIGHT );
#endif

    for( S32 pdeIdx = 0 ; pdeIdx < NUM_DIFFUSIBLE_ELEMS ; pdeIdx++ ) {
	v_gridPhi[ pdeIdx] = A_INIT_AGAR_CONCENTRATION[pdeIdx];
    }


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

void ModelRoutine::updatePDEBufferAdvectionVelocityInPDEBufferRegion( const S32 elemIdx, const S32 dim, const VIdx& startVIdx0, const VIdx& startVIdx1, const VIdx& pdeBufferBoxSize, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferAdvectionVelocityDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferRHSLinear( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, const REAL gridPhi, REAL& gridRHS/* uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	gridRHS = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::adjustPDEBufferRHSTimeDependentLinear( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, const REAL gridPhi, REAL& gridRHS/* INOUT, uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

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

