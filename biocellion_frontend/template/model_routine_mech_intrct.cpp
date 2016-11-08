/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

/* MODEL START */

#include "model_define.h"

/* MODEL END */

using namespace std;

#if HAS_SPAGENT  
void ModelRoutine::initJunctionSpAgent( const VIdx& vIdx0, const SpAgent& spAgent0, const VIdx& vIdx1, const SpAgent& spAgent1, const VReal& dir/* unit direction vector from spAgent1 to spAgent0 */, const REAL& dist, BOOL& link, JunctionEnd& end0/* dummy if link == false */, JunctionEnd& end1/* dummy if link == false */ ) {
	/* MODEL START */

        S32 type0 = spAgent0.state.getType();
        S32 type1 = spAgent1.state.getType();

        REAL R0 = A_AGENT_SHOVING_SCALE[type0]*spAgent0.state.getRadius();
        REAL R1 = A_AGENT_SHOVING_SCALE[type1]*spAgent1.state.getRadius();

        REAL dist_threshold = A_AGENT_BOND_CREATE_FACTOR[type0] *(R0 + R1);
        if (  dist < dist_threshold )  {
          link = true;
          end0.setType(0);
          end1.setType(0);
        }
        else {
          link = false;
        }

	/* MODEL END */

	return;
}

void ModelRoutine::computeMechIntrctSpAgent( const S32 iter, const VIdx& vIdx0, const SpAgent& spAgent0, const UBEnv& ubEnv0, const VIdx& vIdx1, const SpAgent& spAgent1, const UBEnv& ubEnv1, const VReal& dir/* unit direction vector from spAgent1 to spAgent0 */, const REAL& dist, MechIntrctData& mechIntrctData0, MechIntrctData& mechIntrctData1, BOOL& link, JunctionEnd& end0/* dummy if link == false */, JunctionEnd& end1/* dummy if link == false */, BOOL& unlink ) {

    /* MODEL START */
    link = false;
    unlink = false;

    S32 type0 = spAgent0.state.getType();
    S32 type1 = spAgent1.state.getType();

    REAL R0 = A_AGENT_SHOVING_SCALE[type0]*spAgent0.state.getRadius();
    REAL R1 = A_AGENT_SHOVING_SCALE[type1]*spAgent1.state.getRadius();
    REAL dist_threshold = R0 + R1;
    REAL D = R0 + R1 - 0.5*A_AGENT_SHOVING_LIMIT[type0] - 0.5*A_AGENT_SHOVING_LIMIT[type1]; 
    REAL mag = 0.0;
    REAL xij  = D - dist ; 
    
  
    if( dist <= D ) {/* shoving */
        mag = 0.5 * (xij) ;
    }
    else {/* adhesion */
        if( A_AGENT_ADHESION_S[type0][type1] > 0.0 )  {
#if REAL_IS_FLOAT
            mag = 0.5 * xij*expf( -xij*xij / A_AGENT_ADHESION_S[type0][type1] );
#else
            mag = 0.5 * xij*exp( -xij*xij / A_AGENT_ADHESION_S[type0][type1] );
#endif 
        }
    }         

    if ( A_AGENT_BOND_S[type0][type1] > 0.0 ){
        REAL sij = A_AGENT_BOND_S[type0][type1] ;
        if(spAgent0.junctionData.isLinked(spAgent1.junctionData) == true) {
            if( dist > A_AGENT_BOND_DESTROY_FACTOR[type0]* dist_threshold ) {
                unlink = true;/* break junction */
            }
            else{
                // compute elastic force
                REAL D = R0 + R1;
                REAL xij  = D - dist  ;
#if REAL_IS_FLOAT
                REAL Fij = 0.5 * xij * tanhf(FABS(xij)*sij);
#else
                REAL Fij = 0.5 * xij * tanh(FABS(xij)*sij);
#endif
                mag = mag + Fij ;
            }
        }
        else {/* no junction */
            if( dist < A_AGENT_BOND_CREATE_FACTOR[type0]*dist_threshold ) {
                link = true;/* form junction */
                end0.setType(0);
                end1.setType(0);

                // add force rigth away
                REAL D = R0 + R1;
                REAL xij  = D - dist  ;
#if REAL_IS_FLOAT
                REAL Fij = 0.5 * xij * tanhf(FABS(xij)*sij);
#else
                REAL Fij = 0.5 * xij * tanh(FABS(xij)*sij);
#endif
                mag = mag + Fij ;
            }
        }
    }

    mechIntrctData0.setModelReal(CELL_MECH_REAL_FORCE_X,dir[0]*mag);
    mechIntrctData0.setModelReal(CELL_MECH_REAL_FORCE_Y,dir[1]*mag);
    mechIntrctData0.setModelReal(CELL_MECH_REAL_FORCE_Z,dir[2]*mag);

    mechIntrctData1.setModelReal(CELL_MECH_REAL_FORCE_X,-dir[0]*mag);
    mechIntrctData1.setModelReal(CELL_MECH_REAL_FORCE_Y,-dir[1]*mag);
    mechIntrctData1.setModelReal(CELL_MECH_REAL_FORCE_Z,-dir[2]*mag);

    /* MODEL END */

    return;
}
#endif

