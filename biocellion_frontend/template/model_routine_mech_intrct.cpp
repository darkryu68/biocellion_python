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

void ModelRoutine::computeForceSpAgent( const VIdx& vIdx0, const SpAgent& spAgent0, const Vector<REAL>& v_gridPhi0/* [elemIdx] */, const Vector<REAL>& v_gridModelReal0/* [elemIdx] */, const Vector<S32>& v_gridModelInt0/* [elemIdx] */, const VIdx& vIdx1, const SpAgent& spAgent1, const Vector<REAL>& v_gridPhi1/* [elemIdx] */, const Vector<REAL>& v_gridModelReal1/* [elemIdx] */, const Vector<S32>& v_gridModelInt1/* [elemIdx] */, const VReal& dir/* unit direction vector from spAgent1 to spAgent0 */, const REAL& dist, VReal& force/* force on spAgent0 due to interaction with spAgent1 (force on spAgent1 due to interaction with spAgent0 has same magnitude but the opposite direction), if force has the same direction with dir, two cells push each other, if has the opposite direction, two cells pull each other. */ ) {
	/* MODEL START */
        S32 type0 = spAgent0.state.getType();  
        S32 type1 = spAgent1.state.getType(); 

	REAL R0 = A_AGENT_SHOVING_SCALE[type0]*spAgent0.state.getRadius() ;
        REAL R1 = A_AGENT_SHOVING_SCALE[type1]*spAgent1.state.getRadius();
        REAL D = R0 + R1 ;
	REAL mag = 0.0;/* + for repulsive force, - for adhesive force */
        REAL xij  = D - dist;

        //if (  spAgent0.junctionInfo.getCurId() == 48  )  {
        //   cout<<dist<<" "<<spAgent1.junctionInfo.getCurId()<<" "<<R0<<" "<<R1<<endl;
        //} 
      

	if( dist <= D  ) {/* shoving to remove the overlap */
            mag = 0.5 * ( xij );
	}
        else {/* adhesion */
            if ( A_AGENT_ADHESION_S[type0][type1] == 0.0 ) {
                mag = 0.0;
            }
            else {
                //REAL x = dist / D;
#if REAL_IS_FLOAT
                mag = 0.5 * xij * expf(- xij*xij  / A_AGENT_ADHESION_S[type0][type1]  );
#else
                mag = 0.5 * xij * exp( - xij*xij  / A_AGENT_ADHESION_S[type0][type1]  );
#endif
                //mag = 0.0; 
            }
        }

        // Bonds
        if (spAgent0.junctionInfo.isLinked(spAgent1.junctionInfo) == true) {
            REAL sij = A_AGENT_BOND_S[type0][type1] ;
            if ( sij  >  0.0 ) {    
#if REAL_IS_FLOAT
                mag = mag + 0.5 * xij * tanhf( FABS(xij) * sij);
#else
                mag = mag + 0.5 * xij * tanh( FABS(xij) * sij);
#endif
            }
        } 

 
        force  =  VReal::ZERO;
        for( S32 dim = 0 ; dim < SYSTEM_DIMENSION ; dim++ ) {
                force[dim] = mag * dir[dim];
        }

        //force[2] = 0.0 ;
	/* MODEL END */

	return;
}

void ModelRoutine::computeExtraMechIntrctSpAgent( const VIdx& vIdx0, const SpAgent& spAgent0, const Vector<REAL>& v_gridPhi0/* [elemIdx] */, const Vector<REAL>& v_gridModelReal0/* [elemIdx] */, const Vector<S32>& v_gridModelInt0/* [elemIdx] */, const VIdx& vIdx1, const SpAgent& spAgent1, const Vector<REAL>& v_gridPhi1/* [elemIdx] */, const Vector<REAL>& v_gridModelReal1/* [elemIdx] */, const Vector<S32>& v_gridModelInt1/* [elemIdx] */, const VReal& dir/* unit direction vector from spAgent1 to spAgent0 */, const REAL& dist, ExtraMechIntrctData& extraMechIntrctData0, ExtraMechIntrctData& extraMechIntrctData1, BOOL& link, JunctionEnd& end0/* dummy if link == false */, JunctionEnd& end1/* dummy if link == false */, BOOL& unlink ) {
	/* MODEL START */

	link = false;
	unlink = false;

        S32 type0 = spAgent0.state.getType();
        S32 type1 = spAgent1.state.getType();

        REAL R0 = A_AGENT_SHOVING_SCALE[type0]*spAgent0.state.getRadius();
        REAL R1 = A_AGENT_SHOVING_SCALE[type1]*spAgent1.state.getRadius();

        REAL dist_threshold = A_AGENT_BOND_DESTROY_FACTOR[type0] *(R0 + R1);

        if( spAgent0.junctionInfo.isLinked( spAgent1.junctionInfo ) == true ) {
           if( dist > dist_threshold ) {
              unlink = true;/* break junction */
           }
        }
        else {/* no junction */
           if(  dist < dist_threshold  ) {
              link = true;/* form junction */
              end0.setType(0);
              end1.setType(0);
           }
        }

	/* MODEL END */

	return;
}
#endif

