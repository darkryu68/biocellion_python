/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

/* UESR START */

#include "model_define.h"

/* UESR END */

using namespace std;

#if HAS_SPAGENT
void ModelRoutine::updateSpAgentOutput( const VIdx& vIdx, const SpAgent& spAgent, REAL& color, Vector<REAL>& v_extra ) {

   /* MODEL START */
   if(spAgent.state.getModelReal(CELL_MODEL_REAL_UPTAKE_PCT)>0.0 ){/*live */
      color = spAgent.state.getType();
   }
   else {/* dead */
      color = NUM_AGENT_TYPES ;
   }
   CHECK( v_extra.size() == 0 );

   v_extra[0] = REAL( spAgent.junctionInfo.getCurId() )  ;
   v_extra[1] = spAgent.state.getModelReal( CELL_MODEL_REAL_BIOMAS );
   v_extra[2] = spAgent.state.getModelReal( CELL_MODEL_REAL_INERT );

   // display values of ODEs  
   for ( S32 i = 3 ; i < NUM_AGENT_OUTPUTS; i++ ) {
       S32 type = spAgent.state.getType(); 
       if ( AA_INDEX_ODE_OUTPUT[ type ][i-3] < 0 )
           v_extra[i] = -1.0 ;  
       else 
          v_extra[i] = spAgent.state.getODEVal(0, AA_INDEX_ODE_OUTPUT[ type ][i-3] ) ;
   }

   /* MODEL END */
   return;
}
#endif

void ModelRoutine::updateSummaryVar( const VIdx& vIdx, const NbrBox<const UBAgentData*>& ubAgentDataPtrNbrBox, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, Vector<REAL>& v_realVal/* [elemIdx] */, Vector<S32>& v_intVal/* [elemIdx] */ ) {
   /* MODEL START */

   const UBAgentData& ubAgentData = *( ubAgentDataPtrNbrBox.getVal( 0, 0, 0 ) );
   REAL rCellVol[ NUM_AGENT_TYPES ] ;
   for ( S32 type = 0 ; type < NUM_AGENT_TYPES; type++) {
      rCellVol[ type ] = 0.0 ;
   }

   CHECK( v_realVal.size() == NUM_GRID_SUMMARY_REALS );
   CHECK( v_intVal.size() == 0 );

   // This is surely not good practice 
   for (S32 pdeIdx = 0; pdeIdx < NUM_DIFFUSIBLE_ELEMS; pdeIdx++){
       v_realVal[pdeIdx] = v_gridPhiNbrBox[pdeIdx].getVal(0,0,0);
   }

   for( S32 i = 0 ; i < ( S32 )ubAgentData.v_spAgent.size() ; i++ ) {
       const SpAgent& spAgent = ubAgentData.v_spAgent[i];
       S32 type = spAgent.state.getType() ;
 
       if(spAgent.state.getModelReal(CELL_MODEL_REAL_UPTAKE_PCT)>0.0) {
          rCellVol[type] += 1.0 ;
       }	
   }
   
   // assuming order in summary structure, diffusibles first, agents later
   for ( S32 type = 0 ; type < NUM_AGENT_TYPES; type++) {
      v_realVal[type + NUM_DIFFUSIBLE_ELEMS] = rCellVol[type] ;
   }
   /* MODEL END */
   return;
}

