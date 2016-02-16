import os
import sys 
from BiocellionParam  import diffusible_solutes, cell_types, domain_parameters, mechanical_parameters, multigrid_solver_parm, basic_simulation_param 

def write_biocell_config( diffusibles, celltypes, myreactions, myforces, mydomain, mygridsolver, mysimulator, directory ):

 # write info of the agents
 config = open(directory+"/model_routine_config.cpp", 'w')
 input_config  = open('template/model_routine_config.cpp', 'r')

 for line in input_config:
     config.write(line)


 config.write("#if HAS_SPAGENT\n")
 config.write("void ModelRoutine::updateSpAgentInfo( Vector<SpAgentInfo>& v_spAgentInfo ) {\n")
 config.write("\t///// UpdateSpAgentInfo \n")
 config.write("\tSpAgentInfo info; \n")
 config.write("\tv_spAgentInfo.resize(NUM_AGENT_TYPES);\n")

 dMax = 0.0;
 for cell in celltypes:
    idx = int (celltypes[cell]['id'])
    dMaxi = 2*celltypes[cell]['divRadius']*myforces[idx][idx]['shoveFactor']
    dMaxi = dMaxi * myforces[idx][idx]['attachDestroyFactor']
    if ( dMaxi > dMax ) :
       dMax = dMaxi
        

 for cell in celltypes:

    config.write("\tinfo.dMax = " + str( dMax  ) +";\n")
    config.write("\tinfo.hasBool = false;\n")
    config.write("\tinfo.numBoolVars = 0;\n")
    config.write("\tinfo.numStateModelReals = ")
    config.write(" NUM_MODEL_REALS; \n") # his will change cell->
    config.write("\tinfo.numStateModelInts = NUM_MODEL_INTS ;\n")
    config.write("\tinfo.numExtraMechIntrctModelReals = 0;\n")
    config.write("\tinfo.numExtraMechIntrctModelInts = 0;\n")
    config.write("\tinfo.v_odeNetInfo.clear();\n")
    config.write("\tv_spAgentInfo[AGENT_TYPE_")
    config.write(cell + "] = info;\n\n")
 config.write("\treturn;\n")
 config.write("}\n")
 config.write("#endif\n\n")


 config.write("void ModelRoutine::updatePDEInfo( Vector<PDEInfo>& v_pdeInfo ) {\n")
 config.write("\t///// UpdatePDEInfo \n")
 config.write("\tPDEInfo pdeInfo;\n")
 config.write("\tGridPhiInfo gridPhiInfo;\n")
 config.write("\tv_pdeInfo.resize( NUM_DIFFUSIBLE_ELEMS );\n")

 for sol in diffusibles:
    config.write("\tpdeInfo.pdeType = PDE_TYPE_REACTION_DIFFUSION_TIME_DEPENDENT_LINEAR;\n");
    config.write("\tpdeInfo.numLevels = A_NUM_AMR_LEVELS[DIFFUSIBLE_ELEM_")
    config.write(sol + "];\n")
    config.write("\tpdeInfo.numTimeSteps = A_NUM_PDE_TIME_STEPS_PER_STATE_AND_GRID_STEP[DIFFUSIBLE_ELEM_")
    config.write(sol + "];\n")
    config.write("\tpdeInfo.callAdjustRHSTimeDependentLinear=false;\n")
    config.write("\tpdeInfo.mgSolveInfo.numPre = 3;\n")
    config.write("\tpdeInfo.mgSolveInfo.numPost = 3;\n")
    config.write("\tpdeInfo.mgSolveInfo.numBottom = 3;\n")
    config.write("\tpdeInfo.mgSolveInfo.vCycle = true;\n")
    config.write("\tpdeInfo.mgSolveInfo.maxIters = 30;\n")
    config.write("\tpdeInfo.mgSolveInfo.epsilon = 1e-8;\n")
    config.write("\tpdeInfo.mgSolveInfo.hang = 1e-8;\n")
    config.write("\tpdeInfo.mgSolveInfo.normThreshold = A_MG_NORM_THRESHOLD[DIFFUSIBLE_ELEM_")
    config.write(sol + "];\n")
    config.write("\tgridPhiInfo.elemIdx=DIFFUSIBLE_ELEM_"+sol+";\n")
    config.write('\tgridPhiInfo.name = "'+sol+ '";' + "\n")
    config.write("\tgridPhiInfo.aa_bcType[0][0] =")
    config.write(mygridsolver['bcType_x-'] + ";\n")
    config.write("\tgridPhiInfo.aa_bcVal[0][0] = ") 
    config.write( str(mygridsolver['bcVal_x-']) + ";\n")

    config.write("\tgridPhiInfo.aa_bcType[0][0] =")
    config.write(mygridsolver['bcType_x-'] + ";\n")
    config.write("\tgridPhiInfo.aa_bcVal[0][0] = ") 
    config.write( str(mygridsolver['bcVal_x-']) + ";\n")


    config.write("\tgridPhiInfo.aa_bcType[0][1] =")
    config.write(mygridsolver['bcType_x+'] + ";\n")
    config.write("\tgridPhiInfo.aa_bcVal[0][1] = ") 
    config.write( str(mygridsolver['bcVal_x+']) + ";\n")

    config.write("\tgridPhiInfo.aa_bcType[1][0] =")
    config.write(mygridsolver['bcType_y-'] + ";\n")
    config.write("\tgridPhiInfo.aa_bcVal[1][0] = ") 
    config.write( str(mygridsolver['bcVal_y-']) + ";\n")

    config.write("\tgridPhiInfo.aa_bcType[1][1] =")
    config.write(mygridsolver['bcType_y+'] + ";\n")
    config.write("\tgridPhiInfo.aa_bcVal[1][1] = ") 
    config.write( str(mygridsolver['bcVal_y+']) + ";\n")

    config.write("\tgridPhiInfo.aa_bcType[2][0] =")
    config.write(mygridsolver['bcType_z-'] + ";\n")
    config.write("\tgridPhiInfo.aa_bcVal[2][0] = ") 
    config.write( str(mygridsolver['bcVal_z-']) + ";\n")

    config.write("\tgridPhiInfo.aa_bcType[2][1] =")
    config.write(mygridsolver['bcType_z+'] + ";\n")
    config.write("\tgridPhiInfo.aa_bcVal[2][1] = ") 
    config.write( str(mygridsolver['bcVal_z+']) + ";\n")
 
    config.write("\tgridPhiInfo.errorThresholdVal=A_MG_NORM_THRESHOLD[")
    config.write("DIFFUSIBLE_ELEM_"+sol+"];\n")
   
    config.write("\tgridPhiInfo.warningThresholdVal=A_MG_NORM_THRESHOLD[")
    config.write("DIFFUSIBLE_ELEM_"+sol+"];\n")
    
    config.write("\tgridPhiInfo.setNegToZero = true;\n")
    config.write("\tpdeInfo.v_gridPhiInfo.assign(1,gridPhiInfo );\n")
    config.write("\tv_pdeInfo[DIFFUSIBLE_ELEM_")
    config.write(sol + "] = pdeInfo;") 
    config.write("\n\n")
  
 config.write("\n\treturn;\n");
 config.write("}\n\n");
 
 

    
   
 config.write("void  ModelRoutine::updateIfGridModelVarInfo( Vector<IfGridModelVarInfo>& v_ifGridModelRealInfo, Vector<IfGridModelVarInfo>& v_ifGridModelIntInfo ) {\n\n")
 config.write("\tIfGridModelVarInfo info;\n")
 config.write("\tv_ifGridModelRealInfo.resize(NUM_GRID_MODEL_REALS);\n")
 config.write("\tinfo.name = \"colony volume ratio\";\n")
 config.write("\tv_ifGridModelRealInfo[GRID_MODEL_REAL_COLONY_VOL_RATIO] = info;\n")

 for sol in diffusibles:
 
    config.write("\tinfo.name = \"" + sol + "_u_scale\";\n")
    config.write("\tv_ifGridModelRealInfo[GRID_MODEL_REAL_"+sol+"_U_SCALE] = info;\n") 
    config.write("\tinfo.name = \""+sol+"_rhs\";\n")
    config.write("\tv_ifGridModelRealInfo[GRID_MODEL_REAL_"+sol+"_RHS] = info;\n")
    config.write("\tinfo.name = \"" + sol + "_rhs_from_high\";\n")
    config.write("\tv_ifGridModelRealInfo[GRID_MODEL_REAL_"+sol+"_RHS_FROM_HIGH] = info;\n")

 config.write("\tv_ifGridModelIntInfo.clear();\n\n")
 config.write("\treturn ;\n")
 config.write("}\n\n")



 
 config.write("void ModelRoutine::updateSummaryOutputInfo( SummaryOutputInfo& summaryOutputInfo ) {\n\n" )
 config.write("\tsummaryOutputInfo.v_realName.resize( NUM_GRID_SUMMARY_REALS );\n")
 config.write("\tsummaryOutputInfo.v_realType.resize( NUM_GRID_SUMMARY_REALS );\n")

 for sol in diffusibles:
    config.write("\tsummaryOutputInfo.v_realName[GRID_SUMMARY_REAL_")
    config.write( sol+"] = \"" + sol + "\";\n")
    config.write("\tsummaryOutputInfo.v_realType[GRID_SUMMARY_REAL_")
    config.write( sol+"] = SUMMARY_TYPE_MAX;\n")

 for cell in celltypes:
    config.write("\tsummaryOutputInfo.v_realName[GRID_SUMMARY_REAL_LIVE_")
    config.write( cell+"] = \"" +cell+ " cells(live) volume\"; \n")
    config.write("\tsummaryOutputInfo.v_realType[GRID_SUMMARY_REAL_LIVE_")
    config.write( cell+"] = SUMMARY_TYPE_SUM;\n ")
   
 config.write("\n\treturn;\n")
 config.write("}\n\n")

   

    


