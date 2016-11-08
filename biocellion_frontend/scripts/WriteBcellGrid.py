import os
import sys 
from BiocellionParam  import diffusible_solutes, cell_types, domain_parameters, mechanical_parameters, multigrid_solver_parm, basic_simulation_param 

def write_biocell_grid( diffusibles, celltypes, myreactions, myforces, mydomain, mygridsolver, mysimulator, directory ):

 # write info of the agents
 gridf = open(directory+"/model_routine_grid.cpp", 'w')
 input_gridf  = open('template/model_routine_grid.cpp', 'r')

 for line in input_gridf:
     gridf.write(line)

 gridf.write("\n") 
 gridf.write("void ModelRoutine::updateIfSubgridRHSLinear( const S32 elemIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridRHS/* uptake(-) and secretion (+) */ ) {\n\n")  
 
 gridf.write("gridRHS = 0.0 ;\n\n")
 for sol in diffusibles:
   #gridf.write("if ( elemIdx == DIFFUSIBLE_ELEM_"+sol+" ){\n")
    # Search if sol is consummed or secreted by agents
    SolFound={}
    Found = False; 
    for cell in celltypes:
       SolFound[cell] = False 
    for cell in celltypes:
       for rfactor in celltypes[cell]['reactions']:
          for mol in myreactions[rfactor]['yields']:
              if ( mol == sol ) :
                 SolFound[cell] = True  
                 Found= True
                 break
       if  SolFound[cell]  :
            break
    # check if sol  is affected inside cells
    if  not Found :
      continue
 
        
    gridf.write("if ( elemIdx == DIFFUSIBLE_ELEM_"+sol+" ){\n")
    gridf.write("   REAL rhs_up = 0.0; \n")
    gridf.write("   REAL rhs_down = 0.0; \n")
     

    gridf.write("//   for( S32 i = -1 ; i <= 1 ; i++){\n")
    gridf.write("//      for( S32 j = -1 ; j <= 1 ; j++){\n")
    gridf.write("//         for( S32 k = -1 ; k <= 1 ; k++){\n")
    gridf.write("//\t\tconst UBAgentData& ubAgentData = *( ubAgentDataPtrNbrBox.getVal(i,j,k));\n")
    gridf.write("//\t\tVIdx ubVIdxOffset;\n")
    gridf.write("//\t\tubVIdxOffset[0] = i * -1;\n")
    gridf.write("//\t\tubVIdxOffset[1] = j * -1;\n")
    gridf.write("//\t\tubVIdxOffset[2] = k * -1;\n")
    gridf.write("\t\tfor( ubAgentIdx_t l = 0 ; l < ( ubAgentIdx_t )ubAgentData.v_spAgent.size() ; l++ ) {   \n")
    gridf.write("\t\t\tconst SpAgent& spAgent = ubAgentData.v_spAgent[l];\n")
    gridf.write("\t\t\tagentType_t type=spAgent.state.getType();\n")
    gridf.write("//\t\t\tREAL ratio = Util::computeSphereUBVolOvlpRatio( SPHERE_UB_VOL_OVLP_RATIO_MAX_LEVEL, spAgent.vOffset, spAgent.state.getRadius(), ubVIdxOffset );\n")
    gridf.write("//\t\t\tif( ratio > 0.0 ) {\n")
    gridf.write("\t\t\t\t \n")
    gridf.write("\t\t\t\tREAL uptakePct;\n")
    gridf.write("\t\t\t\tREAL secretionPct;\n")
    
    for cell in celltypes:
       if not SolFound[cell] :
          continue
       gridf.write("\t\t\t\tif ( type == AGENT_TYPE_"+cell+"){\n")
       gridf.write("\t\t\t\t\tuptakePct = spAgent.state.getModelReal( CELL_MODEL_REAL_UPTAKE_PCT );\n")
       gridf.write("\t\t\t\t\tsecretionPct = spAgent.state.getModelReal( CELL_MODEL_REAL_SECRETION_PCT );\n")
       gridf.write("\n")
       for mysol in diffusibles: 
          gridf.write("\t\t\t\t\tREAL "+mysol+" = ubEnv.getPhi(DIFFUSIBLE_ELEM_"+mysol+");\n")
       gridf.write("\n")
       for mymol in celltypes[cell]['molecules']:
          gridf.write("\t\t\t\t\tREAL "+mymol+" = spAgent.state.getODEVal(0,ODE_NET_VAR_"+cell+"_"+mymol+");\n")
       gridf.write("\n")

       # Print each reaction factor
       for rfactor in celltypes[cell]['reactions']:

          gridf.write("\t\t\t\t\tREAL r_"+rfactor+"=")

          muMax = myreactions[rfactor]['muMax']
          if ( muMax == 0.0 ) :
             gridf.write( ' 0.0 ; \n ')
             continue
          else :
            gridf.write( str(muMax) )

          for MondEq in myreactions[rfactor]['MonodKinetic']:
             if ( not (MondEq['Ks'] == 0.0) ):
                gridf.write("*MonodEquation("+str(MondEq['Ks'])+","+ MondEq['solute']+")")
          for SimpleInh  in myreactions[rfactor]['SimpleInhibition']:
             gridf.write("*SimpleInhibition("+str(SimpleInh['Ki'])+","+ SimpleInh['solute']+")")
          for binding in myreactions[rfactor]['Binding']:
             gridf.write("*" + binding['solute'] )

          if ( myreactions[rfactor]['catalyzedby'] == "" ) :
             gridf.write(";\n")
          else:
             VolScale = ""
             if ( myreactions[rfactor]['FluxFlag'] ):
                #if myreactions[rfactor]['catalyzedby'] in celltypes[cell]['molecules'] :
                VolScale = "*surface_agent(spAgent.state.getRadius())/(IF_GRID_SPACING*IF_GRID_SPACING*IF_GRID_SPACING)"

             gridf.write( "*"+ myreactions[rfactor]['catalyzedby']+ VolScale+";\n")

       gridf.write("\n")
       # Print the yields 
       for rfactor in celltypes[cell]['reactions']:
          yIdx = -1
          for i in range(0,len(myreactions[rfactor]['yields'])):
             if ( myreactions[rfactor]['yields'][i] == sol ):
                 yIdx = i
       
          if ( yIdx == -1 ):
              continue
    
          yieldf = myreactions[rfactor]['yieldFactors'][yIdx]
          if ( yieldf < 0 ) :
              gridf.write("\t\t\t\t\trhs_down += "+ str(yieldf) + "*r_"+rfactor +"*uptakePct  ;\n" )
          else : 
              gridf.write("\t\t\t\t\trhs_up += "+ str(yieldf) + "*r_"+rfactor +"*secretionPct  ;\n" )
              
        
       gridf.write("\t\t\t\t}\n")
       gridf.write("//\t\t\t}\n")
       gridf.write("\t\t}\n")
 
    
    gridf.write("//         }\n")
    gridf.write("//      }\n")
    gridf.write("//   }\n")
    gridf.write("   gridRHS = rhs_up +  rhs_down;\n") 
    gridf.write("}\n\n")

 gridf.write("return;\n")
 gridf.write("}\n\n")
