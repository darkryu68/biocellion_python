import os
import sys
import math 
from BiocellionParam  import diffusible_solutes, cell_types, domain_parameters, mechanical_parameters, multigrid_solver_parm, basic_simulation_param 
from BiocellionParam import create_reaction, CreateMonodKinetic, CreateSimpleInhibition, create_e_perturbations

def write_biocell_header( diffusibles, celltypes, myreactions, myforces, eperturbations, mydomain, mygridsolver, mysimulator, directory ):

 num_celltypes  = len( celltypes )

 #mechIntrctData0.setModelReal(CELL_MECH_REAL_FORCE_X,dir[0]*mag);
 #mechIntrctData0.setModelReal(CELL_MECH_REAL_FORCE_Y,dir[1]*mag);
 #mechIntrctData0.setModelReal(CELL_MECH_REAL_FORCE_Z,dir[2]*mag);

 # write the files for biocellion
 header = open(directory+"/model_define.h", 'w')
 header.write('#ifndef __MY_DEFINE_H__' +"\n" )
 header.write('#define __MY_DEFINE_H__' +"\n\n")
 header.write('#include "biocellion.h"'+"\n\n")
 header.write('#define NUM_JUNCTION_END_TYPES 1 /* we consider only one junction end type */ '+"\n")
 header.write('#define SYSTEM_DIMENSION ')
 header.write( str(mydomain['nDim']) + "\n\n")
# header.write("const REAL MY_PI = 3.14159265358979323846;\n\n")

 # write functions used by biocellion model
 if ( mydomain['nDim']  == 2  ):
    header.write("inline REAL radius_from_volume( REAL volume ) {\n")
    header.write("\treturn SQRT( volume/( MY_PI*" + str(mydomain['resolution_agent']) +" )); \n")
    header.write("}\n")
    header.write("inline REAL volume_agent(REAL radius){\n")
    header.write("\treturn "+str(mydomain['resolution_agent'])+"*MY_PI*radius*radius; \n")
    header.write("}\n")
    header.write("inline REAL surface_agent(REAL radius){\n")
    header.write("\treturn "+str(mydomain['resolution_agent'])+"*2*MY_PI*radius; \n")
    header.write("}\n")
 else :
    header.write("inline REAL radius_from_volume( REAL volume ) {\n")
    header.write("\treturn CBRT( volume * 3.0 / ( 4.0 * MY_PI ) );\n")
    header.write("}\n\n")
    header.write("inline REAL volume_agent(REAL radius){\n")
    header.write("\treturn (4.0/3.0)*MY_PI*radius*radius*radius; \n")
    header.write("}\n")
    header.write("inline REAL surface_agent(REAL radius){\n")
    header.write("\treturn 4.0*MY_PI*radius*radius; \n")
    header.write("}\n")
 header.write("inline REAL MonodEquation(  REAL Kc , REAL u ) {\n")
 header.write("\treturn  u / ( Kc + u ) ;\n")
 header.write("}\n\n")

 header.write("typedef enum _agent_type_e {   \n")
 for cell in celltypes:
    header.write("\tAGENT_TYPE_" + cell + ",\n")
 header.write("\tNUM_AGENT_TYPES \n"  )
 header.write("} agent_type_e; \n\n"  )

 header.write("typedef enum _alive_cell_model_real_e {   \n")
 header.write("\tCELL_MODEL_REAL_BIOMAS,\n")
 header.write("\tCELL_MODEL_REAL_INERT,\n")
 header.write("\tCELL_MODEL_REAL_UPTAKE_PCT,\n")
 header.write("\tCELL_MODEL_REAL_SECRETION_PCT,\n")
 for sol in diffusibles:
     header.write("\tCELL_MODEL_REAL_"+sol+"_AVG,\n") 

 header.write("\tCELL_NUM_MODEL_REALS\n")
 header.write("} alive_cell_model_real_e;\n\n")

 for cell in celltypes:
     header.write("typedef enum _ode_net_"+ cell +"_var_e {   \n")
     for mol in celltypes[cell]['molecules']: 
        header.write("\tODE_NET_VAR_"+cell+"_"+mol+",\n")
     header.write("\tNUM_ODE_NET_VAR_"+ cell +"\n")
     header.write("} ode_net_"+cell+"_var_e;\n\n")

 header.write("typedef enum _alive_cell_model_int_e { \n")
 header.write("\tCELL_MODEL_INT_BOND_B, \n")
 header.write("\tNUM_MODEL_INTS \n")
 header.write("} alive_cell_model_int_e;\n\n")


 ## 
 header.write("typedef enum _cell_mech_real_e { \n" ) 
 header.write("\tCELL_MECH_REAL_FORCE_X,\n" ) 
 header.write("\tCELL_MECH_REAL_FORCE_Y,\n" ) 
 header.write("\tCELL_MECH_REAL_FORCE_Z,\n" ) 
 header.write("\tNUM_CELL_MECH_REALS\n" ) 
 header.write("} cell_mech_real_e;\n\n" ) 
 ##

 header.write("typedef enum _diffusible_elem_e {\n")
 for solute in diffusibles:
    header.write("\tDIFFUSIBLE_ELEM_"+solute+",\n")
 header.write("\tNUM_DIFFUSIBLE_ELEMS\n")
 header.write("} diffusible_elem_e;\n\n")

 header.write("typedef enum _grid_model_real_e {\n")
 header.write("\tGRID_MODEL_REAL_COLONY_VOL_RATIO,\n")
 for sol in diffusibles:
   header.write("\tGRID_MODEL_REAL_"+sol+"_U_SCALE,\n")
   header.write("\tGRID_MODEL_REAL_"+sol+"_RHS,\n")
   header.write("\tGRID_MODEL_REAL_"+sol+"_RHS_FROM_HIGH,\n")
 header.write("\tNUM_GRID_MODEL_REALS\n")
 header.write("} grid_model_real_e;\n\n")

 header.write("typedef enum _grid_summary_real_e {\n")
 for sol in diffusibles:
    header.write("\tGRID_SUMMARY_REAL_"+sol+",\n")
 for cell in celltypes:
    header.write("\tGRID_SUMMARY_REAL_LIVE_"+cell+",\n")
 header.write("\tNUM_GRID_SUMMARY_REALS\n")
 header.write("} grid_summary_real_e;\n\n")

 header.write("typedef enum _model_rng_type_e {\n")
 header.write("\tMODEL_RNG_UNIFORM,\n")
 header.write("\tMODEL_RNG_UNIFORM_10PERCENT,\n")
 header.write("\tMODEL_RNG_GAUSSIAN,\n")
 header.write("\tNUM_MODEL_RNGS\n")
 header.write("} model_rng_type_e;\n\n")

 header.write("struct IniCellData {\n")
 header.write("\tREAL x_offset;\n")
 header.write("\tREAL y_offset;\n")
 header.write("\tREAL z_offset;\n")
 header.write("\tS32 a_type;\n")
 header.write("\tREAL biomass;\n")
 header.write("\tREAL inert;\n")
 header.write("\tS32  IdxIniCellData;\n")
 header.write("};\n\n")

 header.write("struct ExtConditions {\n")
 header.write("\tREAL xo;\n")
 header.write("\tREAL xf;\n")
 header.write("\tREAL yo;\n")
 header.write("\tREAL yf;\n")
 header.write("\tREAL zo;\n")
 header.write("\tREAL zf;\n")
 header.write("\tS32 TimeStep;\n")
 header.write("\tS32 AgentType;\n")
 header.write("\tS32 Var_Index;\n")
 header.write("\tREAL Var_Value;\n")
 header.write("\tS32 ODE_Index;\n")
 header.write("\tREAL ODE_Value;\n")
 header.write("};\n\n")

 header.write("class UBInitData {\n")
 header.write("public:\n")
 header.write("\tS32 numCells ;  \n")
 header.write("\tS32  IdxIniCellData;\n")
 header.write("};\n\n")

 header.write("extern S32  INI_N_CELLS;\n\n")


 header.write("const REAL IF_GRID_SPACING = ")
 header.write( str(mydomain['resolution']) + ";\n\n")

 header.write("const S32 NUM_AGENT_OUTPUTS = ")
 MolOutput_array = [];
 for cell in celltypes :
    for mol in celltypes[cell]['moloutput']: 
       if mol in MolOutput_array : continue
       else:  MolOutput_array.append( mol  ) 
 header.write( str(3+ len(MolOutput_array)) + ";\n\n")
 

 header.write("const S32 A_NUM_AMR_LEVELS[NUM_DIFFUSIBLE_ELEMS]={") 
 comma = "" 
 for sol in diffusibles:
    header.write(comma + "2" ) 
    comma = ", "
 header.write("};\n\n")   

 header.write("const BOOL A_AGENT_BORDER_DISAPPEAR[3]={")
 if ( mygridsolver['dbtype_x'] == 'DISAPPEAR' ) :  header.write("true" )
 else: header.write("false" )
 if ( mygridsolver['dbtype_y'] == 'DISAPPEAR' ) :  header.write(",true" )
 else: header.write(",false" )
 if ( mygridsolver['dbtype_z'] == 'DISAPPEAR' ) :  header.write(",true" )
 else: header.write(",false" )
 header.write("};\n\n")


 header.write("const S32 AGAR_HEIGHT = ") 
 header.write( str(mydomain['agar_heigth'])+";\n" ) 
 header.write("const S32 AGAR_TOP_BUFFER_HEIGHT = 0;\n\n");

 header.write("const REAL A_INIT_AGAR_CONCENTRATION[NUM_DIFFUSIBLE_ELEMS]={")
 comma = "" 
 for sol in diffusibles:
     header.write(comma +str(diffusibles[sol]['ini_con_agar']))
     comma = ", "
     #header.write(","+str(sol['ini_con_agar'])   )
 header.write("};\n")

 header.write("const REAL A_INIT_IF_CONCENTRATION[NUM_DIFFUSIBLE_ELEMS]={")
 comma = ""
 for sol in diffusibles:
     header.write(comma + '0.0')
     comma = ", "
     #header.write(","+str(sol['ini_con_agar'])   )
 header.write("};\n")

 header.write("const REAL A_DIFFUSION_COEFF_AGAR[NUM_DIFFUSIBLE_ELEMS]={")
 comma = "" ; # looks horrrible
 for sol in diffusibles:
    diff =diffusibles[sol]['Dcoef_agar'] 
    header.write(comma +str(diff))
    comma = ", "
 header.write("};\n")

 header.write("const REAL A_DIFFUSION_COEFF_COLONY[NUM_DIFFUSIBLE_ELEMS]={")
 cont = 0 ; # looks horrrible
 for sol in diffusibles:
    diff =diffusibles[sol]['Dcoef_cd'] * mydomain['biofilmDiffusivity']
    if ( cont == 0 ) :
       header.write(str(diff))
       cont = 1
    else :
       header.write(", "+ str(diff))
 header.write("};\n")

 header.write("const BOOL A_DIFFUSION_COLONY_ONLY[NUM_DIFFUSIBLE_ELEMS]={")
 cont = 0 ; # looks horrrible
 for sol in diffusibles:
    flag = 'false'
    if ( diffusibles[sol]['colonyonly'] ):
        flag= "true"
    if ( cont == 0 ) :
       header.write( flag  )
       cont = 1
    else :
       header.write(", "+ flag  )
 header.write("};\n\n")

 header.write("const REAL A_NUM_ODE_NET_VAR[NUM_AGENT_TYPES]={")
 comma = ""
 for cell in celltypes:
     header.write( comma+ "NUM_ODE_NET_VAR_" +  cell )
     comma = ", "
 header.write("};\n")
 
 header.write("const REAL A_DENSITY_BIOMASS[NUM_AGENT_TYPES]={")
 comma = "" 
 for cell in celltypes:
     header.write( comma +   str( celltypes[cell]['biomass_density']))
     comma = ", "
 header.write("};\n\n")

 header.write("const S32 A_BIOMASS_ODE_INDEX[NUM_AGENT_TYPES]={")
 comma = "" 
 for cell in celltypes: 
     ode_idx = "-1"
     for reaction in celltypes[cell]['molecules']:
         if (reaction == "biomass") :
            ode_idx = "ODE_NET_VAR_"+cell+"_biomass" ; 
     header.write( comma + ode_idx  )
     comma = ", "
 header.write("};\n\n")

 header.write("const REAL A_DENSITY_INERT[NUM_AGENT_TYPES]={")
 comma = "" 
 for cell in celltypes:
     header.write(comma + str( celltypes[cell]['inert_density']))
     comma = ", "
 header.write("};\n\n")
 
 header.write("const REAL A_DIVISION_RADIUS[NUM_AGENT_TYPES]={")
 comma = ""
 for cell in celltypes:
     header.write(comma + str(celltypes[cell]['divRadius']))
     comma = ","
 header.write("};\n")
 
 header.write("const REAL A_MAX_CELL_RADIUS[NUM_AGENT_TYPES]={")
 comma =""
 for cell in celltypes:
     header.write(comma + str(1.1*celltypes[cell]['divRadius'] )) 
     comma = ","
 header.write("};\n")
 
 header.write("const REAL A_MIN_CELL_RADIUS[NUM_AGENT_TYPES]={")
 comma =""
 for cell in celltypes:
     header.write(comma + str(celltypes[cell]['deathRadius']))
     comma = ", "
 header.write("};\n")

 header.write("const REAL A_MAX_CELL_VOL[NUM_AGENT_TYPES]={")
 comma = "" 
 for cell in celltypes:
    rad = 1.1*celltypes[cell]['divRadius']
    if mydomain['nDim']  == 2 : 
       header.write(comma+str(math.pi*rad*rad*mydomain['resolution_agent']))
    else : header.write(comma+str(4.0*math.pi*rad*rad*rad/3.0 ))
    comma= ", "
 header.write("};\n")
 
 header.write("const REAL A_MIN_CELL_VOL[NUM_AGENT_TYPES]={")
 comma = ""
 for cell in celltypes:
    rad = 1.1*celltypes[cell]['deathRadius']
    if mydomain['nDim']  == 2 :
       header.write(comma+str(math.pi*rad*rad*mydomain['resolution_agent']))
    else : header.write(comma+str(4.0*math.pi*rad*rad*rad/3.0 ))
    comma= ", "
 header.write("};\n\n")

 header.write("const REAL A_DIFFUSION_COEFF_CELLS[NUM_AGENT_TYPES]={")
 comma =""
 for cell in celltypes:
     header.write(comma + str(celltypes[cell]['Dcoef'] ))
     comma = ","
 header.write("};\n")

 # here is the link between the agents and \
 header.write("const REAL A_AGENT_ADHESION_S[NUM_AGENT_TYPES][NUM_AGENT_TYPES]={")
 comma = "" 
 for i in range( len(celltypes) ):
    header.write( comma + "{" )
    comma = ""
    for j in range( len(celltypes ) ):
        header.write( comma + str(  myforces[i][j]['adh_strength']  ))
        comma = ", "
    header.write("}") 
 header.write("};\n")

 header.write("const REAL A_AGENT_SHOVING_SCALE[NUM_AGENT_TYPES]={")
 comma = "" 
 for i in range( len(celltypes) ):
    header.write( comma + str(  myforces[i][i]['shoveFactor']  ))
    comma = ", "
 header.write("};\n") 

 header.write("const REAL A_AGENT_SHOVING_LIMIT[NUM_AGENT_TYPES]={")
 comma = ""  
 for cell in celltypes:
    c_id = int( celltypes[cell]['id']  )
    header.write( comma + str(  myforces[c_id][c_id]['shoveLimit']  ))
    comma = ", "
 header.write("};\n")

 header.write("const REAL A_AGENT_BOND_CREATE_FACTOR[NUM_AGENT_TYPES]={")
 comma = ""  
 for i in range( num_celltypes ) :
    header.write( comma + str( myforces[i][i]['attachCreateFactor']))
    comma = ", "
 header.write("};\n")
 
 header.write("const REAL A_AGENT_BOND_DESTROY_FACTOR[NUM_AGENT_TYPES] ={")
 comma = ""  
 for i in range( num_celltypes ) : 
    header.write( comma + str( myforces[i][i]['attachDestroyFactor']))
    comma = ", "
 header.write("};\n")
 
 header.write("const REAL A_AGENT_BOND_S[NUM_AGENT_TYPES][NUM_AGENT_TYPES]={")
 comma = ""
 for i in range(  num_celltypes  ):
    header.write( comma + "{" )
    comma = ""
    for j in range(  num_celltypes  ):
        header.write( comma + str(  myforces[i][j]['bond_strength']  ))
        comma = ", "
    header.write("}")
 header.write("};\n")

 header.write("const REAL A_AGENT_BOND_BOUNDARY_CREATE[NUM_AGENT_TYPES] ={")
 comma = ""
 for i in range( num_celltypes ) :
    header.write( comma + str( myforces[i][i]['attachToBoundaryCreateFactor']))
    comma = ", "
 header.write("};\n")
 
 header.write("const REAL A_AGENT_BOND_BOUNDARY_DESTROY[NUM_AGENT_TYPES] ={")
 comma = ""
 for i in range( num_celltypes ) :
    header.write( comma + str( myforces[i][i]['attachToBoundaryDestroyFactor']))
    comma = ", "
 header.write("};\n")

 header.write("const REAL A_AGENT_BOND_BOUNDARY_S[NUM_AGENT_TYPES] ={")
 comma = ""
 for i in range( num_celltypes ) :
    header.write( comma + str( myforces[i][i]['tightJunctionToBoundaryStrength']))
    comma = ", "
 header.write("};\n\n")

 
 

 header.write("const REAL UB_FULL_COLONY_VOL_RATIO = 0.4;\n\n")

 header.write("const REAL BASELINE_TIME_STEP_DURATION = ")
 header.write( str( mysimulator['baselinetime'])+";\n") 
 header.write("const S32 NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE= ")
 header.write( str( mysimulator['growthtimestep'])+";\n")

 header.write("const S32 A_NUM_PDE_TIME_STEPS_PER_STATE_AND_GRID_STEP[NUM_DIFFUSIBLE_ELEMS]={")
 comma = "" 
 for sol in diffusibles:
    header.write( comma + str( mygridsolver['numofsolversteps']) )
    comma = ", "
 header.write("};\n\n")

 header.write("const S32 SPHERE_UB_VOL_OVLP_RATIO_MAX_LEVEL = 1;\n")

 header.write("const REAL A_MG_NORM_THRESHOLD[NUM_DIFFUSIBLE_ELEMS] = {")
 comma = ""  
 for sol in diffusibles:
     header.write(comma +  "1e-19" )
     comma = ", "
 header.write("};\n\n")


 header.write("const REAL U_SCALE_MAX_INC_RATIO = 1.05;\n")
 header.write("const REAL U_SCALE_LIMIT_THRESHOLD = 0.0005;\n")
 header.write("const REAL UPTAKE_PCT_INC_RATIO = 1.05;\n")
 header.write("const REAL SECRETION_PCT_CHANGE_RATIO = 1.05;\n\n")

 header.write("const REAL A_FLOATING_NUTRIENTS_MAX_DELTA[NUM_DIFFUSIBLE_ELEMS]={")
 comma = ""  
 for sol in diffusibles:
       header.write(comma + "1e-15" )
       comma = ", "
 header.write("};\n\n")

 # here is the link between the agents and \
 header.write("const S32 A_AGENT_GROWTH_SOURCE[NUM_AGENT_TYPES]={")
 comma ="" 
 for cell in celltypes:
    header.write(comma + "0 " )
    comma = ", "
 header.write("};\n\n")

 header.write("const S32 A_BIOMAS_GRID_STATE_DEPENDENCE[NUM_AGENT_TYPES]={")
 comma ="" 
 for cell in celltypes:
    header.write(comma + "4 "  )
    comma = ", "
 header.write("};\n\n") 

 # ODE variables that will printed in VTK file
 header.write("const S32 AA_INDEX_ODE_OUTPUT[NUM_AGENT_TYPES][NUM_AGENT_OUTPUTS -3]={")
 comma = ""
 for cell in celltypes:
    header.write( comma + "{" )
    comma = ""
    for mol in MolOutput_array:
       if mol in celltypes[cell]['moloutput'] :
          header.write( comma + 'ODE_NET_VAR_'+cell+'_'+mol) 
       else: 
          header.write( comma + '-1' )
       comma = ", "
    header.write("}")
    comma = ","
 header.write("};\n\n") 

 # EXTERNAL PERTURBATION 
 NumEpertur = 0
 for epert in eperturbations:
    NumEpertur = NumEpertur + 1 
 header.write("const S32 NUM_E_PERTURBATIONS="+str(NumEpertur)+";\n")
 
 header.write("const ExtConditions A_E_PERTURBATIONS[NUM_E_PERTURBATIONS]={\n")
 comma = ""
 for name in eperturbations:

       agent=eperturbations[name]['AgentType']
       
       header.write(comma + "{" )    
       header.write(str(eperturbations[name]['xo'])+",")
       header.write(str(eperturbations[name]['xf'])+",")
       header.write(str(eperturbations[name]['yo'])+",")
       header.write(str(eperturbations[name]['yf'])+",")
       header.write(str(eperturbations[name]['zo'])+",")
       header.write(str(eperturbations[name]['zf'])+",")
       header.write(str(eperturbations[name]['TimeStep'])+",")
       if ( agent == "" ): header.write("-1,")
       else: header.write("AGENT_TYPE_"+agent+",")

       if ( eperturbations[name]['variable'] == "biomass" ): 
          header.write("CELL_MODEL_REAL_BIOMAS,")
       elif  ( eperturbations[name]['variable'] == "inert" ):
          header.write("CELL_MODEL_REAL_INERT,")
       else:
          header.write("-1,")
       header.write(str(eperturbations[name]['variable_val'])+",")

       if ( eperturbations[name]['molecule'] == ""):
          header.write("-1,")
       else:
          if ( agent == "" ):
             header.write("-1,")
          else:
             header.write("ODE_NET_VAR_"+agent+"_"+str(eperturbations[name]['molecule'])+",")
       header.write(str(eperturbations[name]['molecule_val'])+"}\n")
       
       comma = ","
        
 header.write("};\n\n") 
       
 header.write("#endif\n")
