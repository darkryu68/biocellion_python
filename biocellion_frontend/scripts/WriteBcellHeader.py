import os
import sys 
from BiocellionParam  import diffusible_solutes, cell_types, domain_parameters, mechanical_parameters, multigrid_solver_parm, basic_simulation_param 
from BiocellionParam import create_reaction, CreateMonodKinetic, CreateSimpleInhibition

def write_biocell_header( diffusibles, celltypes, myreactions, myforces, mydomain, mygridsolver, mysimulator, directory ):

 num_celltypes  = len( celltypes )

 # write the files for biocellion
 header = open(directory+"/model_define.h", 'w')
 header.write('#ifndef __MY_DEFINE_H__' +"\n" )
 header.write('#define __MY_DEFINE_H__' +"\n\n")
 header.write('#include "biocellion.h"'+"\n\n")
 header.write('#define NUM_JUNCTION_END_TYPES 1 /* we consider only one junction end type */ '+"\n")
 header.write('#define IF_GRID_SPACING_SEL ' )
 header.write( str(mydomain['resolution']) + "\n")
 header.write('#define SYSTEM_DIMENSION ')
 header.write( str(mydomain['nDim']) + "\n\n")

 header.write("typedef enum _agent_type_e {  \n"  )
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
    if ( len( celltypes[cell]['molecules'] ) ) > 0 :
       header.write("typedef enum _ode_net_"+ cell +"_var_e {   \n")
       for mol in celltypes[cell]['molecules']: 
          header.write("\tODE_NET_VAR_"+cell+"_"+mol+",\n")
       header.write("\tNUM_ODE_NET_VAR_"+ cell +"\n")
       header.write("} ode_net_"+cell+"_var_e;\n\n")

 header.write("typedef enum _alive_cell_model_int_e { \n")
 header.write("\tCELL_MODEL_INT_BOND_B, \n")
 header.write("\tNUM_MODEL_INTS \n")
 header.write("} alive_cell_model_int_e;\n\n")

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

 header.write("class UBInitData {\n")
 header.write("public:\n")
 header.write("\tS32 numCells ;  \n")
 header.write("\tREAL x_offset; \n")
 header.write("\tREAL y_offset; \n")
 header.write("\tREAL z_offset; \n")
 header.write("\tS32 a_type; \n")
 header.write("\tREAL biomass; \n")
 header.write("\tREAL inert; \n")
 header.write("\tS32  IdxIniCellData;\n")
 header.write("};\n\n")

 header.write("extern S32  INI_N_CELLS;\n\n")

 header.write("const REAL MY_PI = 3.14159265358979323846;\n\n")

 header.write("const REAL IF_GRID_SPACING = ")
 header.write( str(mydomain['resolution']) + ";\n\n")

 header.write("const S32 A_NUM_AMR_LEVELS[NUM_DIFFUSIBLE_ELEMS]={") 
 comma = "" 
 for sol in diffusibles:
    header.write(comma + "2" ) 
    comma = ", "
 header.write("};\n\n")   

 header.write("const S32 AGAR_HEIGHT = ") 
 header.write( str(mydomain['agar_heigth'])+";\n" ) 
 header.write("const S32 AGAR_TOP_BUFFER_HEIGHT = 4;\n");

 header.write("const REAL A_INIT_CONCENTRATION[NUM_DIFFUSIBLE_ELEMS]={")
 comma = "" 
 for sol in diffusibles:
     header.write(comma +str(diffusibles[sol]['ini_con_agar']))
     comma = ", "
     #header.write(","+str(sol['ini_con_agar'])   )
 header.write("};\n\n")

 header.write("const REAL A_DIFFUSION_COEFF_AGAR[NUM_DIFFUSIBLE_ELEMS]={")
 comma = "" ; # looks horrrible
 for sol in diffusibles:
    diff =diffusibles[sol]['Dcoef_agar']*mydomain['biofilmDiffusivity'] * 0.5 
    header.write(comma +str(diff))
    comma = ", "
 header.write("};\n\n")

 header.write("const REAL A_DIFFUSION_COEFF_COLONY[NUM_DIFFUSIBLE_ELEMS]={")
 cont = 0 ; # looks horrrible
 for sol in diffusibles:
    diff =diffusibles[sol]['Dcoef_cd']*mydomain['biofilmDiffusivity']
    if ( cont == 0 ) :
       header.write(str(diff))
       cont = 1
    else :
       header.write(", "+ str(diff))
 header.write("};\n\n")


 # this will be dirty code
 header.write("const REAL AA_SECRETION_RATE_PER_CELL[NUM_AGENT_TYPES][NUM_DIFFUSIBLE_ELEMS]={") 
 coma = "" 
 for cell in celltypes: 
    header.write(coma + "{")
    coma = ""
    for sol in diffusibles:
       header.write(coma + "0.0 ")
       coma = ","
    header.write("} ")
    coma = ","
 header.write("};\n") 

 header.write("const REAL AA_UPTAKE_RATE_PER_CELL[NUM_AGENT_TYPES][NUM_DIFFUSIBLE_ELEMS]={")
 coma="" 
 for cell in celltypes:
    header.write(coma + "{")
    coma = ""
    for sol in diffusibles:
       header.write(coma + "0.0 ")
       coma = ","
    header.write("} ")
    coma = ","
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
 
 header.write("const REAL A_BIOMASS_DEGRADATION[NUM_AGENT_TYPES]={")
 comma = ""
 for cell in celltypes:
     header.write( comma +  "0.0"  )
     comma = ", "
 header.write("};\n\n")

 header.write("const REAL DIVISION_RADIUS = ")
 cell =  celltypes.keys()[0]
 header.write( str(celltypes[cell]['divRadius'])+ ";\n")
 header.write("const REAL MAX_CELL_RADIUS = DIVISION_RADIUS * 1.1;\n")

 header.write("const REAL MIN_CELL_RADIUS = ")
 cell =  celltypes.keys()[0]
 header.write( str(celltypes[cell][ 'deathRadius'  ])+ ";\n")

 if ( celltypes[cell]['shape'] == 'cylinder' ):
    header.write("const REAL MAX_CELL_VOL = (  MY_PI ) * MAX_CELL_RADIUS * MAX_CELL_RADIUS * IF_GRID_SPACING   ;\n")
    header.write("const REAL MIN_CELL_VOL = (  MY_PI ) * MIN_CELL_RADIUS * MIN_CELL_RADIUS * IF_GRID_SPACING ;\n\n")
 else:
    header.write("const REAL MAX_CELL_VOL = (  MY_PI ) * MAX_CELL_RADIUS * MAX_CELL_RADIUS * MAX_CELL_RADIUS ;\n")
    header.write("const REAL MIN_CELL_VOL = (  MY_PI ) * MIN_CELL_RADIUS * MIN_CELL_RADIUS * MIN_CELL_RADIUS ;\n\n")


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
 
 if ( celltypes[cell]['shape'] == 'cylinder' ):
    header.write("inline REAL radius_from_volume( REAL volume ) {\n")  
    header.write("\treturn SQRT( volume/( MY_PI*IF_GRID_SPACING )); \n")
    header.write("}\n\n")
 else :
    header.write("inline REAL radius_from_volume( REAL volume ) {\n")
    header.write("\treturn CBRT( volume * 3.0 / ( 4.0 * MY_PI ) );\n")
    header.write("}\n\n")    

 header.write("inline REAL MonodEquation(  REAL Kc , REAL u ) {\n")
 header.write("\treturn  u / ( Kc + u ) ;\n")
 header.write("}\n\n") 
 # Write the biomas rate equation, it depends of cell type
 #header.write("inline REAL BiomasRate(S32 type, REAL Biomass, REAL uscale, REAL uptake) {\n")  
 #header.write("switch (type)  {\n")
 #for cell in celltypes: 
 #   header.write("case "+ str( celltypes[cell]['id'] ) + ":\n" )
    
 #   name = ""
 #   for reaction_name in celltypes[cell]['reactions']:
 #       if ( myreactions[reaction_name]['yields'] == "biomass" ):
 #            name = reaction_name   
 #   if ( name == "") :
 #      header.write("\treturn 0.0; \n")
 #   else :
 #      if ( len( myreactions[name]['MonodKinetic'] ) > 0 ):
 #         header.write("\treturn " + str(myreactions[name]['muMax']) )
 #         for reaction_i in myreactions[name]['MonodKinetic'] :
 #            if ( reaction_i['Ks'] != 0.0 ) :   
 #               header.write("* ( uscale/("+str(reaction_i['Ks'])+" +uscale))")
 #         header.write("  * Biomass * uptake ; \n")
 #      else :       
 #           header.write("\treturn " + str(myreactions[name]['muMax']) + "; \n")

 #   header.write("\tbreak;\n" ) 

 #header.write("default :\n")
 #header.write("\treturn 0.0; \n")
 #header.write("\tbreak; \n")
 #header.write("};\n\n")
 #header.write("}\n")

 header.write("#endif\n")
