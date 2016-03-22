import os
import sys 
from BiocellionParam  import diffusible_solutes, cell_types, domain_parameters, mechanical_parameters, multigrid_solver_parm, basic_simulation_param 

def write_biocell_agent( diffusibles, celltypes, myreactions, myforces, mydomain, mygridsolver, mysimulator, directory ):

 # write info of the agents
 agentf = open(directory+"/model_routine_agent.cpp", 'w')
 input_agentf  = open('template/model_routine_agent.cpp', 'r')

 for line in input_agentf:
     agentf.write(line)

 agentf.write("\n") 
 agentf.write("void ModelRoutine::spAgentCRNODERHS( const S32 odeNetIdx, const VIdx& vIdx, const SpAgent& spAgent, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox, const Vector<double>& v_y, Vector<double>& v_f ){\n\n" )

 agentf.write("agentType_t type = spAgent.state.getType() ;\n\n");

 # Print the diffusibles
 for solute in diffusibles :
     agentf.write("REAL "+solute+" = spAgent.state.getModelReal(CELL_MODEL_REAL_"+solute+"_AVG );\n") 
 agentf.write("\n")

 for cell in celltypes:
    agentf.write("if ( type == AGENT_TYPE_"+cell+"){\n")

    # Print the molecules ODEs
    for mol in  celltypes[cell]['molecules'] :
        agentf.write("\tREAL "+mol+" = v_y[ODE_NET_VAR_"+cell+"_"+mol+"];\n")
    agentf.write("\n") 
   
    # Print each reaction factor
    for rfactor in myreactions:

       agentf.write("\tREAL r_"+rfactor+"=")

       muMax = myreactions[rfactor]['muMax'] 
       if  ( muMax == 0.0 ) :
          agentf.write(' 0.0 ; \n ') 
          continue
       else :
          agentf.write( str(muMax) )
    
       for MondEq in myreactions[rfactor]['MonodKinetic']:
          agentf.write("*MonodEquation("+str(MondEq['Ks'])+","+ MondEq['solute']+")")  
       for SimpleInh  in myreactions[rfactor]['SimpleInhibition']:
          agentf.write("*SimpleInhibition("+str(SimpleInh['Ki'])+","+ SimpleInh['solute']+")")  

       if ( myreactions[rfactor]['catalyzedby'] == "" ) :
          agentf.write(";\n")
       else:
          agentf.write( "*"+ myreactions[rfactor]['catalyzedby']+";\n")  
    
    # Print the eqution for each concentration rate (v_f)        
    agentf.write("\n")
    for mol  in celltypes[cell]['molecules'] :
       equation_f = ""
       sign = ""
       for rfactor in celltypes[cell]['reactions']   : 
           for i in range(0, len( myreactions[rfactor]['yields']) ):
              if ( myreactions[rfactor]['yields'][i] == mol ):
                yieldf = str(myreactions[rfactor]['yieldFactors'][i])
                equation_f += sign + yieldf + "*r_"+rfactor
                sign = "+" 
                break;
 
       if ( equation_f == "" ):
          sys.exit("ERROR : no rate definition for" + mol + "\n")
       else:
          agentf.write("\tv_f[ODE_NET_VAR_"+cell+"_"+mol+"]="+ equation_f+";\n")
    
    agentf.write("\n")
    agentf.write("}\n") 
  
 agentf.write("\nreturn;\n")
 agentf.write("}\n\n")

