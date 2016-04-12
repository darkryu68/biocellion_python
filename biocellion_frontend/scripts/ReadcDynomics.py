import sys 
import xml.etree.ElementTree as ET
import random
from xml.etree.ElementTree import Element, SubElement, Comment
from ElementTree_pretty import prettify
from BiocellionParam  import diffusible_solutes, cell_types, domain_parameters, mechanical_parameters, multigrid_solver_parm, basic_simulation_param 
from BiocellionParam import create_reaction, CreateMonodKinetic, CreateSimpleInhibition, CreateBinding, create_e_perturbations

def read_xml( diffusibles, celltypes, myreactions, myforces, eperturbations,  mydomain, mygridsolver, mysimulator, xmlfilename, directory ):


 ## read parameters from input file.
 tree = ET.parse(xmlfilename)
 root = tree.getroot()

 # Find the number of solutes
 bcell_num_diffusibles = 0;
 for solute in root.findall('solute'):
     name = solute.get('name')

     if ( name == 'pressure' ):
        continue;

     diffusibles[ name ] = diffusible_solutes();
     bcell_num_diffusibles = bcell_num_diffusibles + 1

 # Find the number cell types
 for species in root.findall('species'):
     name = species.get('name')
     celltypes[ name ] = cell_types()
  
 # assign an id to the cell types
 bcell_num_celltypes = 0
 for cell  in celltypes :
     celltypes[cell]['id'] = bcell_num_celltypes
     bcell_num_celltypes = bcell_num_celltypes + 1  
  
 # Initialize the mechanical interactions
 for i in range(0,bcell_num_celltypes):
     myforces.append([])
     for j in range(0,bcell_num_celltypes):
        myforces[i].append( mechanical_parameters()  )

 # Read information of reacions 
 bcell_num_reactions =0
 for reaction in root.findall('reaction'):
     name = reaction.get('name')
     bcell_num_reactions = bcell_num_reactions + 1
     myreactions[name] = create_reaction()
     myreactions[name]['catalyzedby'] = reaction.get('catalyzedBy')
     if ( reaction.get('class')=="ReactionFactorFlux"):
        myreactions[name]['FluxFlag'] = True;

     for param in reaction.iter('param'):
         if ( param.get('name')=="muMax"):
             #myreactions[name]['muMax'] = float(param.text)/3600.0 # s^-1 
             myreactions[name]['muMax'] = float(param.text) # s^-1 
     for yields in reaction.iter('yield'):
         for param in yields.iter('param'):
             myreactions[name]['yields'].append(param.get('name'))
             myreactions[name]['yieldFactors'].append(float(param.text))
     for kinetic in reaction.iter('kineticFactor'):
         if ( kinetic.get('class') == "MonodKinetic") :
            temp = CreateMonodKinetic()
            temp['solute'] = kinetic.get('solute') 
            for param in kinetic.iter('param'):
               if ( param.get('name') == "Ks"  ):
                  temp['Ks'] = float(param.text)
            myreactions[name]['MonodKinetic'].append(temp)
         elif ( kinetic.get('class') == "FirstOrderKinetic"):
            temp = CreateMonodKinetic()
            myreactions[name]['MonodKinetic'].append(temp)
         elif ( kinetic.get('class') == "Binding" ): 
            temp = CreateBinding()
            temp['solute'] = kinetic.get('solute') 
            myreactions[name]['Binding'].append(temp)
         else :
            sys.exit("ERROR : reaction factor not supperted  ")
                 


 # read timestep from cDynomics
 cDynomicTime = 0
 for simulator in root.findall('simulator'):
    for param in simulator.iter('param'):
       if ( param.get('name')=="agentTimeStep"):
          cDynomicTime = round( float(param.text) * 3600)
          mysimulator['baselinetime'] = cDynomicTime
          #print "SEEE  " ,  mysimulator['baselinetime'], "\n"  
    if ( cDynomicTime == 0 ):
       sys.exit("ERROR : reading agentTimeStep")   
     
    for param in simulator.iter('param'):
       if ( param.get('name')=="outputPeriod"):
          temp = float(param.text)*3600 / mysimulator['baselinetime']  
          if ( temp >= 1 ) :
             mysimulator['outputPeriod']=round(temp)
       elif ( param.get('name')=="endOfSimulation"): 
          temp = float(param.text) * 3600
          mysimulator['numbersteps']=round(temp/cDynomicTime)
       elif ( param.get('name') == "restartPreviousRun" ) :
          if ( param.text == "true" ) :
             sys.exit("ERROR : restartPreviousRun not implemented")
         
 # Input section is not implemented
  
 # read information from solutes  
 if (len(diffusibles)  >  0 ):
    for solute in root.findall('solute'):
       name = solute.get('name')
       if ( name == 'pressure' ):
           continue;

       for param in solute.findall('param'):
          if ( param.get('name') == "diffusivity" ):
              factor = 1.0; # change units to um^2 /s
              if ( param.get('unit') == "m2.day-1" ):
                  factor = 1.15741e7 
              value =float(param.text) * factor 
              diffusibles[name]['Dcoef_cd'] = value
              diffusibles[name]['Dcoef_agar'] = value
              if ( param.get('colonyonly') == "true" ):
                  diffusibles[name]['colonyonly'] = True  

 # read information from particles (density)
 density_biomass = -1.0
 density_inert  = -1.0 
 for particle in root.findall('particle'):
    if ( particle.get('name') == "biomass" ) :
       for param in particle.findall('param'):
          if ( param.get('name') == "density" ):
             density_biomass = float(param.text)
    elif ( particle.get('name') == "inert" ):
       for param in particle.findall('param'):
          if ( param.get('name') == "density" ) :
             density_inert = float(param.text)
    else :
       sys.exit("ERROR :particle name unknow ")
       
 if ( density_biomass == -1.0):
    sys.exit("ERROR :density of biomass not found ")
 
 if ( density_inert  == -1.0):
    sys.exit("ERROR :density of inert not found ")
 
 for solute in celltypes: 
    celltypes[solute]['biomass_density'] = density_biomass  
    celltypes[solute]['inert_density'] = density_inert  
   
 # Information from Agar 
 # where this can be specified?? 
 for agar in root.iter('agar'):
    for param in agar.findall('param'):
       if ( param.get('name') == "depth" ):
          mydomain['agar_heigth'] = int( 4*round( int(param.text)/4.0 + 0.25)) 
       elif ( param.get('name') == "width" ):
          mydomain['ny'] = int( 4*round( int(param.text)/4.0 + 0.25))
       elif ( param.get('name') == "height" ):
          mydomain['nz'] = int( 4*round( int(param.text)/4.0 + 0.25))
       elif ( param.get('name') == "concentration"):
          for solute in diffusibles :
             diffusibles[solute]['ini_con_agar']=float(param.text) 
           
 # Read information of Computation  Domain 
 for domain in root.iter('computationDomain'):
    for grid in domain.findall('grid'):
       mydomain['nx'] = int(4*round( int(grid.get('nI'))/4.0  + 0.25))
       mydomain['ny'] = int(4*round( int(grid.get('nJ'))/4.0  + 0.25))
       mydomain['nz'] = int(4*round( int(grid.get('nK'))/4.0  + 0.25))
       mydomain['nDim'] = int(grid.get('nDim'))
      
    for param in domain.findall('param'):  
       
        if ( param.get('name') == "resolution" ):  
           mydomain['resolution'] = float(param.text) 
        elif ( param.get('name') == "biofilmDiffusivity"): 
           mydomain['biofilmDiffusivity'] = float(param.text) 
            
 for agentgrid  in root.iter('agentGrid'): 
    for param in agentgrid.findall('param'):
       if ( param.get('name') == "resolution" ):  
           mydomain['resolution_agent'] = float(param.text) 
 if ( mydomain['resolution_agent'] == 0.0 ):
    mydomain['resolution_agent'] = mydomain['resolution']

 # read boundar conditions 
 for bc in root.iter('boundaryCondition'):
    if (bc.get('name') == "yNz" ):
       if ( bc.get('class') == "BoundaryZeroFlux" ):
          mygridsolver['bcType_x+']='BC_TYPE_NEUMANN_CONST'
          mygridsolver['bcVal_x+']= 0.0
          mygridsolver['dbtype_x'] = 'HARD'
       elif  ( bc.get('class') == "BoundaryZeroFluxDisappear"):
          mygridsolver['bcType_x+']='BC_TYPE_NEUMANN_CONST'
          mygridsolver['bcVal_x+']= 0.0
          mygridsolver['dbtype_x'] = 'DISAPPEAR'
    elif (bc.get('name') == "y0z" ):
       if ( bc.get('class') == "BoundaryZeroFlux" ):
          mygridsolver['bcType_x-']='BC_TYPE_NEUMANN_CONST'
          mygridsolver['bcVal_x-']= 0.0
          mygridsolver['dbtype_x'] = 'HARD'
       elif ( bc.get('class') == "BoundaryZeroFluxDisappear"):
          mygridsolver['bcType_x-']='BC_TYPE_NEUMANN_CONST'
          mygridsolver['bcVal_x-']= 0.0
          mygridsolver['dbtype_x'] = 'DISAPPEAR'
    elif (bc.get('name') == "xNz" ):
       if ( bc.get('class') == "BoundaryZeroFlux" ):
          mygridsolver['bcType_y+']='BC_TYPE_NEUMANN_CONST'
          mygridsolver['bcVal_y+']= 0.0
          mygridsolver['dbtype_y'] = 'HARD'
       elif ( bc.get('class') == "BoundaryZeroFluxDisappear"):
          mygridsolver['bcType_y+']='BC_TYPE_NEUMANN_CONST'
          mygridsolver['bcVal_y+']= 0.0
          mygridsolver['dbtype_y'] = 'DISAPPEAR'
    elif (bc.get('name') == "x0z" ):
       if ( bc.get('class') == "BoundaryZeroFlux" ):
          mygridsolver['bcType_y-']='BC_TYPE_NEUMANN_CONST'
          mygridsolver['bcVal_y-']= 0.0
          mygridsolver['dbtype_y'] = 'HARD'
       elif ( bc.get('class') == "BoundaryZeroFluxDisappear"):
          mygridsolver['bcType_y-']='BC_TYPE_NEUMANN_CONST'
          mygridsolver['bcVal_y-']= 0.0
          mygridsolver['dbtype_y'] = 'DISAPPEAR'
    elif (bc.get('name') == "xNy" ):
       if ( bc.get('class') == "BoundaryZeroFlux" ):
          mygridsolver['bcType_z+']='BC_TYPE_NEUMANN_CONST'
          mygridsolver['bcVal_z+']= 0.0
          mygridsolver['dbtype_z'] = 'HARD'
       elif ( bc.get('class') == "BoundaryZeroFluxDisappear"):
          mygridsolver['bcType_z+']='BC_TYPE_NEUMANN_CONST'
          mygridsolver['bcVal_z+']= 0.0
          mygridsolver['dbtype_z'] = 'DISAPPEAR'
    elif (bc.get('name') == "x0y" ):
       if ( bc.get('class') == "BoundaryZeroFlux"):
          mygridsolver['bcType_z-']='BC_TYPE_NEUMANN_CONST'
          mygridsolver['bcVal_z-']= 0.0
          mygridsolver['dbtype_z'] = 'HARD'
       elif ( bc.get('class') == "BoundaryZeroFluxDisappear"):
          mygridsolver['bcType_z-']='BC_TYPE_NEUMANN_CONST'
          mygridsolver['bcVal_z-']= 0.0
          mygridsolver['dbtype_z'] = 'DISAPPEAR'


 #solver not implemented yet
 cells_file  = open( directory + '/cells.txt', 'w')
 for species in root.findall('species'):
    name = species.get('name')
    type_id = celltypes[ name ]['id']

    for param in species.iter('param'):
       if ( param.get('name') == "divRadius" ):
          celltypes[name]['divRadius']=float(param.text)
       elif ( param.get('name') == "deathRadius" ):
          celltypes[name]['deathRadius']=float(param.text)    
       elif ( param.get('name') == "shoveFactor" ):
          myforces[type_id][type_id]['shoveFactor'] = float(param.text)
       elif ( param.get('name') == "shoveLimit" ):
          myforces[type_id][type_id]['shoveLimit'] = float(param.text)
       elif ( param.get('name') == "attachCreateFactor" ):
          myforces[type_id][type_id]['attachCreateFactor'] = float(param.text)
       elif ( param.get('name') == "attachDestroyFactor" ):
          myforces[type_id][type_id]['attachDestroyFactor'] = float(param.text)
       elif ( param.get('name') == "attachToBoundaryCreateFactor" ):
          myforces[type_id][type_id]['attachToBoundaryCreateFactor'] = float(param.text)
       elif ( param.get('name') == "attachToBoundaryDestroyFactor" ):
          myforces[type_id][type_id]['attachToBoundaryDestroyFactor'] = float(param.text)
       elif ( param.get('name') == "tightJunctionToBoundaryStrength" ):
          myforces[type_id][type_id]['tightJunctionToBoundaryStrength'] = float(param.text)  
          

       if ( mydomain['nDim'] == 2 ) : # in cdynamics , 2 dimnesional siulations uses cylinder
          for agents in celltypes:
            celltypes[agents]['shape']='cylinder'

    # mechanical parameters   
    for adm in species.iter('adhesion'):
        name_j  =  adm.get('withSpecies')
        id_j = celltypes[name_j]['id'] 
        myforces[type_id][id_j]['adh_strength']=float(adm.get('strength'))
   
    for bond in species.iter('tightJunction'):
        name_j = bond.get('withSpecies') 
        ######## TODO #### check if name is available 
        id_j = celltypes[name_j]['id'] 
        myforces[type_id][id_j]['bond_strength']=float(bond.get('stiffness'))
 
    # Aff a list of reactions from xml
    for reaction in species.findall('reaction'):
        ReactionName = reaction.get('name')  
        celltypes[name]['reactions'].append(ReactionName)


    # Initial conditions
    # read information from particles (density)
    biomass_ini = 0.0
    inert_ini  = 0.0
    for particle in species.findall('particle'):
       if ( particle.get('name') == "biomass" ) :
          for param in particle.findall('param'):
             if ( param.get('name') == "mass" ):
                biomass_ini = float(param.text)
       elif ( particle.get('name') == "inert" ):
          for param in particle.findall('param'):
             if ( param.get('name') == "mass" ) :
                inert_ini = float(param.text)
       else :
          sys.exit("ERROR :particle name unknow ")
     
     
    # read initial configuration of cells and write it in cells.txt
    for area in species.findall('initArea'):
       x_coor = [-1.0 , -1.0]
       y_coor = [-1.0 , -1.0]
       z_coor = [-1.0 , -1.0]
       NumCoordinates = 0
       for coordinates in area.findall('coordinates'):
          NumCoordinates = NumCoordinates + 1
          if (NumCoordinates >  2):
             sys.exit("ERROR :  wrong number of coordinate in initArea ") 
          x_coor[NumCoordinates-1] = float(coordinates.get('x'))
          y_coor[NumCoordinates-1] = float(coordinates.get('y'))
          z_coor[NumCoordinates-1] = float(coordinates.get('z'))
       
       Ncellsx = 0;
       Ncellsy = 0;
       Ncellsz = 0;
       for blocks in area.findall('blocks'):
          Ncellsx = int(blocks.get('rows'))
          Ncellsy = int(blocks.get('cols'))
          Ncellsz = int(blocks.get('bars',default=1))
      
       NumCells = Ncellsx * Ncellsy * Ncellsz  
       if ( NumCells > 0 ) :
          
          cells_file.write( str(NumCells) +"\n"  )
          h_x = abs(x_coor[1] - x_coor[0]) / Ncellsx 
          h_y = abs(y_coor[1] - y_coor[0]) / Ncellsy 
          h_z = abs(z_coor[1] - z_coor[0]) / Ncellsz 

          for i in range(0,Ncellsx ):
             for j in range(0,Ncellsy ):
                for k in range(0,Ncellsz ):
                   xx = x_coor[0] +  h_x/2 + i*h_x 
                   yy = y_coor[0] +  h_y/2 + j*h_y
                   zz = z_coor[0] +  h_z/2 + k*h_z
            
                   if ( mydomain['nDim'] == 2 ):
                       zz = 0.5 * mydomain['resolution']  
        

                   value = ( xx, yy, zz, biomass_ini, inert_ini, celltypes[ name ]['id'] )
                   mybuffer = "%f %f %f %f %f %d\n" % value
                   cells_file.write( mybuffer  ) 
       else:  
          NumCells = int( area.get('number') )

          cells_file.write( str(NumCells) +"\n" )
          # how to put cells randomly on region especified by coordinate
          dx =  abs(x_coor[1] - x_coor[0]) 
          dy =  abs(y_coor[1] - y_coor[0]) 
          dz =  abs(z_coor[1] - z_coor[0])
             
          for i in range(NumCells): 
              xx = 0.5*( x_coor[0] + x_coor[1] ) + dx*(random.random() - 0.5)  
              yy = 0.5*( y_coor[0] + y_coor[1] ) + dy*(random.random() - 0.5)  
              zz = 0.5*( z_coor[0] + z_coor[1] ) + dz*(random.random() - 0.5) 
   
              if ( mydomain['nDim'] == 2 ):
                 zz = 0.5 * mydomain['resolution']  
        
              value = ( xx, yy, zz, biomass_ini, inert_ini, celltypes[ name ]['id'] )
              mybuffer = "%f %f %f %f %f %d\n" % value
              cells_file.write( mybuffer  )     
                                     
 cells_file.close( )

 for species in root.findall('species'):
    name = species.get('name')
    type_id = celltypes[ name ]['id']

    for reaction in species.findall('reaction'):
       ReactionName = reaction.get('name')
       for mol in myreactions[ReactionName]['yields'] :

          if ( diffusibles.has_key(mol)  ):
              continue
          found = False
          for inmol in celltypes[name]['molecules'] :
             if ( inmol == mol ):
                found = True
                continue
               
          if ( not found ) : 
              celltypes[name]['molecules'].append(mol)
     
    for output in species.findall('moloutput'):
        mol_name = output.get('name')
        def_value = output.get('default')
        celltypes[name]['moloutput'].append( mol_name)
        celltypes[name]['moloutput_default'].append(float(def_value))
 
    for conditions in species.findall('entryConditions'):
        x_coor = [-1.0 , -1.0]
        y_coor = [-1.0 , -1.0]
        z_coor = [-1.0 , -1.0]
        Time = -1
        variable =""
        var_value = -1.0
        molecule =""
        mol_value = -1.0
        NumCoordinates = 0
        for condition in conditions.findall('entryCondition'):
           if ( condition.get('type')=="location"):
              for coordinates in condition.findall('coordinates'):
                 NumCoordinates = NumCoordinates + 1
                 if (NumCoordinates >  2):
                    sys.exit("ERROR: wrong number of coordinates ")
                 x_coor[NumCoordinates-1] = float(coordinates.get('x'))
                 y_coor[NumCoordinates-1] = float(coordinates.get('y'))
                 z_coor[NumCoordinates-1] = float(coordinates.get('z'))
           elif ( condition.get('type')=="timing"): 
              for param in condition.findall('param'):
                 if ( param.get('name') == "biomass" ):
                    var_value = float(param.text)
                    variable = "biomass"
                 elif ( param.get('name') == "time" ):
                    if ( param.get('unit') == "hour" ):
                       Time = int (float(param.text)*3600/mysimulator['baselinetime'])
                    elif( param.get('unit') == "hour" ):
                       Time = int  (float(param.text)/mysimulator['baselinetime'])
                    else :
                       Time = int(param.text)
           
        if (( NumCoordinates == 2) and ( Time > 0 )):
            eperturbations[name] = create_e_perturbations();
            eperturbations[name]['xo'] = x_coor[0] 
            eperturbations[name]['xf'] = x_coor[1] 
            eperturbations[name]['yo'] = y_coor[0] 
            eperturbations[name]['yf'] = y_coor[1] 
            eperturbations[name]['zo'] = z_coor[0] 
            eperturbations[name]['zf'] = z_coor[1] 
            eperturbations[name]['TimeStep'] = Time 
            eperturbations[name]['AgentType'] = name  
            if ( variable == "biomass") :
               eperturbations[name]['variable'] = 'biomass'
               eperturbations[name]['variable_val'] = var_value 




        

 #  Reabtions[ name ] reactions inside every cell 
 

# This is is commented becuase state and grid will be 
# updated at every baseline step. Biocellion will find the steady-state
# solution of the diffusion equations and will be update every baseline
# time step
 
 # solver used to get time steps tu update cell states
 #for solver in root.findall('solver'):
 #   if ( solver.get('name') == "solutes" ) :
 #      for param in solver.iter('param'):
 #         if ( param.get('name') == "preStep" ):
 #            mysimulator['growthtimestep'] = int(param.text) 



