import sys 

# high level parameters of solute diffusables
def diffusible_solutes():
    '''Returns a dictionary that represent a solute '''
    solute = dict()
    solute['Dcoef_cd'] = 0.0  # Diff. coefficient computation domain
    solute['Dcoef_agar'] = 0.0 # Diff. coefficient agar
    solute['degradation'] = 0.0  # degradation rate
    solute['ini_con_cd'] = 0.0 # initial concentration in computation domain
    solute['ini_con_agar'] = 0.0 # initial concentration  agar
    return solute 

# high level parameters for cell types
def cell_types():
    '''Returns a dictionary that represent a cell type '''
    ctype = dict()
    ctype['muMax'] = 0.0  # Maximum glucose uptake rate // more flexiblo to specify the reaction type
    ctype['Ks'] = 0.0  # Monod Kinetic parameter
    ctype['divRadius'] = 0.0  # Maximum radius (used for cell division) 
    ctype['deathRadius'] = 0.0  # Minimum radius (used for cell death)
    ctype['biomass_density'] = 0.0 # density of cell
    ctype['inert_density'] = 0.0 # density of cell
    ctype['shape'] = 'sphere' # spherical shape (we consider only shperical or cylindrical
    ctype['id'] = 0  # identifier for biocellion 
    ctype['reactions'] = []  # set of names of reactions 
    return ctype 

# high level parameters for intracellular reactions
def create_reaction():
    '''Returns a dictionary that represent a reaction '''
    reactionfactor = dict()
    reactionfactor['catalyzedby'] = ''
    reactionfactor['yield'] = ''
    reactionfactor['muMax'] = 0.0  # factor
    reactionfactor['MonodKinetic']  = []  # a list of Monod reactions
    reactionfactor['SimpleInhibition'] = [] # a lis of simple inhibition 
    return reactionfactor

def CreateMonodKinetic():
    '''Returns a dictionary that represent a MonodKinetic '''
    MonodKinetic = dict()
    MonodKinetic['Ks'] = 0.0 
    MonodKinetic['solute'] = '' # name of the solute,  
    return MonodKinetic

def CreateSimpleInhibition():
    '''Returns a dictionary that represent a SimpleInhibition '''    
    SimpleInhibition = dict()
    SimpleInhibition['Ki'] = 0.0
    MonodKinetic['solute'] = '' # name of the solute,  
    
 

# parameters that determine the reaction domain
def domain_parameters():
    domain = dict()
    domain['nx'] = 0  # Number of voxels in x directions 
    domain['ny'] = 0  # Number of voxels in x directions 
    domain['nz'] = 0  # Number of voxels in x directions 
    domain['agar_heigth'] = 0  # idetermine the size of the agar  
    domain['resolution'] = 0.0 # Grid resolution   
    domain['biofilmDiffusivity'] = 0.0  # multiplier of diffusion coenficcient
    domain['nDim'] = 3 # 2 or 3 dimentional  
    return domain

def mechanical_parameters():
    force  = dict() 
    force['shoveFactor'] = 1.0 # Factor of radius  for shoving force
    force['shoveLimit'] = 0.0 # overalp tolerance of shoving force
    force['adh_strength'] = 0.0 # cell adhesion strength
    force['attachCreateFactor'] = 1.0  # Bonds between cells
    force['attachDestroyFactor'] = 1.0
    force['bond_strength'] = 0.0
    force['attachToBoundaryCreateFactor'] = 1.0 # Bonds with surface
    force['attachToBoundaryDestroyFactor']= 1.0
    force['tightJunctionToBoundaryStrength'] = 0.0
    return force

def multigrid_solver_parm():
    solver = dict()
    solver['bcType_x+'] = 'BC_TYPE_NEUMANN_CONST'  # type boundary condition
    solver['bcVal_x+'] = 0.0                       # values  
    solver['bcType_x-'] = 'BC_TYPE_NEUMANN_CONST'  # type boundary condition
    solver['bcVal_x-'] = 0.0                       # value
    solver['bcType_y+'] = 'BC_TYPE_NEUMANN_CONST'  # type boundary condition
    solver['bcVal_y+'] = 0.0
    solver['bcType_y-'] = 'BC_TYPE_NEUMANN_CONST'  # type boundary condition
    solver['bcVal_y-'] = 0.0
    solver['bcType_z+'] = 'BC_TYPE_NEUMANN_CONST'  # type boundary condition
    solver['bcVal_z+'] = 0.0
    solver['bcType_z-'] = 'BC_TYPE_NEUMANN_CONST'  # type boundary condition
    solver['bcVal_z-'] = 0.0
    solver['numofsolversteps'] = 1 # dont understan how cDynomics handles this
    return solver  

def basic_simulation_param():
   simulation = dict()
   simulation['outputPeriod'] = 1 # print files every time
   simulation['baselinetime'] = 1 #one second
   simulation['numbersteps'] = 1 # end simulation after one baseline
   simulation['growthtimestep'] = 1 # 
   return simulation 
