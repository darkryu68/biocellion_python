#!/usr/bin/python
import os
import sys 
import xml.etree.ElementTree as ET

sys.path.insert(0,'./scripts')

from xml.etree.ElementTree import Element, SubElement, Comment
from ElementTree_pretty import prettify
from BiocellionParam  import diffusible_solutes, cell_types, domain_parameters, mechanical_parameters, multigrid_solver_parm, basic_simulation_param, create_e_perturbations 
from ReadcDynomics import read_xml
from WriteBcellHeader import write_biocell_header
from WriteBcellConfig import write_biocell_config
from WriteBcellAgent import write_biocell_agent
from WriteBcellGrid import write_biocell_grid
from WriteBcellXML import write_bcell_xml  
from WriteBcellOtherFiles import write_biocell_otherfiles  


## read parameters from input file.

xmlfilename = './examples/basic-growth-noAgar-2D.xml'
dirname = './examples/basic_growth_noAgar_2D'

#xmlfilename = './examples/basic-growth-noAgar-3D.xml'
#dirname = './examples/basic_growth_noAgar_3D'

#xmlfilename = './examples/basic-growth-glucose-2D.xml'
#dirname = './examples/basic_growth_glucose_2D' 

#xmlfilename = './examples/QS_Growth_2D.xml'
#dirname = './examples/QS_Growth_2D'

#xmlfilename = './examples/InducedWrinkles_2D.xml'
#dirname = './examples/InducedWrinkles_2D'

if not os.path.exists(dirname):
    os.makedirs(dirname)

if not os.path.exists(dirname+'/model'):
    os.makedirs( dirname +'/model' )

if not os.path.exists(dirname+'/output'):
    os.makedirs( dirname +'/output' )


diffusibles = dict()  # list of solutes (names)
celltypes = dict() # list of cell types (names)
myreactions = dict() # list of reactions  
eperturbations = dict() # list of perturbations
myforces = []  # parameters for force computations
mydomain = domain_parameters() # parameters of the domain 
mygridsolver =  multigrid_solver_parm() # parameters of the multigrid solver
mysimulator = basic_simulation_param() # timings and other basic arameters of simulation



# read xml from cDynomics
read_xml(diffusibles, celltypes, myreactions, myforces, eperturbations, mydomain, mygridsolver, mysimulator, xmlfilename, dirname+'/model')

print " -----  Cell Types -------------- "
print  celltypes
print " -----  Reactions --------------- "
print myreactions
print " -----  Solutes ----------------- "
print diffusibles 
print " -----  Forces ------------------ "
print myforces 
print "------  E Perturbation ---------- "
print eperturbations 

# 2D case
if ( mydomain['nDim'] == 2 ) :
   if ( mydomain['nz'] == 1 ) :
       mydomain['nz'] = 4 
   

write_biocell_header(diffusibles, celltypes, myreactions, myforces, eperturbations, mydomain, mygridsolver, mysimulator, dirname+'/model' )
write_biocell_config(diffusibles, celltypes, myreactions, myforces, mydomain, mygridsolver, mysimulator, dirname+'/model')
write_biocell_agent(diffusibles, celltypes, myreactions, myforces, mydomain, mygridsolver, mysimulator, dirname+'/model')
write_biocell_grid(diffusibles, celltypes, myreactions, myforces, mydomain, mygridsolver, mysimulator, dirname+'/model')
write_bcell_xml(diffusibles, celltypes, myreactions, myforces, mydomain,mygridsolver, mysimulator, dirname+'/model')
write_biocell_otherfiles(diffusibles, celltypes, myreactions, myforces, mydomain, mygridsolver, mysimulator, dirname+'/model')


