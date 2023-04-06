# BioXFEL_Tutorial
This tutorial is intended to refine a protein structure against a difference electron density (DED) map  at 1ps. We will be using reciprocial space refinement against extrapolated strcture factors.
The refined structure will be evaluated by comparing the observed and calculated DED maps using the pearson correlation coefficent (pcc).

The data correspnds to the paper Claesson et al., 2020, https://doi.org/10.7554/eLife.53514

Command to clone the git
git clone git@github.com:madanmx/BioXFEL_Tutorial.git


Tutorial is devided into 2/3? steps.

Step 1: Compute and compare pcc of observed and calculated DED maps. Here the structures (pdbs) where already generated. 
Change directory to the quick_check (cd quick_check)

Files and folders:
dark.pdb           #dark refined pdb
refine.mtz         #refine.mtz with dark phases 
1ps_extra_25.mtz   #extrapolated mtz file with 25 activation factor
1ps.pdb            #File from the literature (Elin et. al.,) 
1ps.pdb_wd.map     #Calculated DED map
calc_dmap_b2.sh    #Script to calculate calculated DED map
mock-dark and mock-dark.f #Fortran programs to be used for calculated DED map
pcc/ 
  1ps.pdb_wd.map   # calculated 1ps DED map 
  1ps_wd.map       # observed 1ps DED map 
  pcc_all.sh       # Script to calculate pcc
  pcc.inp          # input file generted from the script and will be input to the mock-dark
    /pcc_out/      # pcc folder 
      1ps.pdb_wd.map   # calculated 1ps DED map 
      1ps_wd.map       # observed 1ps DED map 
      pcc_all.sh   # pcc script   
      pcc.inp      #input file to the RDiff fortran program
      /pcc_out/    #Directory created from the pcc script    
      pcc_summary.txt #all the calues around selected region of comparision will be stored here
      RDiffCCP4    # Fortran program for pcc calculation
      RDiffCCP4.f
 work/             # All the files from the above script


All the relevant files are in the BioXFEL_directories


Step 2: Refine a better structure

Files and folders required:
CBD_phenix.params        #parameter file for the refinement
/dark/                   #dark structure related files
    dark.pdb             #dark.pdb  
    LBV.cif              #cif used to refine dark
    LBV_modified.cif     #cif with some restratins turned off
    refine.mtz           #dark refine.mtz
/extrapolated_mtz/       #extrapolated map directory  
    1ps_extra_25.mtz
/full_dmaps/             #observed DED maps
    1ps_wd.map        
refine_extra.sh          #script to refine 1ps 
refine_extra_Bfactorless.sh #script to refine 1ps with no B factor refinement 
refine_extra_LBV.sh         #script to refine using LBV.cif without restraints off 
/scripts/                  
    calc_dmap_b2.sh         #See tutorial 1
    mock-dark               #See tutorial 1
    mock-dark.f             #See tutorial 1
    pcc_all.sh              #See tutorial 1 
    RDiffCCP4               #See tutorial 1
    RDiffCCP4.f             #See tutorial 1
