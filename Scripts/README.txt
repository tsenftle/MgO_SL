# Code written by Chun-Yen Liu
# Ph.D. Student
# Department of Chemical and Biomolecular Engineering
# Rice University, Houtson, TX 77005
#
# Adapted from:
# Interaction trends between single metal atoms and oxide supports identified with density functional theory and statistical learning
# Nolan J. O'Connor, A S M Jonayat, Michael J. Janik*, Thomas P. Senftle* 
#
# Distributed as part of the publication - 
# Using Statistical Learning to Predict Interactions Between Single Metal Atoms and Modified MgO(100) Supports
# Chun-Yen Liu, Shijia Zhang, Daniel Martinez, Meng Li*, and Thomas P. Senftle*
#
# Department of Chemical and Biomolecular Engineering and Department of Statistics, Rice University, Houston, TX 77005 (USA)
# *meng@rice.edu
# *tsenftle@rice.edu


The Code is devided into three parts - 


[Step 1]
Directory: ./Step_1
CODE: Main.m

INPUT: BindingEnergyV1.xlsx

OUTPUT: data_set.mat (matlab data file)

[Step 2]
Directory: ./Step_2
CODE: SL script.R
INPUT: ../Step_1/data_set.mat (matlab data file) 
OUTPUT: (index of selected descriptors from each methods)
MgO_dopant_LS.txt (Descriptors selected by LASSO using dopant-modified MgO data)
MgO_dopant_HS.txt (Descriptors selected by Horseshoe prior using dopant-modified MgO data)
MgO_dopant_DL.txt (Descriptors selected by Dirichlet-Laplace prior using dopant-modified MgO data)
MgO_adsorbate_LS.txt (Descriptors selected by LASSO using adsorbate-modified MgO data)
MgO_adsorbate_HS.txt (Descriptors selected by Horseshoe prior using adsorbate-modified MgO data)
MgO_adsorbate_DL.txt (Descriptors selected by Dirichlet-Laplace prior using adsorbate-modified MgO data)

[Step 3]
Directory: ./Step_3
CODE: Analysis_MgO.m

INPUT:
../Step_1/data_set.mat (matlab data file) 
../Step_2/MgO_dopant_LS.txt (Feature space selected from step 2)
OUTPUT: part of Figure 7 and 8 in the manuscript


[Step 4]
Directory: ./Step_4
CODE: Main_CaO.m

INPUT: BindingEnergy_CaO_V1.xlsx

OUTPUT: data_set_CaO.mat (matlab data file) 
CODE: Analysis_CaO.m

INPUT: 
./data_set_CaO.mat
../Step_2/MgO_dopant_LS.txt (Feature space selected from step 2)
OUTPUT: part of Figure 9 and 10 in the manuscript