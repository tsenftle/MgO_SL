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
# npj Comput Mater 6, 102 (2020) 
# https://rdcu.be/b5KXZ
#
# Department of Chemical and Biomolecular Engineering and Department of Statistics, Rice University, Houston, TX 77005 (USA)
# *meng@rice.edu
# *tsenftle@rice.edu


The Code is devided into four parts - 

[Step 1]
Directory: ./Step_1/dopant
CODE: Main_MgO_dopant.m
INPUT: BindingEnergy_dopant_V1.xlsx
OUTPUT: data_set_dopant.mat
----------------------------------------------------------------------------------------------------------------------------------------
Directory: ./Step_1/adsorbate
CODE: Main_MgO_adsorbate.m
INPUT: BindingEnergy_adsorbate_V1.xlsx
OUTPUT: data_set_adsorbate.mat
========================================================================================================================================
[Step 2]
Directory: ./Step_2/dopant
CODE: SL script - dopant.R
INPUT: ../Step_1/dopant/data_set_dopant.mat
OUTPUT: (index of selected descriptors from each methods)
MgO_dopant_LS.txt (Descriptors selected by LASSO using dopant-modified MgO data)
MgO_dopant_HS.txt (Descriptors selected by Horseshoe prior using dopant-modified MgO data)
MgO_dopant_DL.txt (Descriptors selected by Dirichlet-Laplace prior using dopant-modified MgO data)
----------------------------------------------------------------------------------------------------------------------------------------
Directory: ./Step_2/adsorbate
CODE: SL script - adsorbate.R
INPUT: ../Step_1/adsorbate/data_set_adsorbate.mat
OUTPUT: (index of selected descriptors from each methods)
MgO_adsorbate_LS.txt (Descriptors selected by LASSO using adsorbate-modified MgO data)
MgO_adsorbate_HS.txt (Descriptors selected by Horseshoe prior using adsorbate-modified MgO data)
MgO_adsorbate_DL.txt (Descriptors selected by Dirichlet-Laplace prior using adsorbate-modified MgO data)
========================================================================================================================================
[Step 3]
Directory: ./Step_3
CODE: Analysis_MgO_dopant_TM.m
INPUT:
../Step_1/dopant/data_set_dopant.mat
../Step_2/dopant/MgO_dopant_LS.txt (Feature space selected from step 2, LASSO in this case)
OUTPUT: analysis in RMSE and R^2
----------------------------------------------------------------------------------------------------------------------------------------
Directory: ./Step_3
CODE: Analysis_MgO_adsorbate_TM.m
INPUT:
../Step_1/adsorbate/data_set_adsorbate.mat
../Step_2/adsorbate/MgO_adsorbate_LS.txt (Feature space selected from step 2, LASSO in this case)
OUTPUT: analysis in RMSE and R^2
========================================================================================================================================
[Step 4]
Directory: ./Step_4/BaO
CODE: Main_BaO.m
INPUT: BindingEnergy_dopant_BaO_V1.xlsx
OUTPUT: data_set_BaO.mat
----------------------------------------------------------------------------------------------------------------------------------------
Directory: ./Step_4/BaO
CODE: Analysis_BaO.m
INPUT:
./data_set_BaO.mat
../Step_2/dopant/MgO_dopant_LS.txt (Feature space selected from step 2, LASSO in this case)
OUTPUT: analysis in RMSE and R^2
----------------------------------------------------------------------------------------------------------------------------------------
Directory: ./Step_4/CaO
CODE: Main_CaO.m
INPUT: BindingEnergy_dopant_CaO_V1.xlsx
OUTPUT: data_set_CaO.mat
----------------------------------------------------------------------------------------------------------------------------------------
Directory: ./Step_4/CaO
CODE: Analysis_CaO.m
INPUT:
./data_set_CaO.mat
../Step_2/dopant/MgO_dopant_LS.txt (Feature space selected from step 2, LASSO in this case)
OUTPUT: analysis in RMSE and R^2
----------------------------------------------------------------------------------------------------------------------------------------
Directory: ./Step_4/ZnO
CODE: Main_ZnO.m
INPUT: BindingEnergy_dopant_ZnO_V1.xlsx
OUTPUT: data_set_ZnO.mat
----------------------------------------------------------------------------------------------------------------------------------------
Directory: ./Step_4/ZnO
CODE: Analysis_ZnO.m
INPUT:
./data_set_ZnO.mat
../Step_2/dopant/MgO_dopant_LS.txt (Feature space selected from step 2, LASSO in this case)
OUTPUT: analysis in RMSE and R^2
----------------------------------------------------------------------------------------------------------------------------------------
