#This is used for both comvert SUSSING->LGAL , LGAL->NIFTY


# AHF table  -  nifty2014 format
AHFdir = "/home/nifty2014/SemiAnalytics/DATA/125Mpc_Euclid/"
# prefix for AHF halo table
AHFprefix = "sussing"
# tree file - Sussing2013 format
SUSSINGtree = "/home/nifty2014/SemiAnalytics/DATA/125Mpc_Euclid/sussing_tree.list"
# snapshot info - Sussing2013 format
SNAPfile = "/home/nifty2014/SemiAnalytics/DATA/125Mpc_Euclid/sussing_index.lst"
# trees file
FileOut = "./treedata/trees_107.0"
# tree_dbids file
FileOut2 = "./treedata/tree_dbids_107.0"
halostruct_file = "halostruct.py"
# model to estimate spin
# 99 ; Bullock2001 spin parameter = 0.02
spin_model = 0
# nifty2014 - force to use wrong mass definition
nifty_forcemass = 0
# prefix for nifty output
prefix = "nifty_II"
# output folder for nifty format
output_folder = "/mnt/lustre/scratch/cs390/nIFTy/MassCompare/force_M200b/SAM/"
#
lgaltree_output = "/mnt/lustre/scratch/cs390/nIFTy/MassCompare/force_M200b/outputs/"


# AHF table  -  nifty2014 format 
AHFdir = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm_II/SUBFIND/"
# prefix for AHF halo table
AHFprefix = "62.5_dm"   
# tree file - Sussing2013 format
SUSSINGtree = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm_II/MergerTree+SUBFIND.txt" 
# snapshot info - Sussing2013 format 
SNAPfile = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm_II/snapidzred.txt"       
# trees file 
FileOut = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm_II/treedata/trees_061.0"  
# tree_dbids file
FileOut2 = "/mnt/lustre/scratch/cs390/nIFTy/62.5_dm_II/treedata/tree_dbids_061.0" 
# model to estimate spin
# 99 ; Bullock2001 spin parameter = 0.02
spin_model = 0
# nifty2014 
force_mass
