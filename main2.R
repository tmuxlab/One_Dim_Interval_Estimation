# install.packages("tidyverse")
library(tidyverse)
if (!requireNamespace("peakRAM", quietly = TRUE)) {
  install.packages("peakRAM", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("fs", quietly = TRUE)) {
  install.packages("fs", repos = "https://cloud.r-project.org")
}
library(peakRAM)
library(fs)

# Load Hal wavelet estimation module
WSE_Path = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/src/WaveletShrinkageEstimation.R")
print("Load Hal wavelet estimation module")
source(WSE_Path)

WT_Path = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/src/WaveletTransform.R")
print("Load Hal wavelet Transformation module")
source(WT_Path)

Method_Path = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/SimMethodPathSource_r.R")
source(Method_Path)

BS_Path = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/BSSim_r.R")
source(BS_Path)

#DS
DS=4

#最大解像度（ループ用） DS1.DS2:5 DS3.DS4:6
roopJ=6

# Level Dependent Threshold

#H-ldt-----------------------------------------------------------------
 directory_path = paste0("./BS-output/NDT_WSE/DS",DS,"-")
 out_path = paste0("./output/NDT_WSE/DS",DS,"-")
 BSSim(DS,roopJ,directory_path,out_path,"H","h","ldt")
 BSSim(DS,roopJ,directory_path,out_path,"H","s","ldt")

#TI-ldt----------------------------------------------------------------
# directory_path = paste0("./BS-output/NDT_WSE/TI/DS",DS,"-")
# out_path = paste0("./output/NDT_WSE/TI/DS",DS,"-")
# BSSim(DS,roopJ,directory_path,out_path,"TI","h","ldt")
# BSSim(DS,roopJ,directory_path,out_path,"TI","s","ldt")

# Universal Threshold -------------------------------------------------

# #A1-ut-----------------------------------------------------------------
#  directory_path = paste0("./BS-output/DT_Ans_WSE/A1/DS",DS,"-")
#  out_path = paste0("./output/DT_Ans_WSE/A1/DS",DS,"-")
#  BSSim(DS,roopJ,directory_path,out_path,"A1","h","ut")
#  BSSim(DS,roopJ,directory_path,out_path,"A1","s","ut")
#  
# #A2-ut-----------------------------------------------------------------
#  directory_path = paste0("./BS-output/DT_Ans_WSE/A2/DS",DS,"-")
#  out_path = paste0("./output/DT_Ans_WSE/A2/DS",DS,"-")
#  BSSim(DS,roopJ,directory_path,out_path,"A2","h","ut")
#  BSSim(DS,roopJ,directory_path,out_path,"A2","s","ut")
# 
# #A3-ut-----------------------------------------------------------------
  directory_path = paste0("./BS-output/DT_Ans_WSE/A3/DS",DS,"-")
  out_path = paste0("./output/DT_Ans_WSE/A3/DS",DS,"-")
  BSSim(DS,roopJ,directory_path,out_path,"A3","h","ut")
  BSSim(DS,roopJ,directory_path,out_path,"A3","s","ut")
# 
# #B1-ut-----------------------------------------------------------------
#  directory_path = paste0("./BS-output/DT_Bar_WSE/B1/DS",DS,"-")
#  out_path = paste0("./output/DT_Bar_WSE/B1/DS",DS,"-")
 # BSSim(DS,roopJ,directory_path,out_path,"B1","h","ut")
#  BSSim(DS,roopJ,directory_path,out_path,"B1","s","ut")
# 
# #B2-ut-----------------------------------------------------------------
  directory_path = paste0("./BS-output/DT_Bar_WSE/B2/DS",DS,"-")
  out_path = paste0("./output/DT_Bar_WSE/B2/DS",DS,"-")
  BSSim(DS,roopJ,directory_path,out_path,"B2","h","ut")
  BSSim(DS,roopJ,directory_path,out_path,"B2","s","ut")
# 
# #Fi-ut-----------------------------------------------------------------
#  directory_path = paste0("./BS-output/DT_Fit_WSE/DS",DS,"-")
#  out_path = paste0("./output/DT_Fit_WSE/DS",DS,"-")
#  BSSim(DS,roopJ,directory_path,out_path,"Fi","h","ut")
#  BSSim(DS,roopJ,directory_path,out_path,"Fi","s","ut")

#Fr-ut-----------------------------------------------------------------
# directory_path = paste0("./BS-output/DT_Fre_WSE/DS",DS,"-")
# out_path = paste0("./output/DT_Fre_WSE/DS",DS,"-")
# BSSim(DS,roopJ,directory_path,out_path,"Fr","h","ut")
# BSSim(DS,roopJ,directory_path,out_path,"Fr","s","ut")

# Leave-out Half cross validation Threshold ---------------------------

#A1-lht----------------------------------------------------------------
#  directory_path = paste0("./BS-output/DT_Ans_WSE/A1/DS",DS,"-")
#  out_path = paste0("./output/DT_Ans_WSE/A1/DS",DS,"-")
#  BSSim(DS,roopJ,directory_path,out_path,"A1","h","lht")
#  BSSim(DS,roopJ,directory_path,out_path,"A1","s","lht")
# 
# #A2-lht----------------------------------------------------------------
#  directory_path = paste0("./BS-output/DT_Ans_WSE/A2/DS",DS,"-")
#  out_path = paste0("./output/DT_Ans_WSE/A2/DS",DS,"-")
#  BSSim(DS,roopJ,directory_path,out_path,"A2","h","lht")
#  BSSim(DS,roopJ,directory_path,out_path,"A2","s","lht")
# 
# #A3-lht----------------------------------------------------------------
#  directory_path = paste0("./BS-output/DT_Ans_WSE/A3/DS",DS,"-")
#  out_path = paste0("./output/DT_Ans_WSE/A3/DS",DS,"-")
#  BSSim(DS,roopJ,directory_path,out_path,"A3","h","lht")
#  BSSim(DS,roopJ,directory_path,out_path,"A3","s","lht")

#B1-lht----------------------------------------------------------------
# directory_path = paste0("./BS-output/DT_Bar_WSE/B1/DS",DS,"-")
# out_path = paste0("./output/DT_Bar_WSE/B1/DS",DS,"-")
# BSSim(DS,roopJ,directory_path,out_path,"B1","h","lht")
# BSSim(DS,roopJ,directory_path,out_path,"B1","s","lht")

#B2-lht----------------------------------------------------------------
# directory_path = paste0("./BS-output/DT_Bar_WSE/B2/DS",DS,"-")
# out_path = paste0("./output/DT_Bar_WSE/B2/DS",DS,"-")
# BSSim(DS,roopJ,directory_path,out_path,"B2","h","lht")
# BSSim(DS,roopJ,directory_path,out_path,"B2","s","lht")

#Fi-lht----------------------------------------------------------------
# directory_path = paste0("./BS-output/DT_Fit_WSE/DS",DS,"-")
# out_path = paste0("./output/DT_Fit_WSE/DS",DS,"-")
# BSSim(DS,roopJ,directory_path,out_path,"Fi","h","lht")
# BSSim(DS,roopJ,directory_path,out_path,"Fi","s","lht")

#Fr-lht----------------------------------------------------------------
# directory_path = paste0("./BS-output/DT_Fre_WSE/DS",DS,"-")
# out_path = paste0("./output/DT_Fre_WSE/DS",DS,"-")
# BSSim(DS,roopJ,directory_path,out_path,"Fr","h","lht")
# BSSim(DS,roopJ,directory_path,out_path,"Fr","s","lht")

# Leave-out half Universal Threshold ----------------------------------

# #A1-lut----------------------------------------------------------------
#  directory_path = paste0("./BS-output/DT_Ans_WSE/A1/DS",DS,"-")
#  out_path = paste0("./output/DT_Ans_WSE/A1/DS",DS,"-")
#  BSSim(DS,roopJ,directory_path,out_path,"A1","h","lut")
#  BSSim(DS,roopJ,directory_path,out_path,"A1","s","lut")
# 
# #A2-lut----------------------------------------------------------------
#  directory_path = paste0("./BS-output/DT_Ans_WSE/A2/DS",DS,"-")
#  out_path = paste0("./output/DT_Ans_WSE/A2/DS",DS,"-")
#  BSSim(DS,roopJ,directory_path,out_path,"A2","h","lut")
#  BSSim(DS,roopJ,directory_path,out_path,"A2","s","lut")
# 
# #A3-lut----------------------------------------------------------------
#  directory_path = paste0("./BS-output/DT_Ans_WSE/A3/DS",DS,"-")
#  out_path = paste0("./output/DT_Ans_WSE/A3/DS",DS,"-")
#  BSSim(DS,roopJ,directory_path,out_path,"A3","h","lut")
#  BSSim(DS,roopJ,directory_path,out_path,"A3","s","lut")

#B1-lut----------------------------------------------------------------
# directory_path = paste0("./BS-output/DT_Bar_WSE/B1/DS",DS,"-")
# out_path = paste0("./output/DT_Bar_WSE/B1/DS",DS,"-")
# BSSim(DS,roopJ,directory_path,out_path,"B1","h","lut")
# BSSim(DS,roopJ,directory_path,out_path,"B1","s","lut")

#B2-lut----------------------------------------------------------------
# directory_path = paste0("./BS-output/DT_Bar_WSE/B2/DS",DS,"-")
# out_path = paste0("./output/DT_Bar_WSE/B2/DS",DS,"-")
# BSSim(DS,roopJ,directory_path,out_path,"B2","h","lut")
# BSSim(DS,roopJ,directory_path,out_path,"B2","s","lut")

#Fi-lut----------------------------------------------------------------
# directory_path = paste0("./BS-output/DT_Fit_WSE/DS",DS,"-")
# out_path = paste0("./output/DT_Fit_WSE/DS",DS,"-")
# BSSim(DS,roopJ,directory_path,out_path,"Fi","h","lut")
# BSSim(DS,roopJ,directory_path,out_path,"Fi","s","lut")

#Fr-lut----------------------------------------------------------------
# directory_path = paste0("./BS-output/DT_Fre_WSE/DS",DS,"-")
# out_path = paste0("./output/DT_Fre_WSE/DS",DS,"-")
# BSSim(DS,roopJ,directory_path,out_path,"Fr","h","lut")
# BSSim(DS,roopJ,directory_path,out_path,"Fr","s","lut")
