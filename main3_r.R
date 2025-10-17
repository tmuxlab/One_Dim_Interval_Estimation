#UTF-8
library(e1071)
library(tidyverse)

# Load wavelet conversion module
WaveletTransform_Path = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/src/WaveletShrinkageEstimation.R")
print("Load wavelet conversion module")
source(WaveletTransform_Path)

# Load PE.IE
BSintervals_Path = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/BS_intervals_r.R")
print("Load PE+IE module")
source(BSintervals_Path)

DS = 4
roopJ = 6
ConfidenceAll <- function(DS,dir_path,DTName,ThreName,Thre,i){
  #BSintervals(DS,dir_path,"COS",DTName,ThreName,Thre,i)
  #BSintervals(DS,dir_path,"GS-",DTName,ThreName,Thre,i)
  ##BSintervals(DS,dir_path,"GS+",DTName,ThreName,Thre,i)
  #BSintervals(DS,dir_path,"OrdSta",DTName,ThreName,Thre,i)
  #BSintervals(DS,dir_path,"ThiArr",DTName,ThreName,Thre,i)
  BSintervals(DS,dir_path,"ThiInt",DTName,ThreName,Thre,i)
  # BSintervals(DS,dir_path,"TimTra",DTName,ThreName,Thre,i)
  # BSintervals(DS,dir_path,"IntGen",DTName,ThreName,Thre,i)
  # BSintervals(DS,dir_path,"RGS",DTName,ThreName,Thre,i)
}

# s ThresholdMode
for(i in 2:roopJ){
  ConfidenceAll(DS,"NDT_WSE","H","h","ldt",i)
  ConfidenceAll(DS,"DT_Ans_WSE/A3","A3","h","ut",i)
  ConfidenceAll(DS,"DT_Bar_WSE/B2","B2","h","ut",i)
  # ConfidenceAll(DS,"DT_Ans_WSE/A1","A1","s","ut",i)
  # ConfidenceAll(DS,"DT_Ans_WSE/A2","A2","s","ut",i)
  # ConfidenceAll(DS,"DT_Ans_WSE/A3","A3","s","ut",i)
  # ConfidenceAll(DS,"DT_Bar_WSE/B1","B1","s","ut",i)
  # ConfidenceAll(DS,"DT_Bar_WSE/B2","B2","s","ut",i)
  # ConfidenceAll(DS,"DT_Fit_WSE","Fi","s","ut",i)  
  # ConfidenceAll(DS,"DT_Ans_WSE/A1","A1","s","lut",i)
  # ConfidenceAll(DS,"DT_Ans_WSE/A2","A2","s","lut",i)
  # ConfidenceAll(DS,"DT_Ans_WSE/A3","A3","s","lut",i)
  # ConfidenceAll(DS,"DT_Ans_WSE/A1","A1","s","lht",i)
  # ConfidenceAll(DS,"DT_Ans_WSE/A2","A2","s","lht",i)
  # ConfidenceAll(DS,"DT_Ans_WSE/A3","A3","s","lht",i)
}

# h A2 GS- lht i=2 まで保存した
# h ThresholdMode
#for(i in 3:roopJ){
  #ConfidenceAll(DS,"DT_Ans_WSE/A1","A1","h","ut",i)
  #ConfidenceAll(DS,"DT_Ans_WSE/A2","A2","h","ut",i)
  #ConfidenceAll(DS,"DT_Ans_WSE/A3","A3","h","ut",i)
  # ConfidenceAll(DS,"DT_Bar_WSE/B1","B1","h","ut",i)
  #ConfidenceAll(DS,"DT_Bar_WSE/B2","B2","h","ut",i)
  # ConfidenceAll(DS,"DT_Fit_WSE","Fi","h","ut",i)  
  #ConfidenceAll(DS,"DT_Ans_WSE/A1","A1","h","lut",i)
  # ConfidenceAll(DS,"DT_Ans_WSE/A2","A2","h","lut",i)
  # ConfidenceAll(DS,"DT_Ans_WSE/A3","A3","h","lut",i)
  #ConfidenceAll(DS,"DT_Ans_WSE/A1","A1","h","lht",i)
  # ConfidenceAll(DS,"DT_Ans_WSE/A2","A2","h","lht",i)
  # ConfidenceAll(DS,"DT_Ans_WSE/A3","A3","h","lht",i)
#}
