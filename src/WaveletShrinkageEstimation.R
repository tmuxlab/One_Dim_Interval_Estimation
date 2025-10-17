# Load wavelet conversion module
WaveletTransformPath = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/src/WaveletTransform.R")
source(WaveletTransformPath)
# Load data conversion module
DtPath = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/src/DataTransform.R")
source(DtPath)
# Load Threshold Module
ThresholdPath = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/src/Threshold.R")
source(ThresholdPath)
# Load Evaluation Module
EvaluationPath = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/src/EvaluationIndex.R")
source(EvaluationPath)

# Hal wavelet estimation without data transformation
Wse = function(Data, DataTransform, ThresholdName, ThresholdMode, Var, Index, InitThresholdValue)
{
  if(DataTransform == "none" && ThresholdName != "ldt"){
    print("Please check the parameter. If you want to use DataTransform=none, please set ThresholdName=ldt.")
  }
  else if(DataTransform != "none" && ThresholdName == "ldt"){
    print("Please check the parameter. If you want to use ThresholdName=ldt, please set DataTransform=none.")
  }
  else{
    GroupLength = 2^Index
    # Get data length
    DataLength = length(Data)
    if(GroupLength >= GetGroupLength(DataLength)){
      # Get subdata length
      GroupLength = GetGroupLength(DataLength)
    }

    # Cut the original data into a number of sub-data of length 2^J
    Groups = GetGroups(Data,GroupLength)

    if(DataTransform == "Fi"){
      # Transform the sub-data into Gaussian data by Fisz transformation
      Cs1  = GetScalingCoefficientsFromGroups(Groups)
      Ds1  = GetWaveletCoefficientsFromGroups(Cs1)
      Fi1  = FiszTransformFromGroups(Cs1,Ds1,Var)
      FiszGroups = InverseHaarWaveletTransformForGroups(Cs1,Fi1)

      FiszGroups = lapply(FiszGroups, function(x) x/GroupLength**0.5)
      
      # Calculate c
      Cs2  = GetScalingCoefficientsFromGroups(FiszGroups)
      #Calculate d
      Ds2  = GetWaveletCoefficientsFromGroups(Cs2)
      
      # Noise reduction of wavelet coefficients using ThresholdMode noise reduction rule, ThresholdName threshold
      DenoiseDs2 = ThresholdForGroups(Ds2,ThresholdMode,ThresholdName, DataTransform, Groups, InitThresholdValue)
        
      # Perform inverse Fisz data conversion
      InverseGropus = InverseHaarWaveletTransformForGroups(Cs2,DenoiseDs2)
      Cs3  = GetScalingCoefficientsFromGroups(InverseGropus)
      Fs2  = GetWaveletCoefficientsFromGroups(Cs3)
      CDs  = InverseFiszTransformFromGroups(Cs3,Fs2,Var)
      Cs4 = CDs[[1]]
      Ds3 = CDs[[2]]
      
      # Perform inverse wavelet conversion
      ThresholdedGroups = InverseHaarWaveletTransformForGroups(Cs4,Ds3)
      ThresholdedGroups = lapply(ThresholdedGroups, function(x) x*GroupLength**0.5)
      
      # Perform moving average
      ThresholdedData= MovingAverage(ThresholdedGroups,DataLength)

      Result = list(EstimationData=ThresholdedData, Cs=Cs4,Ds=Ds3, denoiseDs=DenoiseDs2)
    }
    else{
      if(DataTransform == "A1" || DataTransform == "A2"|| DataTransform == "A3"){
        #Transform sub-data to Gaussian data by Anscombe
        Groups = AnscombeTransformFromGroups(Groups,Var)
      }
      else if(DataTransform == "B1"){
        #Transform sub-data to Gaussian data by Bartlet
        Groups = BartlettTransformFromGroups(Groups,Var)
      }
      else if(DataTransform == "B2"){
        #Transform sub-data to Gaussian data by Bartlet
        Groups = BartlettTransform2FromGroups(Groups,Var)
      }
      else if (DataTransform == "Fr") {
        Groups = FreemanTransformFromGroups(Groups,Var)
      }
      else{
        Groups = Groups
      }
      Groups = lapply(Groups, function(x) x/(GroupLength**0.5))

      # Calculate c
      Cs = GetScalingCoefficientsFromGroups(Groups)
      #Calculate d
      Ds = GetWaveletCoefficientsFromGroups(Cs)

      DenoisedDs = ThresholdForGroups(Ds,ThresholdMode,ThresholdName, DataTransform, Groups, InitThresholdValue)
      # Perform inverse wavelet conversion
      ThresholdedGroups = InverseHaarWaveletTransformForGroups(Cs,DenoisedDs)
      ThresholdedGroups = lapply(ThresholdedGroups, function(x) x*GroupLength**0.5)

      # Perform moving average
      if(ThresholdName == "none"){
        ThresholdedData = ThresholdedGroups
      }
      else {
        ThresholdedData = MovingAverage(ThresholdedGroups,DataLength)
      }

      if(DataTransform == "A1"){
      # Perform inverse Anscombe data conversion
      ThresholdedData = InverseAnscombeTransformFromGroup(ThresholdedData,Var);
      }
      else if(DataTransform == "A2"){
        # Perform inverse Anscombe data conversion
        ThresholdedData = InverseAnscombeTransform2FromGroup(ThresholdedData,Var);
      }
      else if(DataTransform == "A3"){
        # Perform inverse Anscombe data conversion
        ThresholdedData = InverseAnscombeTransform3FromGroup(ThresholdedData,Var);
      }
      else if(DataTransform == "B1"){
        # Perform inverse Anscombe data conversion
        ThresholdedData = InverseBartlettTransformFromGroup(ThresholdedData,Var);
      }
      else if(DataTransform == "B2"){
        # Perform inverse Anscombe data conversion
        ThresholdedData = InverseBartlettTransform2FromGroup(ThresholdedData,Var);
      }
      else if (DataTransform == "Fr") {
        ThresholdedData = InverseFreemanTransformFromGroup(ThresholdedData,Var)
      }
      else{
        ThresholdedData = ThresholdedData
      }
      i = 1
      while(i < DataLength){
        if(ThresholdedData[i] <= 0){
          ThresholdedData[i] = 0
        }
        i = i + 1
      }
      Result = list(EstimationData=ThresholdedData, Cs=Cs,Ds=Ds, DenoisedDs=DenoisedDs)
    }

    return(Result)
  }
}

# Translation-invariant Hal wavelet estimation without data transformation
Tipsh = function(Data, ThresholdMode, Var, Index)
{
  ThresholdName = "ldt"
  GroupLength = 2^Index
  # Get data length
  DataLength = length(Data)
  if(GroupLength >= GetGroupLength(DataLength)){
          # Get subdata length
          GroupLength = GetGroupLength(DataLength)
  }
  # Cut the original data into a number of sub-data of length 2^J
  Groups = GetGroups(Data,GroupLength)
  Groups = lapply(Groups, function(x) x/GroupLength**0.5)
  ThresholdedGroups = list()
  #Transration-Invariant Denoising
  for(i in 1:(DataLength-GroupLength+1)){
          templist=Groups[[i]]
          shiftgroup=list()
          Lists=list()
          for(h in 1:(GroupLength-1)){
                  shiftgroup = list(c(templist[(h+1):GroupLength],templist[1:h]))
                  Cs = GetScalingCoefficientsFromGroup(as.numeric(shiftgroup[[1]]))
                  Ds = GetWaveletCoefficientsFromGroup(Cs)
                  DenoisedDs = ThresholdForGroup(Ds,ThresholdMode,ThresholdName)
                  ThresholdedGroup = InverseHaarWaveletTransformForGroup(Cs,DenoisedDs)
                  Lists = append(Lists,list(ThresholdedGroup))
          }
          Cs=GetScalingCoefficientsFromGroup(templist)
          DS=GetWaveletCoefficientsFromGroup(Cs)
          DenoisedDs = ThresholdForGroup(Ds,ThresholdMode,ThresholdName)
          ThresholdedGroup = InverseHaarWaveletTransformForGroup(Cs,DenoisedDs)
          Lists = append(Lists,list(ThresholdedGroup))
          for(h in 1:(GroupLength-1)){
                  templist=Lists[[h]]
                  Lists[h]=list(c(templist[(GroupLength-h+1):(GroupLength)],templist[1:(GroupLength-h)]))
          }
          SumList=c()
          for(j in 1:GroupLength){
                  sl=0
                  for(k in 1:GroupLength){
                          sl=sl+Lists[[k]][j]
                  }
                  SumList = c(SumList,sl)
          }
          ThresholdedGroups[[i]]=(SumList/GroupLength)
  }
  
  ThresholdedGroups = lapply(ThresholdedGroups, function(x) x*GroupLength**0.5)
  # Perform moving average
  EstimationData = MovingAverage(ThresholdedGroups, DataLength)
  EstimationData = pmax(EstimationData, 0)
  Result = list(EstimationData=EstimationData, Cs=Cs,Ds=Ds, DenoisedDs=DenoisedDs)
  # Return Results
  return(Result)
}

# Load data from file
LoadData = function(DataPath)
{
  DataPath = paste0(dirname(rstudioapi::getSourceEditorContext()$path),DataPath)
  Ds = read.table(DataPath)[2]
  Ds = as.numeric(Ds$V2)
  return(Ds)
}


# creating file format
CreateFile = function(i, ResultPath, time)
{
    FileNameEstimationData = paste0(time,"_EstimationData_J=",i  ,".csv")
    FileNameCoefficients = paste0(time,"_Coefficients_J=",i  ,".csv")
    FileNameVariable = paste0(time,"_Variable_J=",i  ,".RData")
    EstimationData = paste0(ResultPath, FileNameEstimationData)
    Coefficients = paste0(ResultPath, FileNameCoefficients)
    Variable = paste0(ResultPath, FileNameVariable)
    FilePath = list(EstimationData = EstimationData, Coefficients = Coefficients, Variable = Variable)
    return(FilePath)
}

# creating result
CreateResult = function(Hard, Soft, Index, ResultPath){
  time = Sys.time() %>% format("%H-%M-%S")
  EstimationData = list(Hard = round(Hard$EstimationData, digits = 3), Soft = round(Soft$EstimationData, digits = 3))
  hard_coe= rbind("Cs",as.data.frame(t(sapply(Hard$Cs, unlist))),"Ds",as.data.frame(t(sapply(Hard$Ds, unlist))),"DenoisedDs",as.data.frame(t(sapply(Hard$DenoisedDs, unlist))))
  soft_coe= rbind("Cs",as.data.frame(t(sapply(Soft$Cs, unlist))),"Ds",as.data.frame(t(sapply(Soft$Ds, unlist))),"DenoisedDs",as.data.frame(t(sapply(Soft$DenoisedDs, unlist))))
  Coefficients = rbind("Hard",hard_coe,"Soft",soft_coe)
  FilePath = CreateFile(Index, ResultPath, time)
  write.csv(EstimationData, FilePath$EstimationData, row.names = FALSE)
  write.csv(Coefficients, FilePath$Coefficients, row.names = FALSE)
  save(Hard, Soft, file = FilePath$Variable)
}

create_txt = function(i,j,directory_path,threshold,estimation){
  file_name_hdata = paste0("DS",j,"-k",estimation,"-h-",threshold,"(J=",i,").txt")
  file_name_sdata = paste0("DS",j,"-k",estimation,"-s-",threshold,"(J=",i,").txt")
  hdata = paste0(directory_path,file_name_hdata)
  sdata = paste0(directory_path,file_name_sdata)
  file_path = list(hdata=hdata,sdata=sdata)
  return(file_path)
}