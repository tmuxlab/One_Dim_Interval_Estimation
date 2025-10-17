#Get maximum resolution
GetHighestResolutionLevel = function(GroupLength)
{
  Level = log2(GroupLength)
  Level = as.integer(Level)
  return(Level)
}
#Get a value that satisfies the largest integer power of 2 less than or equal to dataLength.
#example：when DataLength = 62, return 32.
#reason：2^5 = 32 < 62 < 64 = 2^6
GetGroupLength = function(DataLength)
{
  i = 1
  x = i
  while (i <= DataLength)
  {
    x = i
    i = i * 2
  }
  return(x)
}

#Divide Data set Data of arbitrary length into multiple sub-datasets according to groupLength.
GetGroups = function(Data, GroupLength)
{
  i = 0
  DataLength = length(Data)
  TempList = list()
  while (GroupLength + i <= DataLength)
  {
    TempData = Data
    a = i+1
    b = GroupLength + i
    CutData = TempData[a : b]
    TempList = append(TempList,list(CutData))
    i = i + 1
  }
  return(TempList)
}

#Scale coefficients of discrete Hal wavelets expanded for a set of Data
GetScalingCoefficientsFromGroup = function(TimeList)
{
  Lists = list()
  J = GetHighestResolutionLevel( length(TimeList))
  Lists = append(Lists, list(TimeList))
  j = 1
  while(j <= J)
  {
    TempList = c()
    k = 1
    while(k <= 2**(J - j))
    {
      coe = (1/sqrt(2))*(Lists[[j]][2 * k - 1] + Lists[[j]][2 * k])
      TempList = append(TempList,coe)
      k = k + 1
    }
    Lists = append(Lists, list(TempList))
    j = j + 1
  }
  return(Lists)
}

#Wavelet coefficients of discrete Hal wavelets are simultaneously expanded for multiple Data sets
GetWaveletCoefficientsFromGroup = function(CoeList)
{
  Lists = list()
  J = GetHighestResolutionLevel(length(CoeList[[1]]) )
  Lists = append(Lists, list(CoeList[[1]]))
  j = 1
  while(j <= J)
  {
    TempList = c()
    k = 1
    while(k <= 2**(J - j))
    {
      c = (1/sqrt(2))*(CoeList[[j]][2 * k - 1] - CoeList[[j]][2 * k])
      TempList = append(TempList,c)
      k = k + 1
    }
    Lists = append(Lists,list(TempList))
    j = j + 1
  }
  return(Lists)
}

#Scale coefficients of the discrete Hal wavelet are simultaneously expanded for multiple Data sets
GetScalingCoefficientsFromGroups = function(Groups)
{
  Lists = list()
  GroupsLength = length(Groups)
  i = 1
  while(i <= GroupsLength)
  {
    TempList = GetScalingCoefficientsFromGroup(Groups[[i]])
    Lists = append(Lists,list(TempList))
    i = i + 1
  }
  return(Lists)
}

#Wavelet coefficients of discrete Hal wavelets unfolded for a set of Data
GetWaveletCoefficientsFromGroups = function(CS)
{
  Lists = list()
  GroupsLength = length(CS)
  i = 1
  while(i <= GroupsLength)
  {
    TempList = GetWaveletCoefficientsFromGroup(CS[[i]])
    Lists = append(Lists,list(TempList))
    i = i + 1
  }
  return(Lists)
}

#Convert a set of Haar wavelet coefficients with Haar scale coefficients to the original Data
InverseHaarWaveletTransformForGroup = function(ScalingCoefficient,WaveletCoefficient)
{
  GroupLength = length(ScalingCoefficient)
  if (GroupLength != length(WaveletCoefficient))
  {
    return(FALSE)
  }
  J = GroupLength
  k = 0
  j = GroupLength
  while(j > 1)
  {
    k = 1
    while(k <= 2**(J - j))
    {
      ScalingCoefficient[[j - 1]][2 * k - 1] = (1/sqrt(2))*(ScalingCoefficient[[j]][k] + WaveletCoefficient[[j]][k])
      ScalingCoefficient[[j - 1]][2 * k] = (1/sqrt(2))*(ScalingCoefficient[[j]][k] - WaveletCoefficient[[j]][k])
      k = k + 1
    }
    j = j - 1
  }
  return(ScalingCoefficient[[1]])
}

#Convert multiple sets of Haar wavelet coefficients with Haar scale coefficients to the original Data
InverseHaarWaveletTransformForGroups = function(Cs,Ds)
{
  GroupsLength = length(Cs)
  if (GroupsLength != length(Ds))
  {
    return(FALSE)
  }
  i = 1
  Lists = list()
  while(i <= GroupsLength)
  {
    TempList = InverseHaarWaveletTransformForGroup(Cs[[i]],Ds[[i]])
    Lists = append(Lists,list(TempList))
    i = i + 1
  }
  return(Lists)
}

#Average multiple sub-datasets by displacement and combine them into one set
MovingAverage = function(ThresholdedGroups,DataLength)
{
  dataSum = numeric(DataLength)
  counter = numeric(DataLength)
  Result  = numeric(DataLength)
  
  GroupsSum = length(ThresholdedGroups)
  GroupLength = length(ThresholdedGroups[[1]])
  i = 1
  while(i <= GroupsSum)
  {
    j = 1
    while(j <= GroupLength)
    {
      dataSum[i + j - 1] = dataSum[i + j - 1] + ThresholdedGroups[[i]][j]
      counter[i + j - 1] = counter[i + j - 1] + 1
      j = j + 1
    }
    i = i + 1
  }

  k = 1
  while(k <= DataLength)
  {
    Result[k] = dataSum[k] / counter[k]
    if(Result[k] < 0)
    {
     Result[k] = 0
    }
    k = k + 1
  }
  return(Result)
}
