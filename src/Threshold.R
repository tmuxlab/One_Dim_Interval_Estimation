ThresholdForGroups = function(Ds, ThresholdMode, ThresholdName, DataTransform, Groups , InitThresholdValue)
{
	GroupLength = length(Ds)
	Lists = list()
	u = 1
	while(u <= GroupLength)
	{
		TmpList = ThresholdForGroup(Ds[[u]],ThresholdMode,ThresholdName,DataTransform,Groups, InitThresholdValue, u)
		Lists = append(Lists, list(TmpList))
		u = u + 1
	}
	return(Lists)
}

#Apply the soft or hard thresholding method of ThresholdName to a set of wavelet coefficients
ThresholdForGroup = function(GroupWaveletCoefficients,ThresholdMode,ThresholdName,DataTransform,Groups, InitThresholdValue,j)
{
	if(ThresholdName == 'ut'|| ThresholdName == 'ldt' || ThresholdName == 'lut'|| ThresholdName == 'none'){
		DataLength = length(GroupWaveletCoefficients[[1]])
		t = 1000
		if(ThresholdName == 'ut')
		{
			t = GetUniversalThreshold(DataLength)
		}
		else if (ThresholdName == 'none') {
		   t = InitThresholdValue
		}
		Lists = list()
		Lists = append(Lists,list(GroupWaveletCoefficients[[1]]))
		i = 2
		GroupLength = length(GroupWaveletCoefficients)

		if(ThresholdName == 'ldt' || ThresholdName == 'lut')
		{
			C = GetScalingCoefficientsFromGroup(GroupWaveletCoefficients[[1]])
        	lam0 = mean(GroupWaveletCoefficients[[1]])*(DataLength**0.5)
		}

		while(i <= GroupLength)
		{
			if(ThresholdName == 'ldt')
			{
				TmpList = LdtThreshold(GroupWaveletCoefficients[[i]],ThresholdMode,i,DataLength,lam0)
			}
			else
			{
				if(ThresholdName == 'lut')
				{
					ut_dataLength = length(C[[i]])
					t = GetUniversalThreshold(ut_dataLength)
				}
				TmpList = ThresholdForOneLevel(GroupWaveletCoefficients[[i]],ThresholdMode,t)
			}
			Lists = append(Lists,list(TmpList))
			i = i + 1
		}
	}
	else{
			SubGroupLength = length(Groups[[1]])
			if(j != length(Groups))
			{
				NextValue = NextValue = Groups[[j+1]][SubGroupLength]
			}
			else
			{
				NextValue = NextValue = Groups[[1]][SubGroupLength]
			}
			t = LhtThreshold(Groups[[j]], DataTransform, ThresholdName, ThresholdMode, NextValue)
			j = j + 1
			Lists = list()
			Lists = append(Lists,list(GroupWaveletCoefficients[[1]]))
			i = 2
			GroupLength = length(GroupWaveletCoefficients)
			while(i <= GroupLength)
			{
				TmpList = ThresholdForOneLevel(GroupWaveletCoefficients[[i]],ThresholdMode,t)
				Lists = append(Lists,list(TmpList))
				i = i + 1
			}
	}
    return(Lists)
}

# ---------------------------------
# Get ut, lut threshold
# ut:Universal Threshold
# lut:Level-universal-Threshold
# lut is applied to the length of j-th
# level empirical wavelet coefficients.
# ---------------------------------
GetUniversalThreshold = function(GroupLength)
{
	a = log(GroupLength)
	b = 2*a
	c = b**0.5
	return(c)
}

# ---------------------------------
# Get ldt threshold
# ldt:Level-dependent-Threshold
# ---------------------------------
GetLevelDependentThreshold = function(J,NowLevel,Mean)
{
	a = 2 ** (-1 * 0.5 * (NowLevel+1))
	log2j = log(2 ** (J - NowLevel +1))
	b = 2 * log2j
	c = 4 * (log2j) ** 2
	d = 8 * Mean * log2j
	t = a * (b + ((c+d) ** 0.5))
	return(t)
}

# Thresholding the wavelet coefficients of a layer at a threshold value of t
ThresholdForOneLevel = function(WaveletCoefficients,ThresholdMode,t)
{
	coefficientsLength = length(WaveletCoefficients)
	TmpList = c()
	i = 1
	while(i <= coefficientsLength)
	{
		Tmp = Threshold(WaveletCoefficients[i],t,ThresholdMode)
		TmpList = append(TmpList,Tmp)
		i = i + 1
	}
	return(TmpList)
}

#Calculating ldt and thresholding the Data
LdtThreshold = function(Data,ThresholdMode,loop_level,DataLength,lam0)
{
	#Highest Resolution
	J = GetHighestResolutionLevel(DataLength)
	#Thresholding the Data one by one
	i = 1
	TmpList = c()
	Mean = lam0/length(Data)
	while (i <= length(Data))
	{
		#Get ldt threshold
		t = GetLevelDependentThreshold(J,loop_level,Mean)
		#Threshold processing
		DenoisedData = Threshold(Data[[i]],t,ThresholdMode)
		TmpList = append(TmpList,DenoisedData)
		i = i + 1
	}
	return(TmpList)
}

LhtThreshold = function(OriginalGroups, DataTransform, ThresholdName, ThresholdMode ,NextValue) {
	subgroup_len = length(OriginalGroups)
	Minimum = optim(par = 0, fn = LossFunction, OriginalGroup = OriginalGroups, DataTransform = DataTransform, ThresholdName = ThresholdName, ThresholdMode = ThresholdMode,NextValue = NextValue, method = "Brent",lower = -5,upper = 5)$par
	if (Minimum < 0) {
	Minimum = 0
	}
	ThresholdValue = ((1 - log(2) / log(subgroup_len)) ^ (-0.5)) * Minimum
	return(ThresholdValue)
}

LossFunction = function(t, OriginalGroup, DataTransform, ThresholdName, ThresholdMode, NextValue) {
  # Separate even-numbered and odd-numbered
  OddGroup = OriginalGroup[seq(1, length(OriginalGroup), by = 2)]
  EvenGroup = OriginalGroup[seq(2, length(OriginalGroup), by = 2)]

  OddIndex = log(length(OddGroup), base = 2)
  EvenIndex = log(length(EvenGroup), base = 2)

  # Perform WSE for even and odd numbers
 ThresholdedOddGroup = Wse(OddGroup, DataTransform, "none", ThresholdMode, 1, OddIndex, t)
 ThresholdedEvenGroup = Wse(EvenGroup, DataTransform, "none", ThresholdMode, 1, EvenIndex, t)

 OriginalGroup = append(OriginalGroup,NextValue)

	OddAverageList = list()
	EvenAverageList = list()


	for (i in 1:length(ThresholdedOddGroup$EstimationData)) {
	if (i != length(ThresholdedOddGroup$EstimationData)) {
		OddAverage = (ThresholdedOddGroup$EstimationData[i] + ThresholdedOddGroup$EstimationData[i + 1]) * 0.5
	} else {
		OddAverage = (ThresholdedOddGroup$EstimationData[i] + ThresholdedOddGroup$EstimationData[1]) * 0.5
	}
	OddAverageList= c(OddAverageList, OddAverage)
	}

	for (i in 1:length(ThresholdedEvenGroup$EstimationData)) {
	if (i != length(ThresholdedEvenGroup$EstimationData)) {
		EvenAverage = (ThresholdedEvenGroup$EstimationData[i] + ThresholdedEvenGroup$EstimationData[i + 1]) * 0.5
	} else {
		EvenAverage = (ThresholdedEvenGroup$EstimationData[i] + ThresholdedEvenGroup$EstimationData[1]) * 0.5
	}
	EvenAverageList= c(EvenAverageList, EvenAverage)
	}
  SquaredError = 0
  for (i in 1:length(ThresholdedOddGroup$EstimationData)) {
    OddSquaredError = (OddAverageList[[i]][1] - OriginalGroup[2 * i]) ^ 2
    EvenSquaredError = (EvenAverageList[[i]][1] - OriginalGroup[2 * i + 1]) ^ 2
    SquaredError = SquaredError + OddSquaredError + EvenSquaredError
  }
  return(SquaredError)
}

#Thresholding of the value Coefficient according to the threshold t
Threshold = function(Coefficient,t,ThresholdMode)
{
	if(ThresholdMode == 'h'){
		if(abs(Coefficient) <= t){
	    return(0)
		}
		else{
	    return(Coefficient)
		}
	}
	else{
		if(abs(Coefficient) <= t){
	    return(0)
		}
		else{
			if(Coefficient > 0){
	        	return(Coefficient - t)
	    	}
	    	else{
	        	return(Coefficient + t)
	    	}
		}
	}
}