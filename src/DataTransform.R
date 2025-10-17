# -----------------------------------------------
# Anscombe transformation
# -----------------------------------------------
AnscombeTransformFromGroups = function(Groups,Var)
{
    GroupsLength = length(Groups)
    Lists = list()
    i = 1
    while(i <= GroupsLength)
    {
        Lists = append(Lists, list(AnscombeTransformFromGroup(Groups[[i]],Var)))
        i = i + 1
    }
    return(Lists)
}

#Applying the Anscombe transformation to a data set with a variance of Var after transformation
AnscombeTransformFromGroup = function(Groups,Var)
{
    Anscombelist = c()
    GroupLength = length(Groups)
    i = 1
    while(i <= GroupLength)
    {
        a = Groups[[i]] + 3/8
        b = a**0.5
        c = b * 2 * (Var**0.5)
        Anscombelist = append(Anscombelist, c)
        i = i + 1
    }
    return(Anscombelist)
}

#The inverse Anscombe transformation is applied to multiple data sets simultaneously, and the variance before transformation is Var
InverseAnscombeTransformFromGroups = function(AnscombeData,Var)
{
        #copy.deepcopy() は Pythonのコードで、Rでは存在しません 代わりに deepcopy の処理は不要か、明示的に AnscombeData <- lapply(AnscombeData, identity) などにすべき。
    #AnscombeData = copy.deepcopy(AnscombeData)
    AnscombeData = lapply(AnscombeData, identity)
    GroupsLength = length(AnscombeData)
    i = 1
    Lists = list()
    while(i <= GroupsLength)
    {
        Lists = append(Lists, list(InverseAnscombeTransformFromGroup(AnscombeData[[i]],Var)))
        i = i + 1
    }
    return(Lists)
}

#Applying the inverse Anscombe transformation to a dataset with a variance of Var before transformation
InverseAnscombeTransformFromGroup = function(AnscombeData,Var)
{
    GroupsLength = length(AnscombeData)
    i = 1
    Lists = c()
    while(i <= GroupsLength)
    {
        a = AnscombeData[[i]]
        b = a * a
        d = (2 * (Var**0.5)) ** -2
        c = d*b - 3/8
        c = round(c, 11)
        Lists = append(Lists, c)
        i = i + 1
    }
    return(Lists)
}

# -----------------------------------------
# The inverse Anscombe transformation 2
# ((si/2)^2)-1/8
# -----------------------------------------
InverseAnscombeTransform2FromGroups = function(AnscombeData,Var)
{
    GroupsLength = length(AnscombeData)
    i = 1
    Lists = list()
    while(i <= GroupsLength)
    {
            Lists = append(Lists, list(InverseAnscombeTransform2FromGroup(AnscombeData[[i]],Var)))
            i = i + 1
    }
    return(Lists)
}

#Applying the inverse Anscombe transformation 2 to a dataset with a variance of Var before transformation
InverseAnscombeTransform2FromGroup = function(AnscombeData,Var)
{
    GroupsLength = length(AnscombeData)
    i = 1
    Lists = c()
    while(i <= GroupsLength)
    {
            a = AnscombeData[[i]]
            b = a * a
            d = (2 * (Var**0.5)) ** -2
            c = d*b - 1/8
            c = round(c, 11)
            Lists = append(Lists, c)
            i = i + 1
    }
    return(Lists)
}


# -----------------------------------------
# The inverse Anscombe transformation 3
# (si^2)/4+sqrt(3/2)/(4*si)-11/(8*(si^2))+5*sqrt(3/2)/(8*(si^3))-1/8
# -----------------------------------------
InverseAnscombeTransform3FromGroups = function(AnscombeData,Var)
{
    GroupsLength = length(AnscombeData)
    i = 1
    Lists = list()
    while(i <= GroupsLength)
    {
            Lists = append(Lists, list(InverseAnscombeTransform3FromGroup(AnscombeData[[i]],Var)))
            i = i + 1
    }
    return(Lists)
}


#Applying the inverse Anscombe transformation 3 to a dataset with a variance of Var before transformation
InverseAnscombeTransform3FromGroup = function(AnscombeDataList,Var)
{
    AnscombeData = unlist(AnscombeDataList)
    GroupsLength = length(AnscombeData)
    i = 1
    Lists = c()
    while(i <= GroupsLength)
    {
            a = AnscombeData[[i]]
            b = a * a
            d = (2 * (Var**0.5)) ** -2
            e = a**(-1)
            f = a**(-2)
            g = a**(-3)
            c = d*b + (d**-0.5)*((3/2)**(0.5))*e - (d**-1)*11*f/2 + (d**-1.5)*5*((3/2)**(0.5))*g/4- 1/8
            if(a < 2*(3/8**(0.5)))
            {
                c = 0
            }
            c = round(c, 11)
            Lists = append(Lists, c)
            i = i + 1
    }
    return(Lists)
}

# -----------------------------------------------
# Bartlet
# -----------------------------------------------
BartlettTransformFromGroups = function(Groups,Var)
{
    GroupsLength = length(Groups)
    Lists = list()
    i = 1
    while(i <= GroupsLength)
    {
        Lists = append(Lists, list(BartlettTransformFromGroup(Groups[[i]],Var)))
        i = i + 1
    }
    return(Lists)
}

#Applying a Bartlett transformation to a data set with a variance of Var after transformation
BartlettTransformFromGroup = function(Groups,Var)
{
    Lists = c()
    GroupLength = length(Groups)
    i = 1
    while(i <= GroupLength)
    {
        a = Groups[[i]] + 0.5
        b = a**0.5
        c = b * 2 * (Var**0.5)
        Lists = append(Lists, c)
        i = i + 1
    }
    return(Lists)
}

#The inverse Bartlett transformation is applied simultaneously to multiple data sets, and the variance before transformation is Var
InverseBartlettTransformFromGroups = function(Groups,Var)
{
    GroupsLength = length(Groups)
    Lists = list()
    i = 1
    while(i <= GroupsLength)
    {
        Lists = append(Lists, list(InverseBartlettTransformFromGroup(Groups[[i]],Var)))
        i = i + 1
    }
    return(Lists)
}

#Applying an inverse Bartlett transformation to a dataset with a pre-transformation variance of Var
InverseBartlettTransformFromGroup = function(BartlettData,Var)
{
    GroupsLength = length(BartlettData)
    i = 1
    Lists = c()
    while(i <= GroupsLength)
    {
        a = BartlettData[[i]] * BartlettData[[i]]
        b = (2 * (Var**0.5)) ** -2
        c = b*a - 0.5
        c = round(c, 11)
        Lists = append(Lists, c)
        i = i + 1
    }
    return(Lists)
}

# -----------------------------------------------
# Applying Bartlett transformation 2
# bi=2*sqrt(yi)
# -----------------------------------------------
BartlettTransform2FromGroups = function(Groups,Var)
{
    GroupsLength = length(Groups)
    Lists = list()
    i = 1
    while(i <= GroupsLength)
    {
            Lists = append(Lists, list(BartlettTransform2FromGroup(Groups[[i]],Var)))
            i = i + 1
    }
    return(Lists)
}


#Applying a Bartlett transformation 2 to a data set with a variance of Var after transformation
BartlettTransform2FromGroup = function(Groups,Var)
{
    Lists = c()
    GroupLength = length(Groups)
    i = 1
    while(i <= GroupLength)
    {
            a = Groups[[i]]
            b = a**0.5
            c = b * 2 * (Var**0.5)
            Lists = append(Lists, c)
            i = i + 1
    }
    return(Lists)
}


# -----------------------------------------------
# The inverse Bartlett transformation 2
# (bi^2)/4
# -----------------------------------------------
InverseBartlettTransform2FromGroups = function(Groups,Var)
{
    GroupsLength = length(Groups)
    Lists = list()
    i = 1
    while(i <= GroupsLength)
    {
            Lists = append(Lists, list(InverseBartlettTransform2FromGroup(Groups[[i]],Var)))
            i = i + 1
    }
    return(Lists)
}


#Applying an inverse Bartlett transformation 2 to a dataset with a pre-transformation variance of Var
InverseBartlettTransform2FromGroup = function(BartlettData,Var)
{
    GroupsLength = length(BartlettData)
    i = 1
    Lists = c()
    while(i <= GroupsLength)
    {
            a = BartlettData[[i]] * BartlettData[[i]]
            b = (2 * (Var**0.5)) ** -2
            c = b*a
            c = round(c, 11)
            Lists = append(Lists, c)
            i = i + 1
    }
    return(Lists)
}


# -----------------------------------------------
# Fisz
# -----------------------------------------------
FiszTransformFromGroups = function(ScalingCoefficients,WaveletCoefficients,Var)
{
    GroupsLength = length(ScalingCoefficients)
    Lists = list()
    i = 1
    while(i <= GroupsLength)
    {
        Lists = append(Lists, list(FiszTransformFromGroup(ScalingCoefficients[[i]],WaveletCoefficients[[i]],Var)))
        i = i + 1
    }
    return(Lists)
}

#Applying the Fisz transformation to a data set, the variance after transformation is Var
FiszTransformFromGroup = function(ScalingcCoefficient,waveletCoe,Var)
{
    Lists = list()
    GroupLength = length(ScalingcCoefficient)
    j = 1
    while(j <= GroupLength)
    {
        i = 1
        levelLength = length(ScalingcCoefficient[[j]])
        coeList = c()
        while(i <= levelLength)
        {
            if(ScalingcCoefficient[[j]][i] == 0)
            {
                coeList = append(coeList, 0.0)
            }
            else
            {
                if(ScalingcCoefficient[[j]][i] < 0)
                {
                    print("FiszTransformFromGroup")
                }
                coeList = append(coeList, (Var**0.5) * waveletCoe[[j]][[i]]/(ScalingcCoefficient[[j]][i]**0.5))
            }
            i = i + 1
        }
        Lists = append(Lists, list(coeList))
        j = j + 1
    }
    return(Lists)
}

#The inverse Fisz transformation is applied simultaneously to multiple data sets, and the variance before transformation is Var
InverseFiszTransformFromGroups = function(ScalingCoefficients,FiszCoes,Var)
{
    GroupsLength = length(ScalingCoefficients)
    CsList = list()
    DsList = list()
    Lists = list()
    i = 1
    while(i <= GroupsLength)
    {
        a = InverseFiszTransformFromGroup(ScalingCoefficients[[i]],FiszCoes[[i]],Var)
        CsList = append(CsList, list(a[[1]]))
        DsList = append(DsList, list(a[[2]]))
        i = i + 1
    }

    Lists = append(Lists, list(CsList))
    Lists = append(Lists, list(DsList))
    return(Lists)
}

#Apply the Fisz transformation to a data set with a variance of Var before transformation
InverseFiszTransformFromGroup = function(ScalingcCoefficient,FiszCoefficient,Var)
{
    Lists = list()
    CsLists = list()
    DsLists = list()
    GroupLength = length(ScalingcCoefficient)
    j = GroupLength
    while(j > 0)
    {
        levelLength = length(ScalingcCoefficient[[j]])
        DsList = c()
        i = 1
        while(i <= levelLength)
        {
            DsList = append(DsList, FiszGetDs(ScalingcCoefficient[[j]][i],FiszCoefficient[[j]][i],Var))
            i = i + 1
        }
        i = 1
        while(i <= levelLength && j > 1)
        {
            ScalingcCoefficient[[j - 1]][2 * i - 1]     = ScalingcCoefficient[[j]][i] + DsList[i];
            ScalingcCoefficient[[j - 1]][2 * i] = ScalingcCoefficient[[j]][i] - DsList[i];
            if(ScalingcCoefficient[[j - 1]][2 * i - 1] < 0)
            {
                ScalingcCoefficient[[j - 1]][2 * i - 1] = 0
            }
            if(ScalingcCoefficient[[j - 1]][2 * i] < 0)
            {
                ScalingcCoefficient[[j - 1]][2 * i] = 0
            }
            i = i + 1;
        }
        DsLists = append(DsLists, list(DsList))
        j = j - 1
    }
    DsListx = list()
    i = GroupLength;
    while(i >= 1)
    {
        DsListx = append(DsListx, list(DsLists[[i]]))
        i = i - 1;
    }
    Lists = append(Lists, list(ScalingcCoefficient))
    Lists = append(Lists, list(DsListx))
    return( Lists)
}
#Wavelet coefficients in Poisson space are calculated from scale coefficients and wavelet coefficients in Gaussian space
FiszGetDs = function(c,f,Var)
{
    f = f/(Var ** 0.5)
    res = f*(c**0.5)
    res = round(res, 11)
    return(res)
}

# -----------------------------------------------
# Freeman
# -----------------------------------------------
FreemanTransformFromGroups = function(Groups,Var)
{
    GroupsLength = length(Groups)
    Lists = list()
    i = 1
    while(i <= GroupsLength)
    {
            Lists = append(Lists, list(FreemanTransformFromGroup(Groups[[i]],Var)))
            i = i + 1
    }
    return(Lists)
}

#Applying a Freeman transformation to a data set with a variance of Var after transformation
FreemanTransformFromGroup = function(Groups,Var)
{
    Lists = c()
    GroupLength = length(Groups)
    i = 1
    while(i <= GroupLength)
    {
            a = Groups[[i]] + 1
            b = a**0.5
            d = Groups[[i]]
            e = d**0.5
            c = b * (Var**0.5) + e * (Var**0.5)
            Lists = append(Lists, c)
            i = i + 1
    }
    return(Lists)
}

#The inverse Freeman transformation is applied simultaneously to multiple data sets, and the variance before transformation is Var
InverseFreemanTransformFromGroups = function(Groups,Var)
{
    GroupsLength = length(Groups)
    Lists = list()
    i = 1
    while(i <= GroupsLength)
    {
            Lists = append(Lists, list(InverseFreemanTransformFromGroup(Groups[[i]],Var)))
            i = i + 1
    }
    return(Lists)
}

#Applying an inverse Freeman transformation to a dataset with a pre-transformation variance of Var
InverseFreemanTransformFromGroup = function(FrTransformData,Var)
{
    GroupsLength = length(FrTransformData)
    i = 1
    Lists = c()
    while(i <= GroupsLength)
    {
            a = FrTransformData[[i]] * FrTransformData[[i]]
            b = (2 * (Var**0.5)) ** -2
            d = a**(-1)
            e = b**(-1)
            c = b*a +e*d - 0.5
            c = round(c, 11)
            Lists = append(Lists, c)
            i = i + 1
    }
    return(Lists)
}