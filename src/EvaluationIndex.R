# MSE(Mean Squared Error) 実際はRMSEでは？
MSE = function(originData, estimationList)
{
  if(length(originData)!= length(estimationList))
  {
    return(FALSE)
  }
  index = 1
  total = 0
  while(index <= length(originData))
  {
    total = total + (originData[[index]] - estimationList[[index]])**2
    index = index + 1
  }
  total = sqrt(total)
  return(total / length(originData))
}