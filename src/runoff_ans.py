from pcraster import *
from pcraster.framework import *
import math

# notes #
# out.map: location where discharge was measured
# observed.tss: observed discharge (same unit as discharge variable in script)

# helper functions to get values from a map
def getCellValue(Map, Row, Column):
  Value, Valid=cellvalue(Map, Row, Column)
  if Valid:
    return Value
  else:
    print 'missing value in input of getCellValue'

def getCellValueAtBooleanLocation(location,map):
  # map can be any type, return value always float
  valueMap=mapmaximum(ifthen(location,scalar(map)))
  value=getCellValue(valueMap,1,1)
  return value

#calculate the m, b, RSquared coefficient and the root mean squared error
def calculateRSquaredAndError(myModel):
  SE = 0
  SEymean = 0
  Error = 0
  
  ##Least-Squareds
  Xmean = myModel.sumX / nrOfTimeSteps
  Ymean = myModel.sumY / nrOfTimeSteps
  XYmean = myModel.sumXY / nrOfTimeSteps
  SqXmean = myModel.sumSqX / nrOfTimeSteps

  m = ((Xmean * Ymean) - (XYmean)) / (Xmean ** 2 - SqXmean) #slope
  b = Ymean - m * Xmean #intersection
  
  for j in range(0,181,1):
    SE = SE + (myModel.y[j] - (m * myModel.x[j] + b)) ** 2
    SEymean = SEymean + (myModel.y[j] - Ymean) ** 2
    Error = Error + ((myModel.x[j] - myModel.y[j]) ** 2)

  r2 = 1 - SE/SEymean
  Error = math.sqrt(Error / nrOfTimeSteps)

  return m, b, r2, Error

class MyFirstModel(DynamicModel):
  def __init__(self):
    DynamicModel.__init__(self)
    setclone('dem.map')

  def initial(self):
    # the measurement location
    self.observationLocation=self.readmap('out')

    dem = self.readmap('dem')
    elevationMeteoStation = 208.1
    elevationAboveMeteoStation = dem - elevationMeteoStation
    temperatureLapseRate = 0.005
    self.temperatureCorrection = elevationAboveMeteoStation * temperatureLapseRate
    self.report(self.temperatureCorrection,'tempCor')

    self.snow=0.0

    self.ldd=lddcreate(dem,1e31,1e31,1e31,1e31)
    self.report(self.ldd,'ldd')

    # example how to calculate total precipitation
    self.totPrecip=scalar(0)

    self.sqrError = 0
    self.sumX = 0
    self.sumSqX = 0
    self.sumY = 0
    self.sumXY = 0

    self.x = []
    self.y = []

  def dynamic(self):
    precipitation = timeinputscalar('precip.tss',1)
    self.report(precipitation,'pFromTss')
    temperatureObserved = timeinputscalar('temp.tss',1)
    self.report(temperatureObserved,'tempObs')
    temperature= temperatureObserved - self.temperatureCorrection
    self.report(temperature,'temp')

    freezing=temperature < FTemp
    self.report(freezing,'freez')
    snowFall=ifthenelse(freezing,precipitation,0.0)
    self.report(snowFall,'snowFall')
    rainFall=ifthenelse(pcrnot(freezing),precipitation,0.0)
    self.report(rainFall,'rain')

    self.snow = self.snow+snowFall
    #self.report(self.snow,'snow')

    potentialMelt = ifthenelse(pcrnot(freezing),temperature * Kd,0)
    self.report(potentialMelt,'pmelt')
    actualMelt = min(self.snow, potentialMelt)
    self.report(actualMelt,'amelt')

    self.snow = self.snow - actualMelt
    self.report(self.snow,'snow')

    runoffGenerated = actualMelt + rainFall
    self.report(runoffGenerated,'rg')

    discharge = accuflux(self.ldd,runoffGenerated*cellarea())
    self.report(discharge,'q')

    # reading values from a map at the observation location 
    runoffAtOutflowPoint = getCellValueAtBooleanLocation(self.observationLocation,discharge)
    #print 'modelled runoff at observation location: ', runoffAtOutflowPoint

    observedDischarge = timeinputscalar('observed.tss',1)
    observedDischargePoint = getCellValueAtBooleanLocation(self.observationLocation, observedDischarge)

    #Collecting the discharge model values and observed discharge values
    self.x.append(runoffAtOutflowPoint)
    self.y.append(observedDischargePoint)

    #Sums for the different equation parameters
    x = runoffAtOutflowPoint
    y = observedDischargePoint
    
    self.sumX = self.sumX + x
    self.sumY = self.sumY + y
    self.sumSqX = self.sumSqX + x**2
    self.sumXY = self.sumXY + (x * y)

##Default values
Kd = 0.01
FTemp = 0.0

KdLimits = [1, 20] #Degree-day factor limits [0.001, 0.02]

FTempLimits = [-10, 10] #Freezing temperature limits [-0.01, 0.01]

parameters = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #[m, b, r2, Kd, FTemp, error]
firstTime = True

#all possible combinations values for the Degree-day factor and Freezing Temperature parameters
for i in range(KdLimits[0], KdLimits[1], 1): #'Kd' (degree-day factor) calibration
  Kd = float(i) / 1000.0 
  for j in range(FTempLimits[0], FTempLimits[1], 1): #Freezing temperature calibration
    FTemp = float(j) / 1000.0

    print '\n'
    print 'FTemp = ', FTemp
    print 'Kd = ', Kd

    nrOfTimeSteps = 181
    myModel = MyFirstModel()
    dynamicModel = DynamicFramework(myModel, nrOfTimeSteps)
    dynamicModel.run()

    m, b, r2, Error = calculateRSquaredAndError(myModel)

    print '\n'
    print 'm = ', m
    print 'b = ', b
    print 'r2 = ', r2
    print 'Error = ', Error


    if(firstTime):
      parameters[0] = m
      parameters[1] = b
      parameters[2] = r2
      parameters[3] = Kd
      parameters[4] = FTemp
      parameters[5] = Error
      firstTime = False
    else:
    ##      if(r2 >= r2Min):
    ##        if(r2 > parameters[2]):
        if(Error < parameters[5]):
          parameters[0] = m
          parameters[1] = b
          parameters[2] = r2
          parameters[3] = Kd
          parameters[4] = FTemp
          parameters[5] = Error

print 'Parameters calibrated: '

print '  m = ', parameters[0]
print '  b = ', parameters[1]
print '  r2 = ', parameters[2]
print '  Kd = ', parameters[3]
print '  FTemp = ', parameters[4]
print '  Error = ', parameters[5]
