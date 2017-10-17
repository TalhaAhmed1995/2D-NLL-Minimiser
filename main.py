import numpy as np
from scipy.integrate import trapz

from data import Data
from minimiser import Minimiser
from functions import Functions
import plotting as pt

"""This is the main module, where everything is run. It imports all of the other 4 modules:
1) data.py
2) minimiser.py
3) functions.py
4) plotting.py

All plots and calculated values as featured in the report are outputted here.
Simply clicking run will present the results.
NOTE: AFTER RUNNING, IT CAN TAKE UP TO 35 SECONDS FOR RESULTS TO APPEAR."""
        

### RETRIEVING DATA FROM FILE ###
times, errors = Data('lifetime.txt', 'r').readData() #Decay times and errors in picoseconds
meanSigma = np.mean(errors) #Calculating mean of errors from data
sortedTimes = np.sort(times) #Sorting times in ascending order

### CREATING FUNCTION AND MINIMISER OBJECTS ###
funcs = Functions(times, errors, meanSigma)
minim = Minimiser(times, errors, meanSigma) 

### CALCULATING VALUES AND PLOTTING LOCATIONS FOR 1D FIT (WITHOUT BACKGROUND) ###
tauMin, nllMin, x, y = minim.minimiseNll([0.3, 0.4, 0.5], 1e-5) #Minimising NLL
taus = np.linspace(tauMin-0.1, tauMin+0.1, 100) #Tau values in region of minimum
likelihoods = np.array([funcs.nll(tau) for tau in taus]) #Calculating NLL over these tau values
pdf = funcs.fitFunction(tauMin, sortedTimes, meanSigma) #Calculating pdf function without background
area = trapz(pdf, sortedTimes) #Calculating area of pdf using trapezium rule
#Calculating errors
posError = minim.posError(tauMin, 1e-5)
negError = minim.negError(tauMin, 1e-5)
meanError = np.mean([posError, negError])
parabError = minim.parabError(x, y)
    
### CALCULATING LOG OF ERROR IN AVERAGE LIFETIME WITH LOG OF DATA SUBSET SIZE ###
sizes, sizeErrors = minim.errorVSReadings()

### CALCULATING VALUES AND PLOTTING LOCATIONS FOR 2D FIT (WITH BACKGROUND) ###
bkgNllMin, bkgTauMin, bkgAMin = minim.bkgMinimiseNll(0.4, 0.9)#Minimising 2D NLL
bkgTaus = np.linspace(bkgTauMin - 0.05, bkgTauMin + 0.05, 100)  #Tau values in region of minimum
bkgAs = np.linspace(bkgAMin - 0.05, 0.999, 100)  #a values in region of minimum
fractionBkg = 1 - bkgAMin #Calculating fraction of false readings in sample
pdf2 = funcs.bkgFitFunction(bkgTauMin, bkgAMin, sortedTimes, meanSigma) #Calculating pdf function with background
area2 = trapz(pdf2, sortedTimes) #Calculating area of new pdf using trapezium rule
#Calculating errors
TauAErrors = minim.bkgError(bkgTaus, bkgAs, bkgTauMin, bkgAMin, bkgNllMin, funcs.bkgNll)
meanTauError = np.mean([TauAErrors[0], TauAErrors[1]])
meanAError = np.mean([TauAErrors[2], TauAErrors[3]])
 

### CREATING AND DISPLAYING PLOTS ###
pt.plotContour(bkgTaus, bkgAs, funcs.bkgNll, levels = np.arange(bkgNllMin+0.5, bkgNllMin+0.5+100, 10))
pt.plot3D(bkgTaus, bkgAs, funcs.bkgNll)
intersect, m, c = pt.plotErrorsVSReadings(sizes, sizeErrors)
pt.plotNLL(taus, likelihoods, tauMin, nllMin)
pt.plotHist(times, bins=100, sortedTimes=sortedTimes, pdfs=[pdf, pdf2])

### PRINTING RESULTS ###
print "Area under initial PDF (without background):", area
print "Area under refined PDF (with background):", area2
print "NLL at Minimum 1D:", nllMin
print "Tau at Minimum 1D:", tauMin, "picoseconds"
print "Positive Tau Error 1D:", posError, "picoseconds"
print "Negative Tau Error 1D:", negError, "picoseconds"
print "Mean Tau Error 1D:", meanError, "picoseconds"
print "Last Parabolic Estimate Error of Tau:", parabError, "picoseconds"
print "Intersect of extrapolation of gradient", m, "and intercept", c,  "gives approximately", int(10**intersect[0]), "readings required for accuracy of 0.001ps"  
print "NLL at Minimum 2D:", bkgNllMin
print "Tau at Minimum 2D:", bkgTauMin, "picoseconds"
print "Mean Tau Error from contour at NLLMin + 0.5:", meanTauError, "picoseconds"
print "a at Minimum 2D:", bkgAMin
print "Mean a Error from contour at NLLMin + 0.5:", meanAError
print "Number of Background Readings out of 10,000:", int(np.round(fractionBkg*10000))