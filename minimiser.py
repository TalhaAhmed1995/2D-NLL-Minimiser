import numpy as np
from functions import Functions 
from data import Data
from matplotlib import pyplot as plt

class Minimiser(object):
    
    """A class containting all minimisation and error calculation algorithms.
    This is for both with and without background adjustments made."""
    
    
    def __init__(self, times, errors, sigma):
        self.times = times
        self.errors = errors
        self.sigma = sigma
        self.funcs = Functions(times, errors, sigma)


    def minimiseNll(self, x=[0.3, 0.4, 0.5], tol = 1e-5):
        """Calculates the 1D NLL minimum in the absence of background effects using the parabolic method.
        3 initial points and a tolerance level need to be specificed."""
        
        x = np.asarray(x) # Converting the 3 input x points on the parabola into an array
        y = []
        
        for value in x:
            y.append(self.funcs.nll(value)) #Calculating the 1D NLL at each of these points
            
        np.asarray(y)
        
        difference = 2 * tol 
        
        while difference > tol: # The process is stopped once the difference between consecutive x points is below a specified tolerance
            
            numerator = (x[2]**2 - x[1]**2)*y[0] + (x[0]**2 - x[2]**2)*y[1] + (x[1]**2 - x[0]**2)*y[2]
            denominator = (x[2] - x[1])*y[0] + (x[0] - x[2])*y[1] + (x[1] - x[0])*y[2]
            
            newX = 0.5 * float(numerator)/float(denominator) #Calculating the minimum of the parabola (equation (7) in submitted report)
            
            newY = self.funcs.nll(newX) #Calculating NLL value are new parabola minimum
            
            difference = abs(newX - x[np.argmin(y)]) #If the difference is below the tolerance, the process is stopped
            
            if newY < y[np.argmax(y)]: #Replacing highest y (and corresponding x) value with new one
                x[np.argmax(y)] = newX
                y[np.argmax(y)] = newY
            
                
        xMin = x[np.argmin(y)] #Minimisation of NLL has been achieved 
        yMin = y[np.argmin(y)]
        
        return xMin, yMin, x, y #Returning x,y values are the minimum, as well as the list of the last parabolic estimate to be used for error calculations
        
        
    def parabError(self, x, y):
        """Calculates the error in the optimum value from the minimised NLL using the curvature of the last parabolic estimate from the parabolic algorithm (above).
        The last parabolic estimate x,y values need to be input."""
        
        numerator = (x[1] - x[0]) * (x[2] - x[0]) * (x[2] - x[1])
        
        denominator = (x[2] - x[1])*y[0] + (x[0] - x[2])*y[1] + (x[1] - x[0])*y[2]
        
        error = np.sqrt(0.5 * float(numerator)/float(denominator)) #Calculating the error from the parabola (equation (8) in submitted report)
        
        return error #Returning the error calculated
        
    def posError(self, tauMin, tol=1e-5):
        """Calculates the positive error in the optimum tau value from the minimised NLL by moving at small increments until surpassing NLLMin + 0.5
        The increments are specified by the tolerance."""
        
        nllMin = self.funcs.nll(tauMin) #Calculating NLL from the optimum tau value
        
        shiftNll = nllMin + 0.5 #Calculating NLL 0.5 above the minimum
        
        taus = []
        for i in np.arange(tauMin, 2*tauMin, tol): #Looping over increasing tau values from the optimal estimate
            
            nllValue = self.funcs.nll(i) #Calculating nll at each tau value, one at a time
            
            if nllValue > shiftNll: #If the nll calculation surpasses the 0.5 shift, stop iterating
                break
            else:
                taus.append(i) #If the nll value is below NLLMin + 0.5, keep the tau value
                
        positiveError = abs(taus[-1] - tauMin) #Calculating the abs difference between the last tau value that didn't surpass NLLMin + 0.5 and mean lifeitme estimate
                
        return positiveError #Returning the positive tau error calculated
        
        
    def negError(self, tauMin, tol=1e-5):
        """Calculates the negative error in the optimum tau value from the minimised NLL by moving at small increments until surpassing NLLMin + 0.5
        The increments are specified by the tolerance."""
        
        nllMin = self.funcs.nll(tauMin) #Calculating NLL from the optimum tau value
        
        shiftNll = nllMin + 0.5 #Calculating NLL 0.5 above the minimum
        
        taus = []
        for i in np.arange(tauMin, 0.5*tauMin, -tol): #Looping over decreasing tau values from the optimal estimate
            
            nllValue = self.funcs.nll(i)
            
            if nllValue > shiftNll: #If the nll calculation surpasses the 0.5 shift, stop iterating
                break
            else:
                taus.append(i) #If the nll value is below NLLMin + 0.5, keep the tau value
                
        negativeError = abs(taus[-1] - tauMin) #Calculating the abs difference between the last tau value that didn't surpass NLLMin + 0.5 and mean lifeitme estimate
        
        return negativeError #Returning the negative tau error calculated
        
    def errorVSReadings(self, lowerRange=1000, upperRange=10000):
        """Calculates the error in the optimum tau value for varying subsets of the data provided.
        The range of data sizes required is specified."""
        
        sizes = []
        sizeErrors = []
        
        for i in range(lowerRange, upperRange+100, 100): #Iterating of varying number of readings from the data
            
            times, errors = Data('lifetime.txt', 'r').readData(size=i) #Reading text file up to specified size limit
        
            meanSigma = np.mean(errors) #Calculating the mean sigma from all the errors of the data subset
            minim = Minimiser(times, errors, meanSigma) #Instantiating Minimiser object for the required data subset
            
            tauMin, nllMin, x, y = minim.minimiseNll([0.3, 0.4, 0.5], 1e-5) #Calculating the NLL minimum, along with the last parabolic estimate
            parabError = minim.parabError(x, y) #Calculating error in best estimate of tau using parabolic approach
            
            sizes.append(i)
            sizeErrors.append(parabError)
            
        np.asarray(sizes) 
        np.asarray(sizeErrors)
        
        return sizes, sizeErrors #Returning an array of the varying data sizes, along with an array of associated errors
        
    def bkgMinimiseNll(self, tauStart, aStart, h=1e-5, alpha=1e-5, tol=1e-5):
        """Calculates the 2D NLL minimum including background effects using the gradient method.
        Starting tau and a are input. Step size, alpha and tolerance level also need to be specificed."""
        
        x0 = np.array([tauStart, aStart]) #Starting coordinates in array
        y0 = self.funcs.bkgNll(x0[0], x0[1]) #2D NLL evalated at initial coordinates
        
        tauGrad = 1/float(h) * (self.funcs.bkgNll(x0[0]+h, x0[1]) - self.funcs.bkgNll(x0[0], x0[1])) #Gradient in tau and a calculated using FDS
        aGrad = 1/float(h) * (self.funcs.bkgNll(x0[0], x0[1]+h) - self.funcs.bkgNll(x0[0], x0[1])) 
        
        grads = np.array([tauGrad, aGrad]) #Gradients stored in array
        
        x1 = x0 - alpha*grads #New coordinates found by shifting in oppposite direction to gradient on contour
        y1 = self.funcs.bkgNll(x1[0], x1[1]) #2D NLL evalated at new coordinates
        
        difference = abs(y1 - y0) #Difference between new and previous NLL value calculated
        
        while difference > tol: #Entire process continuously repeated until difference in NLL is below tolerance
            
            x0 = x1
            y0 = self.funcs.bkgNll(x0[0], x0[1])
        
            tauGrad = 1/float(h) * (self.funcs.bkgNll(x0[0]+h, x0[1]) - self.funcs.bkgNll(x0[0], x0[1]))
            aGrad = 1/float(h) * (self.funcs.bkgNll(x0[0], x0[1]+h) - self.funcs.bkgNll(x0[0], x0[1]))
            
            grads = np.array([tauGrad, aGrad])
            
            x1 = x0 - alpha*grads
            y1 = self.funcs.bkgNll(x1[0], x1[1])
            
            difference = abs(y1 - y0)
            
        return y1, x1[0], x1[1] #NLL at minimum and (tau, a) coordinates at minimum returned
        
    def bkgError(self, tau, a, tauMin, aMin, funcMin, function):
        """Calculates the errors in the optimum tau and a values from the minimised 2D NLL by analysing contour at NLLMin + 0.5."""
        
        #Calculating NLL values for each (tau, a) coordinate pair
        TAU, A = np.meshgrid(tau, a)
        ls = np.array([function(tau, a) for tau,a in zip(np.ravel(TAU), np.ravel(A))])
        LS = ls.reshape(TAU.shape)
        
        cp = plt.contour(TAU, A, LS, levels = [funcMin + 0.5]) #Creating contour level at NLLMin + 0.5
        p = cp.collections[0].get_paths()[0]
        v = p.vertices
        #Lists of all tau  and a values on this contour
        tauPoints = v[:,0]
        aPoints = v[:,1]
        
        plt.close()
        
        plusTau = []
        minTau = []
        plusA = []
        minA = []
        
        #Collecting tau values above and below parameter estimate
        for i in tauPoints:
            if i > tauMin:
                plusTau.append(i)
            else:
                minTau.append(i) 
        #Collecting a values above and below parameter estimate
        for i in aPoints:
            if i > aMin:
                plusA.append(i)
            else:
                minA.append(i)
                
        #Positive errors are calculated by taking the average of all values above parameter estimate and taking difference
        #Negative errors are calculated by taking the average of all values below parameter estimate and taking away from this optimum value
        posTauError = np.mean(plusTau) - tauMin
        negTauError = tauMin - np.mean(minTau)
        posAError = np.mean(plusA) - aMin
        negAError = aMin - np.mean(minA)
        
        return posTauError, negTauError, posAError, negAError # Returning positive and negative errors for both tau and a
        
            