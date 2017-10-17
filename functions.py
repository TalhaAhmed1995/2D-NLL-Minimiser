import numpy as np
import scipy.special as sp

class Functions(object):
    
    """A class where all the mathematical fits required for analysis have been created.
    This includes Probability Density Functions (PDFs) and NLL fits."""
    
    
    def __init__(self, times, errors, sigma):
        self.times = times
        self.errors = errors
        self.sigma = sigma
        
        
    def fitFunction(self, tau, times, sigma):
        """Calculates the PDF in the absence of background effects (equation (3) in submitted report).
        Parameters of tau, times and sigma need to be input."""
        
        firstTerm = 1/(2.*tau)*np.exp(0.5*((sigma/tau)**2)-(times/tau))#Exponential term of pdf
        secondTerm = sp.erfc(1/np.sqrt(2)*((sigma/tau)-times/sigma)) #erfc term of pdf
        
        pdf = firstTerm*secondTerm #Multiplying both terms to get overall pdf value
        
        return pdf #Returning the computed pdf value to be used in the nll function
        
    def nll(self, tau):
        """Calculates the NLL (equation (5) in report) in the absence of background effects. 
        The pdf (directly above) without background is used here.
        The value is dependent on the tau parameter used."""
            
        likelihood = -np.sum(np.log(self.fitFunction(tau, self.times, self.errors))) #Taking the negative sum of pdf calculations from all measurements in raw data
        
        return likelihood #Returning the computed NLL value for plotting (later)
        
    def bkgFitFunction(self, tau, a, times, sigma):
        """Calculates the PDF with inclusion of background effects (equation (6) in submitted report).
        Parameters of tau, times and sigma need to be input. 
        This time, the fraction of signal in the sample, a, is also required"""
        
        bkgFunction = 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-0.5 * ((times/sigma)**2)) #Calculating the extra Gaussian term for false readings
        
        pdfWithBkg = a * (self.fitFunction(tau, times, sigma)) + (1 - a) * bkgFunction #Calculating the pdf with background adjustment
        
        return pdfWithBkg #Returning the computed pdf value to be used in the bkgNll function
        
    def bkgNll(self, tau, a):
        """Calculates the NLL (equation (5) in report) with background readings. 
        The pdf (directly above) with background is used here.
        The value is dependent on both the tau and a parameters used."""
            
        likelihood = -np.sum(np.log(self.bkgFitFunction(tau, a, self.times, self.errors))) #Taking the negative sum of pdf calculations from all measurements in raw data
        
        return likelihood #Returning the computed NLL value for plotting (later)