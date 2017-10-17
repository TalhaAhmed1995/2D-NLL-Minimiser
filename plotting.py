from matplotlib import pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

"""This module is strictly used for plotting various graphs.
This includes creating and plotting figures of histograms, 2D plots, 3D surface plots and contour plots"""

def plotHist(times, bins=100, sortedTimes=None, pdfs=None,):
    """Plots a histogram of the raw data. Also includes pdf plots, if specified"""
    
    plt.figure()
    if pdfs is None and sortedTimes is None:
        plt.hist(times, bins=bins, normed=True) #Plotting just the histogram of raw if no pdfs input
    else:
        plt.hist(times, bins=bins, color= 'orange', normed=True, label='Raw Data') #Plotting normalised histogram with specified number of bins
        plt.plot(sortedTimes, pdfs[1], color='cyan', linewidth=2, label='Tau = 0.4097ps') #Plotting pdf with background correction
        plt.plot(sortedTimes, pdfs[0], color='green', linewidth=2, label='Tau = 0.4045ps') #Plotting pdf with no background correction
    plt.title('Histrogram of decay times for 10,000 D0 Samples with Sigma=0.282')
    plt.xlabel('Decay Time (picoseconds)')
    plt.ylabel('Probability Density')
    plt.grid()
    plt.legend()
    plt.show()
    
    
def plotPDF(times, pdf):
    """Plots a PDF over time."""
    
    # Plotting the input PDF function against time
    plt.figure()
    plt.plot(times, pdf)
    plt.title('Fit Function')
    plt.xlabel('Decay Time (picoseconds)')
    plt.ylabel('PDF')
    plt.grid()
    plt.show()
    
def plotNLL(taus, likelihoods, tauMin = None, nllMin = None):
    """Plots a 1D NLL function against varying parameter values ie tau."""
    
    #Plotting NLL function against tau
    plt.figure()
    plt.plot(taus, likelihoods, label='NLL fit')
    if tauMin is not None and nllMin is not None:
        plt.plot(tauMin, nllMin, 'o', label='NLL minimum') # The minimum of the NLL is also indicated on the plot, if requested
    plt.title('Negative Log Likelihood')
    plt.xlabel('Mean Lifetime, Tau (picoseconds)')
    plt.ylabel('NLL')
    plt.grid()
    plt.legend()
    plt.show()
    
def plotErrorsVSReadings(sizes, sizeErrors, lowerLimit=2.5, upperLimit=6.1):
    """Plots variation of error in tau value against number of readings in data set (on log plot).
    Also extroplates the data to find the number of readings required for accuracy of 1fs."""
    
    sizes = np.log10(sizes) #Taking log of number of readings in data
    sizeErrors = np.log10(sizeErrors) #Taking log of the error on the tau estimate from each data subset
    
    A = np.vstack([sizes, np.ones(len(sizes))]).T
    m, c = np.linalg.lstsq(A, sizeErrors)[0] #Obtaining the gradient and intercept values for a least sqaures linear regression fit to the plot
    
    xValues = []
    yValues = []
    
    for x in np.arange(lowerLimit, upperLimit, 1e-5): #Extrapolating the linear fit to data sizes beyond that provided
        y = m*x + c
        
        xValues.append(x)
        yValues.append(y)
        
    yValues = np.asarray(yValues)
    xValues = np.asarray(xValues)
    
    differences = abs(yValues - (-3.0)) #Calculating difference between each point on linear fit and the required accuracy
    intersectIndex = np.argmin(differences) #The smallest difference between the linear fit and the required accuracy gives a good approximation of the intersect
    intersect = np.array([xValues[intersectIndex], yValues[intersectIndex]]) #Determining the intersect coordinates on the plot
    
    #Plotting log-log graph of error against number of readings, along with approximate fit
    plt.figure()
    plt.plot(sizes, sizeErrors, 'ro', label='Data subsets')
    plt.plot(xValues, yValues, 'b', label='Fitted, extrapolated line')
    plt.axhline(-3.0, color='g', linestyle='--', label='Required Accuracy')
    plt.plot(intersect[0], intersect[1], 'bo', label='Intersect')
    plt.title('Error in Optimum Tau value against Number of Readings')
    plt.xlabel('Base 10 log(Number of Data Readings)')
    plt.ylabel('Base 10 log(Error (picoseconds))')
    plt.grid()
    plt.legend()
    plt.show()
    
    return intersect, m, c
    
    
def plot3D(tau, a, function):
    """Plots a 2D NLL function (in 3D space) against varying parameter values ie tau and a."""
    
    #Calculating NLL values for each (tau, a) coordinate pair
    TAU, A = np.meshgrid(tau, a)
    ls = np.array([function(tau, a) for tau,a in zip(np.ravel(TAU), np.ravel(A))])
    LS = ls.reshape(TAU.shape)
    
    #Plotting the NLL surface in 3D space
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(TAU, A, LS, rstride=10, cstride=10, cmap=cm.coolwarm, linewidth=0.1, antialiased=True) #cmap argument provides visual representation of relative surface height
    plt.title('3D surface plot of NLL for varying tau and a values')
    ax.set_xlabel('tau (picoseconds)')
    ax.set_ylabel('a')
    ax.set_zlabel('NLL')
    plt.show()
    
def plotContour(tau, a, function, levels =None):
    """Plots contours of the 2D NLL function against varying parameter values ie tau and a."""
    
    #Calculating NLL values for each (tau, a) coordinate pair
    TAU, A = np.meshgrid(tau, a)
    ls = np.array([function(tau, a) for tau,a in zip(np.ravel(TAU), np.ravel(A))])
    LS = ls.reshape(TAU.shape)
    
    #Plotting the 2D NLL surface contours
    plt.figure()
    if levels is None:
        cp = plt.contour(TAU, A, LS) #Automatically plots contours over range of tau and a values inputted
    else:
        cp = plt.contour(TAU, A, LS, levels=levels) #The contour levels required can also be specified
    plt.clabel(cp, inline=True, 
            fontsize=10)
    plt.title('Contour Plot of NLL in proximity of minimum')
    plt.xlabel('tau (picoseconds)')
    plt.ylabel('a')
    plt.show()      