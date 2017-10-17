import numpy as np

class Data(object):
    
    """A class solely for extracting data from text files."""
    
    def __init__(self, filename, mode):
        self.filename = filename
        self.mode = mode
        
    def readData(self, size=None):
        """Reads the input text file and extracts the data in the form of arrays."""
    
        with open(self.filename, self.mode) as f: #Opening the file so it can be interpreted
            
            data = f.readlines() #Scanning all distinct lines in the file
    
            times = []
            errors = []
            
            if size is None: #If no size is specified, all the time and error values are extracted
                for line in data:
                    measurement = line.split()
                    times.append(float(measurement[0])) 
                    errors.append(float(measurement[1]))
            else:
                for line in data: #If size is specified, a limited amount of time and error values are extracted
                    if len(times) < size:
                            measurement = line.split()
                            times.append(float(measurement[0]))
                            errors.append(float(measurement[1]))
                            
                    else:
                        break
            
            #Converting the time and error lists into numpy arrays for easier manipulation in the future
            times = np.asarray(times)
            errors = np.asarray(errors)
            
            return times, errors #Returning these arrays to be used elsewhere