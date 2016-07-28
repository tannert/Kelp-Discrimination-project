import numpy as np
from scipy import linalg as la
from matplotlib import pyplot as plt
import time
import os


def loadPlotASCII(filename):
    """loads mean ROI spectra from an ENVI ASCII-type spectrum output file
    
    Parameters:
        filename (str): name of input file
            Must be a text file
            Support is only guaranteed for files generated using ENVI's ROI tool -> Options -> Compute statistics from ROIs -> Export -> ASCII
    
    Returns:
        data (ndarray):
            a two-dimensional array where:
                - each row is a band
                - the first column is the band numbers
                - the second column has the wavelengths of each band
                - each of the rest of the columns is a pixel spectrum
    """
    with open('bands.txt') as file:
        bands = map(float,file.read().split(', '))
        band_num = dict(zip(bands,range(1,1+len(bands))))
    
    with open(filename) as file:
        lines = file.read().split('\n')
    
	M = len(str.split(lines[-2]))
	
    N = len(lines) - M - 2
    data = np.zeros((N,M+1))
    for i in range(N):
            line = str.split(lines[i+M+1])
            data[i] = [band_num[float(line[0])]] + line

    headers = map(lambda s: s[s.index(':')+2:],lines[1:1+M])
    headers = ['Band Number'] + headers
    return data, headers
    
def loadPixelsASCII(filename, coords=False):
    """loads pixel spectra from an ENVI ASCII-type image file
    
    Parameters:
        filename (str): name of input file.
            Must be a text file
            Support is only guaranteed for files generated using ENVI's 'Save As Type\ASCII' feature
        coords (bool) (optional):
            if True, the program also returns vectorcoords (see below)
            
    Returns:
        vectors (ndarray):
            a two-dimensional array where each column is a pixel spectrum, and each row is a band
        vectorcoords (list) (optional):
            a list of tuples indicating the image coordinates corresponding to the pixels being loaded
    """
    # with open('bands.txt') as file:
        # bands = map(float,file.read().split(', '))
        # band_num = dict(zip(bands,range(1,1+len(bands))))
        
    with open(filename) as file:
        lines = map(str.split, file.readlines())
        
    num_samples = int(lines[2][3])
    num_lines = int(lines[2][6])
    num_bands = int(lines[2][9])
    
    lines = np.array(map(lambda l: map(int, l), lines[5::2]))
    assert lines.shape == (num_lines, num_bands*num_samples)
    
    vectors = []
    vectorcoords = []
    for i, line in enumerate(lines):
        for j in range(num_samples):
            vectors.append(line[j*num_bands:(j+1)*num_bands])
            vectorcoords.append((i,j))
            
    if coords:
        return np.array(vectors).T, vectorcoords
    else:
        return np.array(vectors).T
            
        

def spectral_angle(v1,v2):
    """computes the spectral angle between two vectors of identical length
    
    vectors must be numpy arrays (type ndarray)
    """
    return np.arccos(v1.dot(v2)/(la.norm(v1)*la.norm(v2)))*180./np.pi
    
def plot_data(data, column):
    """plots, but does not show, one spectrum from a data array.
    To show the plot, use plt.show() after one or more calls to this function
    
    Parameters:
        data (ndarray): first output element from loadPlotASCII
        column (int): the index of the column you want to plot
    """
    plt.plot(data[:,1],data[:,column])
    
def plot_mean(data):
    """same as plot_data above, but plots column 4. DEPRECATED"""
    plt.plot(data[:,1],data[:,4])
    
def plot_normalized_mean(data):
    """same as plot_mean, but normalizes the vector before plotting"""
    plt.plot(data[:,1],data[:,4]/la.norm(data[:,4]))
    
if __name__ == '__main__':
    """this code runs when you run the file from the command line, as opposed to importing it as a module
    I use it for testing
    Right now, it prints the name of each text file in the folder along with the dimensions of the array that is formed when importing it
    """
    filenames = os.listdir('.')
    filenames = filter(lambda s: s.endswith('.txt'),filenames)
    filenames.remove('bands.txt')
    
    
    for filename in filenames:
        print filename[:-4], ':', loadASCII(filename)[0].shape
        