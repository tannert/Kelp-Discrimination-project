import numpy as np
from scipy import linalg as la
from matplotlib import pyplot as plt
import time
import os


def loadASCIIplot(filename):
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
    
def loadASCIIpixels(filename, coords=False):
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
    return np.arccos(v1.dot(v2)/(la.norm(v1)*la.norm(v2)))*180./np.pi
    
def plot_data(data, column):
    plt.plot(data[:,1],data[:,column])
    
def plot_mean(data):
    plt.plot(data[:,1],data[:,4])
    
def plot_normalized_mean(data):
    plt.plot(data[:,1],data[:,4]/la.norm(data[:,4]))
    
if __name__ == '__main__':
    filenames = os.listdir('.')
    filenames = filter(lambda s: s.endswith('.txt'),filenames)
    filenames.remove('bands.txt')
    
    
    for filename in filenames:
        print filename[:-4], ':', loadASCII(filename)[0].shape