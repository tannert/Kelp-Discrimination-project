#kelp.py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import lines as mlines
from matplotlib import markers
from matplotlib import colors
from scipy.linalg import qr, norm, inv
from ascii_tools import loadPlotASCII, loadPixelsASCII, spectral_angle
import os



def import_endmembers():
    """imports text files that have reference endmembers
    
    Returns:
        giant_endmembers (ndarray): array that has giant kelp endmembers as columns
            all text files that begin with 's_' (for south) are considered giant kelp reference endmembers
            c_monterey.txt is also considered reference endmembers
            
        bull_endmembers (ndarray): array that has bull kelp endmembers as columns
            all text files that begin with 'n_' (for north) are considered bull kelp reference endmembers
            
        giant_names (list): list of strings representing each of the giant kelp endemembers in giant_endmembers
        
        bull_names (list): list of strings representing each of the bull kelp endemembers in bull_endmembers
    """
    filenames = filter(lambda s: s.endswith('.txt'), os.listdir('.'))
    filenames.remove('bands.txt')
    n_filenames = filter(lambda s: s.startswith('n_'), filenames)
    s_filenames = filter(lambda s: s.startswith('s_'), filenames)
    giant_endmembers = []
    bull_endmembers = []
    giant_names = []
    bull_names = []
    
    for nfile in n_filenames:
        data, headers = loadPlotASCII(nfile)
        bull_endmembers += list(data.T[2:])
        bull_names += map(lambda s: nfile[:-4] + ' ' + s[s.index(':') + 2:s.index('~')] ,headers[2:])
        
    for sfile in s_filenames:
        data, headers = loadPlotASCII(sfile)
        giant_endmembers += list(data.T[2:])
        giant_names += map(lambda s: sfile[:-4] + ' ' + s[s.index(':') + 2:s.index('~')], headers[2:])
        
    data, headers = loadPlotASCII('c_monterey.txt')
    giant_endmembers += list(data.T[2:])
    giant_names += map(lambda s: 'c_monterey ' + s[s.index(':') + 2:s.index('~')], headers[2:])
    
    return giant_endmembers, bull_endmembers, giant_names, bull_names

def import_one_endmember(filename, roi_num):
    """imports just one ROI from filename
    
    Parameters:
        filename (str): name of input file
            Must be a text file
            Support is only guaranteed for files generated using ENVI's ROI tool -> Options -> Compute statistics from ROIs -> Export -> ASCII
            
        roi_num (float): yeah, that's right. it's a float.  It will work if you pass it a string instead, though
            Ex: 11.2
            
    Returns:
        vector (ndarray): the average spectrum of the ROI in question
    """
    data, header = loadPlotASCII(filename)
    i = map(lambda s: s[s.index('kelp')+5:s.index('~')], header[2:]).index(str(roi_num))
    return data.T[i+2]
    
def normalize(v):
    """accepts a vector v as an ndarray and returns a scaled version of v with a norm of 1"""
    return v.astype(float)/norm(v)

def decompose(pixel, endmembers):
    """Decomposes a pixel into a linear combination of normalized endmembers
    Parameters:
        pixel (ndarray): a 1-d numpy array of length n representing the pixel to be decomposed
        endmembers (ndarray): a list of m 1-d numpy arrays of length n, representing the endmembers
    
    Returns:
        a (ndarray): a length-m 1-d array of coefficients
        pr (ndarray): the 1-d length-n remainder vector
    """
    endmembers = np.array(endmembers).T
    E = np.array(map(normalize, endmembers.T)).T
    p = normalize(pixel)
    
    b,r = qr(E, mode = 'economic')
    
    pe = sum([p.dot(bi)*bi for bi in b.T])
    pr = p - pe
    
    R = np.array([[Ei.dot(Ej) for Ej in E.T] for Ei in E.T])
    s = np.array([pe.dot(ei) for ei in E.T])
    a = inv(R).dot(s)
    
    # print norm(a[0]*E[:,0] + a[1]*E[:,1] + pr - p)
    return a, pr
    
def quasi_mesma(pixel, set1, set2):
    """Iterates through two lists of endmembers and finds the pair that best decomposes a pixel
    
    Parameters:
        pixel (ndarray): the spectrum of the pixel to be decomposed
        set1 (list): a list of endmembers (ndarrays) all belonging to the same class
        set2 (list): a second list of endmembers (ndarrays) all belonging to the same class
        
    Returns:
        smallest_res (float): the norm of the residual resulting from unmixing with the optimal endmembers
        best_a (ndarray): a length-m 1-d array of coefficients corresponding to unmixing with the optimal endmembers
        i (int): index of the optimal endmember from set1
        j (int): index of the optimal endmember from set2
    """
    smallest_res = np.inf
    for i, u in enumerate(set1):
        for j, v in enumerate(set2):
            a, res = decompose(pixel, [u, v])
            if norm(res) < smallest_res:
                best_a = a.copy()
                best_ij = (i,j)
                smallest_res = norm(res)

    # print norm(best_a[0]*normalize(set1[best_ij[0]]) + best_a[1]*normalize(set2[best_ij[1]]) + smallest_res)
    
    return (smallest_res, best_a) + best_ij
    
def classify(filename):
    """Classifies each vector in an ROI export file according to the endmembers in import_endmembers using the projection decomposition method
    
    Parameters:
        filename (str): name of input file
            Must be a text file
            Support is only guaranteed for files generated using ENVI's ROI tool -> Options -> Compute statistics from ROIs -> Export -> ASCII
            
    Returns:
        nothing
        
    Prints: for each ROI
        - the ROI name
        - the norm of the smallest residual
        - the optimal endmembers for unmixing
        - the unmixing coefficients
        - the classification (bull kelp, giant kelp, or inconclusive)
        
    Criteria for conclusiveness:
        - the larger coefficient must be greater than .5
        - the smaller coefficient must be less than .5
       
    """
    gset, bset, gnames, bnames = import_endmembers()
    
    
    data, headers = loadPlotASCII(filename)
    bands = data[:,1]
    for k, v in enumerate(data.T[2:],2):
        res, a, i, j = quasi_mesma(v,gset,bset)
        s = headers[k]
        print filename[:-4], s[s.index(':') + 2:s.index('~')]
        print 'smallest residual was', res
        print 'it was generated by', gnames[i], 'and', bnames[j]
        print 'a =', a
        if a[0] > .5 and a[1] < .5:
            print 'giant kelp'
        elif a[1] > .5 and a[0] < .5:
            print 'bull kelp'
        else:
            print 'INCONCLUSIVE'
        print
    
def classify_by_spectrum(filename):
    """Classifies each vector in an ROI export file according to the endmembers in import_endmembers using spectral angle comparison
    
    Parameters:
        filename (str): name of input file
            Must be a text file
            Support is only guaranteed for files generated using ENVI's ROI tool -> Options -> Compute statistics from ROIs -> Export -> ASCII
            
    Returns:
        nothing
        
    Prints: for each ROI
        - the ROI name
        - the smallest angle when compared to all giant kelp endmembers
        - the smallest angle when compared to all bull kelp endmembers
        - the classification (bull kelp, giant kelp, or inconclusive)
        
    Criteria for conclusiveness:
        - the difference between the lowest bull kelp angle and the lowest giant kelp angle has to be more than 4 degrees       
    """
    gset, bset, gnames, bnames = import_endmembers()
    data, headers = loadPlotASCII(filename)
    bands = data[:,1]
    for k, pixel in enumerate(data.T[2:],2):
    
        s = headers[k]
        print filename[:-4], s[s.index(':') + 2:s.index('~')]
        
        lowest_g_dif = np.inf
        for i, u in enumerate(gset):
            dif = spectral_angle(u,pixel)
            if dif < lowest_g_dif:
                lowest_g_dif = dif
                best_i = i
                
        lowest_b_dif = np.inf
        for j, v in enumerate(bset):
            dif = spectral_angle(v,pixel)
            if dif < lowest_b_dif:
                lowest_b_dif = dif
                best_j = j        
                
        print 'smallest angle with giant kelp was', lowest_g_dif, 'with', gnames[best_i]
        print 'smallest angle with bull kelp was', lowest_b_dif, 'with', bnames[best_j]
        if abs(lowest_b_dif - lowest_g_dif) < 4:
            print 'INCONCLUSIVE'
        elif lowest_b_dif < lowest_g_dif:
            print 'bull kelp'
        elif lowest_b_dif > lowest_g_dif:
            print 'giant kelp'
        print
       
def classify_pixels(filename):
    """Classifies each vector in an ASCII image file according to the endmembers in import_endmembers using the projection decomposition method
        
        Parameters:
            filename (str): name of input file
                Must be a text file
                Support is only guaranteed for files generated using ENVI's 'Save As Type\ASCII' feature
                
        Returns:
            nothing
            
        Prints: for each pixel
            - the image coordinates (i,j)
            - the norm of the smallest residual
            - the optimal endmembers for unmixing
            - the unmixing coefficients
            - the classification (bull kelp, giant kelp, or inconclusive)
            
        Plots:
            A raster image showing classifications of the input pixels
            
        Criteria for conclusiveness:
            - the larger coefficient must be greater than .5
            - the smaller coefficient must be less than .5
            
        Raster key:
            -1: black: inconclusive
             0: white: no data
             1: blue:  bull kelp
             2: green: giant kelp
           
    """
    gset, bset, gnames, bnames = import_endmembers()
    
    data, coords = loadPixelsASCII(filename, True)
    a, b = coords[-1]
    raster = np.zeros((a+1, b+1)).astype(int)
    # bands = data[:,0]
    for k, v in enumerate(data.T[1:]):
        if np.any(v):
            res, a, i, j = quasi_mesma(v,gset,bset)
            print filename[:-4], coords[k]
            print 'smallest residual was', res
            print 'it was generated by', gnames[i], 'and', bnames[j]
            print 'a =', a
            if a[0] > .5 and a[1] < .5:
                print 'giant kelp'
                raster[coords[k]] = 2
            elif a[1] > .5 and a[0] < .5:
                print 'bull kelp'
                raster[coords[k]] = 1
            else:
                print 'INCONCLUSIVE'
                raster[coords[k]] = -1
            print
        
    print raster
    
    cmap = colors.ListedColormap(['black','white','blue','green'])
    bounds = [-1.5,-.5,.5,1.5,2.5]
    cmap_norm = colors.BoundaryNorm(bounds, cmap.N)
    plt.imshow(raster, interpolation = 'nearest', cmap=cmap, norm=cmap_norm)
    plt.show()
    
def classify_pixels_by_spectrum(filename):
    """Classifies each vector in an ASCII image file according to the endmembers in import_endmembers using spectral angle comparison
    
    Parameters:
        filename (str): name of input file
            Must be a text file
            Support is only guaranteed for files generated using ENVI's ROI tool -> Options -> Compute statistics from ROIs -> Export -> ASCII
            
    Returns:
        nothing
        
    Prints: for each pixel
        - the image coordinates (i,j)
        - the smallest angle when compared to all giant kelp endmembers
        - the smallest angle when compared to all bull kelp endmembers
        - the classification (bull kelp, giant kelp, or inconclusive)
    
    Plots:
        A raster image showing classifications of the input pixels
        
    Criteria for conclusiveness:
        - the difference between the lowest bull kelp angle and the lowest giant kelp angle has to be more than 4 degrees
        
    Raster key:
            -1: black: inconclusive
             0: white: no data
             1: blue:  bull kelp
             2: green: giant kelp
    """
    gset, bset, gnames, bnames = import_endmembers()
    data, coords = loadPixelsASCII(filename, True)
    a, b = coords[-1]
    raster = np.zeros((a+1,b+1))
    # bands = data[:,0]
    for k, pixel in enumerate(data.T[1:]):
        if np.any(pixel):
            print filename[:-4], coords[k]
            
            lowest_g_dif = np.inf
            for i, u in enumerate(gset):
                dif = spectral_angle(u,pixel)
                if dif < lowest_g_dif:
                    lowest_g_dif = dif
                    best_i = i
                    
            lowest_b_dif = np.inf
            for j, v in enumerate(bset):
                dif = spectral_angle(v,pixel)
                if dif < lowest_b_dif:
                    lowest_b_dif = dif
                    best_j = j        
                    
            print 'smallest angle with giant kelp was', lowest_g_dif, 'with', gnames[best_i]
            print 'smallest angle with bull kelp was', lowest_b_dif, 'with', bnames[best_j]
            if abs(lowest_b_dif - lowest_g_dif) < 4:
                print 'INCONCLUSIVE'
                raster[coords[k]] = -1
            elif lowest_b_dif < lowest_g_dif:
                print 'bull kelp'
                raster[coords[k]] = 1
            elif lowest_b_dif > lowest_g_dif:
                print 'giant kelp'
                raster[coords[k]] = 2
            print
            
    print raster
    
    cmap = colors.ListedColormap(['black','white','blue','green'])
    bounds = [-1.5,-.5,.5,1.5,2.5]
    cmap_norm = colors.BoundaryNorm(bounds, cmap.N)
    plt.imshow(raster, interpolation = 'nearest', cmap=cmap, norm=cmap_norm)
    plt.show()


    
   # all following functions are comparatively useless   

  
  
def test_quasi_mesma():
    """test function for quasi_mesma above"""
    gset, bset, gnames, bnames = import_endmembers()
    data, headers = loadPlotASCII('c_pebble_beach.txt')
    bands = data[:,1]
    pixel = data[:,2]
    res, a, i, j = quasi_mesma(pixel,gset,bset)
    
    print 'smallest residual was', res
    print 'it was generated by', gnames[i], 'and', bnames[j]
    print 'a =', a
    
    print norm(normalize(pixel) - a[0]*normalize(gset[i]) - a[1]*normalize(bset[j]))
    
    plt.plot(bands, normalize(gset[i]), color = 'green')
    plt.plot(bands, normalize(bset[j]), color = 'blue')
    plt.plot(bands, normalize(pixel), color = 'orange')
    plt.plot(bands, a[0]*normalize(gset[i]) + a[1]*normalize(bset[j]), color = 'purple')
    plt.legend([gnames[i], bnames[j], headers[2], 'best fit'])
    plt.show()
    
def compare_endmembers():
    """generates a plot of one bull kelp and one giant kelp endmember"""
    v = import_one_endmember('c_carmel.txt',1.1)
    u = import_one_endmember('n_fort_ross.txt',1.1)
    with open('bands.txt') as file:
        bands = map(float,file.read().split(', '))
        bands = bands[:138]
    plt.plot(bands,normalize(v),color='green')
    plt.plot(bands,normalize(u),color='blue')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Normalized Reflectance')
    plt.legend(['giant kelp','bull kelp'])
    plt.show()
    
def compare_lots_of_endmembers():
    """plots lots of endmembers at once
    right now, this is set up to plot the band 45/band 15 ratio
    """
    size = 1
    filenames = filter(lambda s: s.endswith('.txt'), os.listdir('.'))
    filenames.remove('bands.txt')
    n_filenames = filter(lambda s: s.startswith('n_'), filenames)
    s_filenames = filter(lambda s: s.startswith('s_'), filenames)
    
    blue_line = mlines.Line2D([],[],color='blue', label = 'bull kelp')
    green_line = mlines.Line2D([],[],color='green', label = 'giant kelp')
    orange_line = mlines.Line2D([],[],color='orange', label = 'monterey kelp')
    
    bands = loadPlotASCII(n_filenames[0])[0][:,1]
    print 'working'
    for nfile in n_filenames:
        data, headers = loadPlotASCII(nfile)
        for v in data.T[2:]:
            # plt.scatter(bands,v, color = 'blue', s = size)
            plt.scatter(0,1.*v[44]/v[14], color = 'blue')
            
    for sfile in s_filenames:
        data, headers = loadPlotASCII(sfile)
        for i, v in enumerate(data.T[2:],2):
            # if normalize(v)[4] > .08:
                # plt.scatter(data.T[1],v, color = 'red', s = size)
                # print sfile
                # print headers[i]
            # else:
                # plt.scatter(bands,v, color = 'green', s = size)
            plt.scatter(1,1.*v[44]/v[14], color = 'green')

    data, headers = loadPlotASCII('c_monterey.txt')
    for v in data.T[2:]:
        # plt.scatter(bands,v, color = 'green', s = size)
        plt.scatter(.5,1.*v[44]/v[14], color = 'orange')
    
    plt.ylabel('Rrs(773)/Rrs(501)')
    plt.tick_params(
        axis='x',
        which='both',
        bottom='off',
        top='off',
        labelbottom='off')
    plt.grid(axis='y')
    plt.legend(handles = [blue_line, green_line, orange_line])
    plt.show()

def plot_kelp_against_endmemebers(filename): #BROKEN NOW since I changed all the filenames
    """this one is currently broken, but only because I changed the way all the filenames are formatted"""
    filenames = os.listdir('.')
    filenames = filter(lambda s: s.endswith('.txt'),filenames)
    filenames.remove('bands.txt')
    bullkelpfilenames = filter(lambda s: s.startswith('041013'),filenames)
    giantkelpfilenames = filter(lambda s: s.startswith('041113'),filenames)
    
    bands = loadPlotASCII(bullkelpfilenames[0])[0][:,1]
    
    bullkelpdata = np.array([loadPlotASCII(e_filename)[0][:,4] for e_filename in bullkelpfilenames])
    giantkelpdata = np.array([loadPlotASCII(e_filename)[0][:,4] for e_filename in giantkelpfilenames])
    
    data, header = loadPlotASCII(filename)
    for h in header: print h
    
    legend = []
    for i,v in enumerate(data.T[2:]):
        if 'kelp 8' in header[i+2] or 'kelp 9' in header[i+2]:
            plt.plot(bands, normalize(v), linewidth = 2)
            legend.append(header[i+2])
            
    plt.legend(map(lambda s: s[6:s.index('~')],legend))
    
    blue_line = mlines.Line2D([],[],color='blue', label = 'bull kelp')
    green_line = mlines.Line2D([],[],color='green', label = 'giant kelp')
    
    for bkelp in bullkelpdata:
        plt.plot(bands, normalize(bkelp), color = 'blue')
    for gkelp in giantkelpdata:
        plt.plot(bands, normalize(gkelp), color = 'green')
    # plt.legend(handles = [blue_line, green_line])
        
    plt.show()
    
def concrete_grass_test():
    """tests the decompose function using the concrete and grass endmembers"""
    data, header = loadPlotASCII('grass_concrete.txt')
    concrete1 = data[:,2]
    concrete2 = data[:,3]
    grass = data[:,4]
    
    a, pr = decompose(concrete2,[grass, concrete1])
    print 'a =', a
    print '||pr|| =', norm(pr)
    print
    print 
    # plt.plot(data[:,1],abs(pr))
    # plt.plot(data[:,1],normalize(lagoon))
    # plt.plot(data[:,1],normalize(grass))
    # plt.plot(data[:,1],normalize(deepwater))
    # plt.legend(['residual','lagoon','grass','deepwater'])
    # plt.show()
    # print
    print 'angle between concretes:', spectral_angle(concrete1, concrete2)
    print 'angle between concrete1 and grass:', spectral_angle(concrete1, grass)
    print 'angle between concrete2 and grass:', spectral_angle(concrete2, grass)

if __name__ == '__main__':
    """this code runs when you run the file from the command line, as opposed to importing it as a module
    I use it for testing
    """
    classify_pixels('pebble beach pixel group 11.txt')
    classify_pixels_by_spectrum('pebble beach pixel group 11.txt')
    
# Randii Wessen was our JPL tour guide
    
    