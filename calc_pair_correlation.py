import numpy as np 
import matplotlib.pyplot as plt 
from ctypes import CDLL, c_int, c_double,  c_float, byref, c_char_p 

# load the dynamic module:

dcd_lib = CDLL("./libdcd.so")

pair_corr_lib = CDLL("./lib_pair_correlation.so")

def string_to_ctypes_string(string):

    strlength = c_int(len(string))

    string_to_bytes = string.encode()

    string = c_char_p(string_to_bytes)

    return string, strlength

def call_read_dcd_header(dcdfile):

    # declare c types varibles:

    dcdfile, strlength = string_to_ctypes_string(dcdfile)

    total_atoms = c_int()

    total_frames = c_int()

    dcd_lib.call_dcd_header(dcdfile,
                            byref(strlength),
                            byref(total_atoms),
                            byref(total_frames))

    return total_atoms.value, total_frames.value


def call_read_dcd_xyz_box(dcdfile, current_frame, total_atoms, return_numpy):

    # declare variables:
    
    dcdfile, strlength = string_to_ctypes_string(dcdfile)
    
    current_frame = c_int(current_frame)

    total_atoms = c_int(total_atoms)

    # box = (c_double*3)()
    # "np.ctypeslib.as_ctypes" converts a numpy array into ctypes array
    box = np.ctypeslib.as_ctypes(np.zeros(3, dtype=np.float64))

    # Don't worry about ctypes array dimensions.  Only array size matters!
    xyz = np.ctypeslib.as_ctypes(np.zeros(total_atoms.value*3,
                                          dtype=np.float32))

    # current_frame = c_int(current_frame)
    dcd_lib.call_dcd_traj(dcdfile,
                          byref(strlength),
                          byref(total_atoms),
                          byref(current_frame),
                          xyz,
                          box)

    # return as Numpy for downstream analysis in Python
    if (return_numpy):

        # switch from column-major order in Fortran to row-major order in C/C++/Python
        # conver the ctypes array into a numpy array
        xyz = np.ctypeslib.as_array(xyz).reshape((total_atoms.value, 3))

        box = np.ctypeslib.as_array(box)

        return xyz, box

    # return as c data type as an intermdiate (more efficient) 
    else:

        # convert to double precision
        xyz = np.ctypeslib.as_ctypes(np.ctypeslib.as_array(xyz)
                                     .astype(np.float64))

        return xyz, box

def compute_pair_correlation(total_atoms, xyz, box, cutoff, num_bins):  

    rdf_histogram = np.ctypeslib.as_ctypes(np.zeros(num_bins, dtype=np.float64)) 

    total_atoms = c_int(total_atoms)

    cutoff = c_double(cutoff)

    num_bins = c_int(num_bins)

    volume = c_double()

    pair_corr_lib.call_homo_pair_dist_his(byref(total_atoms),
                                          byref(cutoff),
                                          byref(num_bins),
                                          xyz,
                                          box,
                                          rdf_histogram,
                                          byref(volume))

    return np.ctypeslib.as_array(rdf_histogram), volume.value

def compute_r_bins(cutoff, num_bins):

    r_interval = cutoff/num_bins 

    r_mid_point = np.zeros(num_bins)

    for i in range(num_bins): 

        r_mid_point[i] = r_interval*0.5 + i*r_interval 

    return r_mid_point  

def normalize_histogram(histogram, bulk_density,num_bins, total_atoms, cutoff, num_configs): 

    histogram = np.ctypeslib.as_ctypes(histogram.astype(np.float64))

    r = compute_r_bins(cutoff, num_bins) 

    gr = np.ctypeslib.as_ctypes(np.zeros(num_bins, dtype=np.float64))

    cutoff = c_double(cutoff)

    num_bins = c_int(num_bins) 

    total_atoms = c_int(total_atoms) 

    num_configs = c_int(num_configs) 

    bulk_density = c_double(bulk_density) 

    pair_corr_lib.call_normalize_histogram(histogram,
                                           byref(num_bins), 
                                           byref(cutoff),
                                           byref(total_atoms),
                                           byref(num_configs),
                                           byref(bulk_density),
                                           gr)


    return r, np.ctypeslib.as_array(gr) 

dcdfile = "traj.dcd"

total_atoms, total_frames = call_read_dcd_header(dcdfile)

cutoff = 9.0

num_bins = 150

cum_rdf_histogram = np.zeros(num_bins) 

rho = 0 

# loop over each frame 
for current_frame in range(1,total_frames+1):

    # extract the coordinates and box
    xyz, box = call_read_dcd_xyz_box(dcdfile, current_frame, total_atoms, return_numpy=False)
    
    # compute distance histogram:
    rdf_histogram, volume = compute_pair_correlation(total_atoms, xyz, box, cutoff, num_bins)

    # accumulate distance histogram 
    cum_rdf_histogram += rdf_histogram  

    # accumlate density 
    rho += total_atoms/volume 

# compute density
bulk_density = rho/total_frames

# normalize the distance histogram
r, gr = normalize_histogram(cum_rdf_histogram, bulk_density,num_bins, total_atoms, cutoff, total_frames)    

# compare with Figure 3 in https://pubs.acs.org/doi/pdf/10.1021/jp805227c 

plt.plot(r, gr) 

plt.title("mW 298K 1bar")

plt.ylim([0,2.5]) 

plt.xlim([2,9])

plt.ylabel("g(r)")

plt.xlabel(r"r($\AA$)")

plt.savefig("gr_mW_298K.png")

plt.show() 


