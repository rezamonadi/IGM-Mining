'''
set_paramters.py : set gp_dla_detection default parameters for cddf module
'''

# physical constants
speed_of_light = 299792458

# metal lines
civ_1548_wavelength = 1548.2049;		 # CIV transition wavelength  Å
civ_1550_wavelength =  1550.77845; 		 # CIV transition wavelength  Å
# oscillator strengths
all_oscillator_strengths = [
    0.189900,
    0.094750
]

# all transition wavelengths in A
all_transition_wavelengths = [
    1.5482040e3,
    1.5507810e3
] 

loading_min_lambda = 1310; 
# convert relative velocity to redshift
kms_to_z = lambda kms : kms * 1000 / speed_of_light
z_to_kms = lambda z   : z * speed_of_light / 1000

# utility functions for redshifting
emitted_wavelengths  = lambda observed_wavelengths, z : observed_wavelengths / (1 + z)
observed_wavelengths = lambda emitted_wavelengths, z  : emitted_wavelengths  * (1 + z) 

# DLA model parameters: absorber range and model
num_lines = 3                                 # number of members of the Lyman series to use

# determines maximum z_DLA to search
max_z_cut = kms_to_z(3000)                    # max z_DLA = z_QSO - max_z_cut
max_z_dla = lambda wavelengths, z_qso : min(
    (max(wavelengths) / lya_wavelength - 1) - max_z_cut,
    z_qso - max_z_cut)

# determines minimum z_DLA to search
min_z_cut = kms_to_z(3000)                    # min z_DLA = z_Ly∞ + min_z_cut
min_z_dla = lambda wavelengths, z_qso : max(
        min(wavelengths) / lya_wavelength - 1,
        observed_wavelengths(lyman_limit, z_qso) / lya_wavelength - 1 +
        min_z_cut)

# base directory for all data
base_directory = 'data'

# utility functions for directories
processed_directory = lambda release : "{}/{}/processed".format(base_directory, release)
spectra_directory   = lambda release : "{}/{}/spectra".format(base_directory, release)
