A set of functions to compute a fiber orientation distribution function (ODF)
 from a volumetric image

Adam Rauff
Musculoskeletal Research Laboratory, University of Utah

February, 20th, 2021

The code is structured to contain a set of functions, used in the 
core computation, listed in the parent direcotry and additional sets of 
further specialized scripts listed in subdirectories.

Main Function - Fiber3D.m - Takes an image and computes and ODF
fft_Radial_filter3.m - Filters the power spectrum
getSpherPts.m - Obtain a sampling of points on the unit sphere 
MinMax.m - Performs a modified min-max normalization on ODF.
nearestpoint.m - Functions that aids in the computation of the 
	condensation of the power spectrum to the unit sphere
ReduceAmp.m - Converts the power spectrum to the unit sphere by summing along the frequencies

Description of subdirectories

Analyze_ODFs - set of functions to visualize and analyze Orientation 
distribution functions

Biomedical_Image_Data - Folder and code containing image data of biomedical 
specimen along with code that processes the data

Consistency_Analysis - Set of functions used to determine the precision of the approximation processes

Noise_Analysis - Set of functions used to determine the robustness to noise of the approximation

ODF_Resolution_Dirac_Delta - Functions used to determine the angular resolution 
afforded by the approximation by attempting to qunatify the dirac delta ODF

QBall_Algorithm - functions used in the computation of the QBall algorithm

S2 Sampling Toolbox - Functions used to sample points at a unit distance (unit sphere).
