
# Intro

## What is cryoET and how is it different from CryoEM


Microscopes (TEM, SEM) are used to image the sample in 2D. The sample is then tilted and imaged again. The images are combined in reciprocal space, 

The field is kind of an evolution of single particle cryoEM, but instead of imaging a single particle, slices (whole areas are imaged). 

The basic idea is "tilting" the sample while it is being showered by electrons .

The are multiple variations on the technique. 
In terms of the target being imaged -- its scale: from hosts of individual particles to whole cells, montages.
The preparation method: cryo-fib, cryo-liftout
The identification of the results: tagging, templating, combinations with light microscopy

There are various tradeoffs in resolution for these variations, but the most refined results come very close ~5-7A (within the same order of magnitude) of the SPA results, but with a much wider "coverage".




Tilt Series

# Process

Projections in real space are Slices in Fourier Space

only about -60 to 60 degrees is possible


reconstructing a regular 3d grid from tilted projections : algorithms to interpolate/estimate 


# Tomography at different scales
large targets:
- serial sectioning 
- montage tomography 
(golgi examples)

small targets:

- phmc complex conformations


# Collection and reconstruction challenges

- focusing, targeting, dose

- x,y shifts
- rotation (position of tilt axis)
- tilt angle
- magnification
- defocus


# Limitations
Resolution:
- radiation damage

to overcome - sub-tomogram averaging

choice of tilt increment (step and range) contribute to resolution
defocus
image alignment

# Identifying objects of interest


CLEM correlated light and electron microscopy
perturing the structure
template matching (visual proteomics)
    - some limitations here: computationally very expensive to due to the need for CCF
    every template. Furthermore the SNR is not high enough in the tomograms yet to image match for any proteins lighter than ~300kDa


# ----------- Downstream analysis

# The data itself

They are still in very early beta, it seems, so i think all of this might change tomorrow:

Most of the data comes from the Jensen's Lab and collaborators. Most of it is single bacterial/mammalian cells which are mutants or in which either some organelles are of interest. There are some viruses. 


The annotations are available for very few of the entries (less than 10%): although they advertise that they are working on a framework to add new annotations to older depositions. 


There is some overlap with EMPIAR, but not a lot. 


# Terms

STA
TEM
SEM


Cryo liftout



