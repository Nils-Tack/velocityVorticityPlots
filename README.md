# velocityVorticityPlots
This repository contains a suite of Matlab scripts and functions to plot velocity and vorticity fields from PIV data. In addition, the script performs basic image processing to extract black and white masks used to mask data when an object is present in the field of view.

1) Getting started
To work with PIV data, you must first create a folder to store the raw PIV images (preferably .jpg) and the PIV data (in .txt. format).
To work with the data, create the main folder called 'Data' and subfolders called 'Figures', 'images', 'masksBW','outlines', 'PIVclean', and 'PIVraw'.
Load the PIV images in the 'images' subfolder. Load the .txt PIV data in the 'PIVraw' subfolder.
The rest of the subfolders will serve to store the output of the master script (for plots and masks).

2) Run the master script
The master script to produce masks and plot velocity and vorticity plots is called 'masterVeloVortPlot'.
All the custom functions used to pre-process and threshold the images for mask-making must be accessible to run the master script. 
