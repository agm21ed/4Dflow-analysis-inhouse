# 4Dv2D-Flow-Analysis
Code for slicing 4D flow data with a 2D plane and extracting flow information from this slice

1) moveFlowImages_4Dv2D.m takes the subject DICOMs from the study folder, and moves them to my own folder - amorganPhD. It keeps the 4D data separated into Magnitude, Velocity (Foot-head), Velocity (Anterior-posterior), and Velocity (Left-right) for that subject, while separating the images within these folders into individual timeframes. i.e. every brain slice for TF1, TF2 and so on. These DICOMs are then converted into a nifti slab within each TF folder.

This file also extracts header information from the first image of each set of velocity data in the study folder. This includes venc, nominal interval (RR interval) and time resolution - for later use.

2) process_4Dv2D.m is the code used to extract a slice of interest from the 4D slab of data. At the moment, the slice needs a name inputting (e.g. 'ACAs'), radius (which seems to work consistently at 150), and the slice centrepoint coordinates and normal vector stored in the vesselcoords+normal excel file.

With this information, the code loads in the slab of data (first velocity 1,2,3, then mag) and uses the extractSlice function (CITE) to extract the relavent slice based on the coordinates and normal, this is then saved for each timeframe for each velocity direction/mag.

The slice across all TFs is saved as a merged nifti and derived images are created for mask drawing later (MeanAbs, Mean, Max and MaxAbs for all 3 velocity directions; Mean and Max for magnitude)

3) processFlowData_4Dv2D_AM.m is a script that runs setFlowParameters_4Dv2D_AM.m and then carries out pipeline_extractFlow_4Dv2D_AM to produce flow graphs for each vessel in the slice, for each velocity direction
