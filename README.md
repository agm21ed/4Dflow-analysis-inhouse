# 4Dv2D-Flow-Analysis
Code for slicing 4D flow data with a 2D plane and extracting flow information from this slice

1) moveFlowImages_4Dv2D.m takes the subject DICOMs from the study folder, and moves them to your own research folder. It keeps the 4D data separated into Magnitude, Velocity (Foot-head), Velocity (Anterior-posterior), and Velocity (Left-right) for that subject, while separating the images within these folders into individual timeframes. i.e. every brain slice for TF1, TF2 and so on. These DICOMs are then converted into a nifti slab within each TF folder.

This file also extracts header information from the first image of each set of velocity data in the study folder. This includes venc, nominal interval (RR interval) and time resolution - for later use.

2) Use FSLeyes volume viewer to find centrepoint coordinates for vessel locations of interest, and distal coordinates along the length of the vessel, enter these into a spreadsheet (i.e. 3 columns for centrepoint x,y,z, 3 columns for distal x,y,x, then repeat for all vessels)

3) Use extractSlice_4Dv2D_v3.m to extract all the 2D slices from the subject's 4D volume, for however many vessels are of interest. A radius of 150 worked for me but this may need changing

4) Draw the masks/ROIs on all 2D images, save these and use move_masks.m to move the masks to the ROIs folder

5) processFlowData_4Dv2D_AM.m is a script that runs setFlowParameters_4Dv2D_AM.m and then carries out pipeline_extractFlow_4Dv2D_AM to produce flow graphs for each vessel in the slice, for each velocity direction. Dot product calculation is used to take the velocities of a pixel in all 3 directions and calculate the velocity along the vessel length/vector

6) Coment the extractFlow line and uncomment the pipeline_finalFlowCalculations line to produce flow waveforms and flow/pulsatility variables (in a spreadsheet) for the subjects of interest
