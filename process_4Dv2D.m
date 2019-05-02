% this code takes the vessel of interest's name, slice central coordinates and
% normal vector, and extracts this slice for all 3 velocities and
% magnitude (from the 4D slab), at every timeframe - for 2D flow analysis later
% advised to use ImageJ's volumetric viewer for finding centrepoint and
% normal


clear all
close all

procRoot='/home/s1005152/W/amorgan_PhD/';
addpath([procRoot 'MATLAB/spm12']);
addpath([procRoot 'MATLAB/NIfTI_tools']);
spm('defaults', 'FMRI'); spm_jobman('initcfg');

image_dir = [procRoot '4DFlowProject/4Dv2Danalysis/'];
SubjectSpreadsheet = [procRoot '4DFlowProject/SubjectDatabase.xlsx'];  
Subjectdata = readtable(SubjectSpreadsheet);

VectorSpreadsheet = [image_dir 'Vessel_centres+normals.xlsx'];
VesselLocations = readtable(VectorSpreadsheet);

%% input values here
% x=90;y=100;z=20 ;
vessel = 'ACAs';
% centre = [y x z]; % centre, a point on the plane - [coronal (1:224), sagittal (1:180), axial (1:40)]
% vec = [0 0 1]; % normal - [coronal, sagittal, axial]; imagej: [y, z, x]


%%
radius = 150;
% nSlices = 40;
nTimeFrames = 25;
velDirections = {'VelF-H', 'VelL-R', 'VelA-P', 'Mag'} ;
velDirections2 = {'F-H', 'L-R', 'A-P', 'Mag'}; % this is because the 3d nifti file saved for each 40 slice slab has these terms in the name,
                                                %so are required when referring to them with spm_vol


for iSubject=2
    
    SubjectNum=iSubject;
    
    flowDir = ([image_dir char(Subjectdata{SubjectNum,1}) '/FLOWIMAGES']); % directory containing niftis
    vesselDir = ([flowDir '/' vessel]); % location of eventual oblique slices
    if ~exist(vesselDir,'dir'); mkdir(vesselDir); else delete([vesselDir '/*']); end % create directory for output or delete files in existing directory
    
    nSlices = Subjectdata{SubjectNum,7};
%     centre = [VesselLocations{SubjectNum,3} VesselLocations{SubjectNum,2} VesselLocations{SubjectNum,4}];
    centre = [83 91 29];
%     vec = [VesselLocations{SubjectNum,6} VesselLocations{SubjectNum,5} VesselLocations{SubjectNum,7}];
    vec = [-4 0 2];
    
    for scanNo = 1:4 % the numbers of the 4Dflow scans (velx, vely, velz, mag)
        
            for TF = 1:nTimeFrames % for each timeframe
                V_vel=spm_vol([flowDir '/' velDirections{scanNo} '/TF_' num2str(TF) '/slab_' velDirections2{scanNo} '_' num2str(TF) '.nii']); % load F-H velocity data for this timeframe
                [displayImage,xyz]=spm_read_vols(V_vel); % loads phase-contrast images
                displayImage2=(rot90(displayImage));
                displayImage2=squeeze(displayImage2(:,:,:));
                
                % extract slice specified by centre coordinates and plane normal at top of file
                [slice, sliceInd,subX,subY,subZ] = extractSlice(displayImage2,centre(1),centre(2),centre(3),vec(1),vec(2),vec(3),radius);
                surf(subX,subY,subZ,slice,'FaceColor','texturemap','EdgeColor','none');
                xlabel('x'); ylabel('y'); zlabel('z');
                drawnow;
                
                directionDir = ([vesselDir '/' velDirections2{scanNo}]); % within the vessel folder, create vel/mag folders
                if ~exist(directionDir,'dir'); mkdir(directionDir); end
                
                %adjust the slice so it can be read as normal 2D PC MRI slices are
                slice = imrotate(slice,180); slice = flipud(slice);
                nii_img = make_nii(slice);
                save_nii(nii_img, [directionDir '/slice_TF' num2str(TF) '.nii']); % save nifti for this timeframe  
                
            end
            
            
            %% Merge nii files to create "4D" files
            if scanNo == 1 || scanNo == 2 || scanNo ==3 % velocity images
                % merge niftis
                spm_file_merge(sort(getMultipleFilePaths([directionDir '/slice_TF*.nii'])),[vessel '_' velDirections2{scanNo} '_Vel.nii'],0);
                delete([directionDir '/slice_TF*.nii']); % delete separate niftis
                
                % read in slice across time
                V_vel=spm_vol([directionDir '/' vessel '_' velDirections2{scanNo} '_Vel.nii']);
                [vel,xyz]=spm_read_vols(V_vel);
                
                % Create derived images
                fileSuffix={'Vel_MeanAbs' 'Vel_Mean' 'Vel_Max' 'Vel_MaxAbs'};
                fileData={ mean(abs(vel),4) mean(vel,4) max(vel,[],4) max(abs(vel),[],4)};
                fileHeader={ V_vel V_vel V_vel V_vel};
                
                for iFile=1:size(fileSuffix,2)
                    V=fileHeader{iFile}(1); V.fname=[directionDir '/' vessel '_' velDirections2{scanNo} '_' fileSuffix{iFile} '.nii'];
                    spm_write_vol(V,fileData{iFile});
                end
                
            

            else % magnitude images
                spm_file_merge(sort(getMultipleFilePaths([directionDir '/slice_TF*.nii'])),[vessel '_Mag.nii'],0);
                delete([directionDir '/slice_TF*.nii']);
                
                V_mag=spm_vol([directionDir '/' vessel '_Mag.nii']);
                [mag,xyz]=spm_read_vols(V_mag);
                
                % Create derived images
                fileSuffix={'Mag_Mean' 'Mag_Max'};
                fileData={ mean(mag,4) max(mag,[],4) };
                fileHeader={ V_mag V_mag };
                
                
                for iFile=1:size(fileSuffix,2)
                    V=fileHeader{iFile}(1); V.fname=[directionDir '/' vessel '_' fileSuffix{iFile} '.nii'];
                    spm_write_vol(V,fileData{iFile});
                end
            end
            
          
    end
       scanNo
    
    
end
    
        
        
    