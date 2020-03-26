% this code takes the vessel of interest's name, slice central coordinates and
% coordinated further along the vessel, and calculates the vector along the vessel
% (i.e. normal the slice plane). It then extracts this slice for all 3 velocities and
% magnitude (from the 4D slab), at every timeframe - for 2D flow analysis later
% advised to use ImageJ's volumetric viewer for finding coordinates

% this is version 3, which now produces a series of images across
% time with combined velocities across all 3 directions

clear all
close all

setFlowParameters_4Dv2D_AM;
procRoot='/DSTORE/BRICIA/amorgan_PhD/';
addpath([procRoot 'MATLAB/NIfTI_tools']);
addpath('/DSTORE/BRICIA/resources/Temp_MS/Matlab/MJT_Matlab_Path/utilities');
spm('defaults', 'FMRI'); spm_jobman('initcfg');

image_dir = [procRoot '4DFlowProject/4Dv2Danalysis/'];
SubjectSpreadsheet = [procRoot '4DFlowProject/SubjectDatabase.xlsx'];
Subjectdata = readtable(SubjectSpreadsheet);

VectorSpreadsheet = [image_dir 'Vessel_centres+normals.xlsx'];
VesselLocations = readtable(VectorSpreadsheet);

vessels = { 'RMCA' 'LMCA' 'RACA' 'LACA' 'RPCA' 'LPCA' 'SSS' 'StS' 'RTS' 'LTS' 'RICA' 'LICA' 'BA'};
radius = 150;
velDirections = {'VelF-H', 'VelR-L', 'VelA-P', 'Mag'} ;
velDirections2 = {'F-H', 'R-L', 'A-P', 'Mag'}; % this is because the 3d nifti file saved for each 40 slice slab has these terms in the name,
%so are required when referring to them with spm_vol




for iSubject=[11]
    
    SubjectNum=iSubject;  
    flowDir = ([image_dir char(Subjectdata{SubjectNum,1}) '/FLOWIMAGES']); % directory containing niftis
    nSlices = Subjectdata{SubjectNum,7};
    nTimeFrames = Subjectdata{SubjectNum,28};

    for iFlowScan = [9 10]  % '1.RMCA 2.LMCA 3.RACA 4.LACA 5.RPCA 6.LPCA 7.SSS 8.StS 9.RTS 10.LTS 11.RICA 12.LICA 13.BA
        
        vesselDir = ([flowDir '/' vessels{iFlowScan}]);
        if ~exist(vesselDir,'dir'); mkdir(vesselDir); else delete([vesselDir '/*']); end % create directory for output or delete files in existing directory
               
        centre = [str2double(VesselLocations{SubjectNum+1, iFlowScan*7 -4}) (str2double(VesselLocations{SubjectNum+1, iFlowScan*7 -5}))  str2double(VesselLocations{SubjectNum+1, iFlowScan*7 -3})];
        dist = [str2double(VesselLocations{iSubject+1, iFlowScan*7 -1}) (str2double(VesselLocations{iSubject+1, iFlowScan*7 -2})) str2double(VesselLocations{iSubject+1, iFlowScan*7})];        
        vec = dist - centre;
        
        
        velMeans = cell(1,3); velAbsMeans = cell(1,3); velMax = cell(1,3); velAbsMax = cell(1,3); % for producing merged images
        
        for scanNo = 1:4 % the numbers of the 4Dflow scans (velx, vely, velz, mag)
            disp(['Processing subject ' char(Subjectdata{SubjectNum,1}) ', ' vessels{iFlowScan} ', ' velDirections{scanNo} ]);
            
            for TF = 1:nTimeFrames % for each timeframe
                V_vel=spm_vol([flowDir '/' velDirections{scanNo} '/vol' num2str(TF-1,'%04.f') '.nii']); % load F-H velocity data for this timeframe
                [displayImage,xyz]=spm_read_vols(V_vel); % loads phase-contrast images
                displayImage2=flipud(rot90(displayImage));
                displayImage2=squeeze(displayImage2(:,:,:));
                
                % extract slice specified by centre coordinates and plane normal at top of file
                [slice, sliceInd,subX,subY,subZ] = extractSlice(displayImage2,centre(1),centre(2),centre(3),vec(1),vec(2),vec(3),radius);
                surf(subX,subY,subZ,slice,'FaceColor','texturemap','EdgeColor','none');
                xlabel('x'); ylabel('y'); zlabel('z');
                drawnow;
                
                directionDir = ([vesselDir '/' velDirections2{scanNo}]); % within the vessel folder, create vel/mag folders
                if ~exist(directionDir,'dir'); mkdir(directionDir); end
                
                nii_img = make_nii(slice);
                save_nii(nii_img, [directionDir '/slice_TF' num2str(TF,'%02.f') '.nii']); % save nifti for this timeframe
                
            end
            
            
            %% Merge nii files to create "4D" files
            disp(['Merging nii files for ' char(Subjectdata{SubjectNum,1}) ', ' vessels{iFlowScan} ]);
            
            if scanNo == 1 || scanNo == 2 || scanNo ==3 % velocity images
                % merge niftis
                spm_file_merge(sort(getMultipleFilePaths([directionDir '/slice_TF*.nii'])),[vessels{iFlowScan} '_' velDirections2{scanNo} '_Vel.nii'],0);
                delete([directionDir '/slice_TF*.nii']); % delete separate niftis
                
                % read in slice across time
                V_vel=spm_vol([directionDir '/' vessels{iFlowScan} '_' velDirections2{scanNo} '_Vel.nii']);
                [vel,xyz]=spm_read_vols(V_vel);
                
                % Create derived images
                fileSuffix={'Vel_MeanAbs' 'Vel_Mean' 'Vel_Max' 'Vel_MaxAbs'};
                fileData={ mean(abs(vel),4) mean(vel,4) max(vel,[],4) max(abs(vel),[],4)}; % takes the absolute or normal velocity values across time (4th dimension) and then takes means/max of these
                fileHeader={ V_vel V_vel V_vel V_vel};
                
                for iFile=1:size(fileSuffix,2)
                    V=fileHeader{iFile}(1); V.fname=[directionDir '/' vessels{iFlowScan} '_' velDirections2{scanNo} '_' fileSuffix{iFile} '.nii'];
                    spm_write_vol(V,fileData{iFile});
                    if iFile == 1; velAbsMeans{scanNo} = mean(abs(vel),4);end
                    if iFile == 2; velMeans{scanNo} = mean(vel,4);end
                    if iFile == 3; velMax{scanNo} = max(vel,[],4);end
                    if iFile == 4; velAbsMax{scanNo} = max(abs(vel),[],4);end
                end
                
                
             %%%%%%%%%%% combined velocity images %%%%%%%%%%%
                if scanNo == 3
                    fileSuffix = {'Combined_Vel_MeanAbs' 'Combined_Vel_Mean' 'Combined_Vel_Max' 'Combined_Vel_MaxAbs'};
                    a = sqrt(velAbsMeans{1}.^2 + velAbsMeans{2}.^2 + velAbsMeans{3}.^2);
                    b = sqrt(velMeans{1}.^2 + velMeans{2}.^2 + velMeans{3}.^2);
                    c = sqrt(velMax{1}.^2 + velMax{2}.^2 + velMax{3}.^2);
                    d = sqrt(velAbsMax{1}.^2 + velAbsMax{1}.^2 + velAbsMax{1}.^2);
                    fileData = {a b c d};
                    fileHeader = V_vel;
                    for iFile=1:size(fileSuffix,2)
                        V=fileHeader(1); V.fname=[vesselDir '/' vessels{iFlowScan} '_' fileSuffix{iFile} '.nii'];
                        spm_write_vol(V,fileData{iFile});
                    end
                    
                    % create a series of velocity images that combines all 3 directions
                    V_vel1 = spm_vol([vesselDir '/' velDirections2{1} '/' vessels{iFlowScan} '_'  velDirections2{1} '_Vel.nii']);
                    [vel1,~]=spm_read_vols(V_vel1);
                    vel1 = squeeze(vel1);
                    V_vel2 = spm_vol([vesselDir '/' velDirections2{2} '/' vessels{iFlowScan} '_'  velDirections2{2} '_Vel.nii']);
                    [vel2,~]=spm_read_vols(V_vel2);
                    vel2 = squeeze(vel2);
                    V_vel3 = spm_vol([vesselDir '/' velDirections2{3} '/' vessels{iFlowScan} '_'  velDirections2{3} '_Vel.nii']);
                    [vel3,~]=spm_read_vols(V_vel3);
                    vel3 = squeeze(vel3);
             
                    vel4 =sqrt(vel1.^2 + vel2.^2 +vel3.^2);
                    fileHeader = V_vel1;
                    V=fileHeader; % V.fname=[vesselDir '/' vessels{iFlowScan} '_Combined_Vel.nii'];
                    
                    for n=1:25
                        V(n,1).fname=[vesselDir '/' vessels{iFlowScan} '_Combined_Vel.nii'];
                        spm_write_vol(V(n,1),vel4(:,:,n));
                    end
                    
                end
                 
                
            else % magnitude images
                spm_file_merge(sort(getMultipleFilePaths([directionDir '/slice_TF*.nii'])),[vessels{iFlowScan} '_Mag.nii'],0);
                delete([directionDir '/slice_TF*.nii']);
                
                V_mag=spm_vol([directionDir '/' vessels{iFlowScan} '_Mag.nii']);
                [mag,xyz]=spm_read_vols(V_mag);
                
                % Create derived images
                fileSuffix={'Mag_Mean' 'Mag_Max'};
                fileData={ mean(mag,4) max(mag,[],4) };
                fileHeader={ V_mag V_mag };
                
                
                for iFile=1:size(fileSuffix,2)
                    V=fileHeader{iFile}(1); V.fname=[directionDir '/' vessels{iFlowScan} '_' fileSuffix{iFile} '.nii'];
                    spm_write_vol(V,fileData{iFile});
                end
            end
            
        end
    end
    
    
end



