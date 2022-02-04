% this code takes the 4D flow dicoms across all 3 directions, as well as the magnitude images, and moves them into the analyis folder
% also: the DICOMs are converted to niftis for and the scan header info is saved for late


dicomDir =  %point to 4D flow dicom directory
targetDir = % point to 4D analysis directory
addpath([targetDir]);


SubjectSpreadsheet = %location of spreadsheet containing DICOM directory IDs and subject IDsSubjectdata = readtable(SubjectSpreadsheet);

for iSubject=12
    for iSes = 1
        
        SubjectNum=iSubject;
        disp(['Acquiring timeframe slabs for HV0' num2str(iSubject,'%02.f')]);
        flowDir = ([targetDir char(Subjectdata{SubjectNum,1}) '/FLOWIMAGES/v' num2str(iSes)]);
        if ~exist(flowDir,'dir'); mkdir(flowDir); else delete([flowDir '/*']); end % create directory for output or delete files in existing directory
        
        if iSes == 1
            a=Subjectdata{SubjectNum,15}; b=Subjectdata{SubjectNum,16}; c=Subjectdata{SubjectNum,17}; d=Subjectdata{SubjectNum,18};
        elseif iSes == 2
            a=Subjectdata{SubjectNum,27}; b=Subjectdata{SubjectNum,28}; c=Subjectdata{SubjectNum,29}; d=Subjectdata{SubjectNum,30};            
        end
        date = num2str(Subjectdata{SubjectNum,23});
        nSlices = Subjectdata{SubjectNum,7};
        nTimeFrames = Subjectdata{SubjectNum,24};
        
        
        for scanNo = [a b c d] % b c d] % the numbers of the 4Dflow scans (velx, vely, velz, mag)
            if scanNo == a % vel-FH
                disp('Direction Foot-Head');
                velDir = ([flowDir '/VelF-H/']);
                system(['mkdir ' velDir]); % make a new directory Vel-FH
                
                cd([dicomDir char(Subjectdata{SubjectNum,1}) '/v' num2str(iSes) '/'  num2str(scanNo) '_4DFlow_WIP_p3_seg4_1mm_P']); % change to velocity directory
                system('dcm2niix .'); % convert dicoms 4D nifti slab
                filename = dir(['__4DFlow_WIP_p3_seg4_1mm_' date '*.nii']);
                system(['fslsplit ' filename(1).name ' -t']) % then use fslsplit to split the slab into N timeframes
                system('gunzip *.nii.gz'); % convert slabs to niftis
                system(['cp ', '*.nii', ' ', velDir]);  % copy all the slabs into the appropriate directory in amorgan_PhD folder
                delete('vol*.nii'); % remove slabs from study folder
                
                cd(targetDir)
                %
                %
                % retrieve the relevant header information for later use
                
                dicom = ([dicomDir char(Subjectdata{SubjectNum,1}) '/v' num2str(iSes) '/' num2str(scanNo) '_4DFlow_WIP_p3_seg4_1mm_P/1.dcm']);
                x=dicominfo(dicom);
                
                x=SiemensCsaParse(x);
                info.TR=x.RepetitionTime; info.TE=x.EchoTime; info.FA=x.FlipAngle; info.HR=round((1/(x.NominalInterval/1000))*60);
                info.TrigTime=x.(dicomlookup('0018','1060')); info.NomInterval=x.(dicomlookup('0018','1062'));
                info.PixelSpacing=x.(dicomlookup('0028','0030'));
                info.FlowVenc=x.csa.FlowVenc;
                info.orientation=x.ImageOrientationPatient;
                
                save([flowDir '/AcquisitionInfo_VelF-H'],'info');
                %
                
            elseif scanNo == b
                disp('Direction Right-Left');
                velDir = ([flowDir '/VelR-L/']);
                system(['mkdir ' velDir]); % make a new directory Vel-RL
                             
                cd([dicomDir char(Subjectdata{SubjectNum,1}) '/v' num2str(iSes) '/'  num2str(scanNo) '_4DFlow_WIP_p3_seg4_1mm_P']); % change to velocity directory
                system('dcm2niix .'); % convert dicoms 4D nifti slab
                filename = dir(['__4DFlow_WIP_p3_seg4_1mm_' date '*.nii']);
                system(['fslsplit ' filename(1).name ' -t']) % then use fslsplit to split the slab into N timeframes
                system('gunzip *.nii.gz'); % convert slabs to niftis
                system(['cp ', '*.nii', ' ', velDir]);  % copy all the slabs into the appropriate directory in amorgan_PhD folder
                delete('vol*.nii'); % remove slabs from study folder
                
                cd(targetDir)
                
                
                % retrieve the relevant header information for later use
                dicom = ([dicomDir char(Subjectdata{SubjectNum,1}) '/v' num2str(iSes) '/' num2str(scanNo) '_4DFlow_WIP_p3_seg4_1mm_P/1.dcm']);
                x=dicominfo(dicom);
                
                x=SiemensCsaParse(x);
                info.TR=x.RepetitionTime; info.TE=x.EchoTime; info.FA=x.FlipAngle; info.HR=round((1/(x.NominalInterval/1000))*60);
                info.TrigTime=x.(dicomlookup('0018','1060')); info.NomInterval=x.(dicomlookup('0018','1062'));
                info.PixelSpacing=x.(dicomlookup('0028','0030'));
                info.FlowVenc=x.csa.FlowVenc;
                info.orientation=x.ImageOrientationPatient;
                
                save([flowDir '/AcquisitionInfo_VelR-L'],'info');
                
            elseif scanNo == c
                disp('Direction Anterior-Posterior');
                velDir = ([flowDir '/VelA-P/']);
                system(['mkdir ' velDir]); % make a new directory Vel-AP
                
                
                cd([dicomDir char(Subjectdata{SubjectNum,1}) '/v' num2str(iSes) '/'  num2str(scanNo) '_4DFlow_WIP_p3_seg4_1mm_P']); % change to velocity directory
                system('dcm2niix .'); % convert dicoms 4D nifti slab
                filename = dir(['__4DFlow_WIP_p3_seg4_1mm_' date '*.nii']);
                system(['fslsplit ' filename(1).name ' -t']) % then use fslsplit to split the slab into N timeframes
                system('gunzip *.nii.gz'); % convert slabs to niftis
                system(['cp ', '*.nii', ' ', velDir]);  % copy all the slabs into the appropriate directory in amorgan_PhD folder
                delete('vol*.nii'); % remove slabs from study folder
                
                cd(targetDir)
                
                
                %             retrieve the relevant header information for later use
                dicom = ([dicomDir char(Subjectdata{SubjectNum,1}) '/v' num2str(iSes) '/' num2str(scanNo) '_4DFlow_WIP_p3_seg4_1mm_P/1.dcm']);
                x=dicominfo(dicom);
                
                x=SiemensCsaParse(x);
                info.TR=x.RepetitionTime; info.TE=x.EchoTime; info.FA=x.FlipAngle; info.HR=round((1/(x.NominalInterval/1000))*60);
                info.TrigTime=x.(dicomlookup('0018','1060')); info.NomInterval=x.(dicomlookup('0018','1062'));
                info.PixelSpacing=x.(dicomlookup('0028','0030'));
                info.FlowVenc=x.csa.FlowVenc;
                info.orientation=x.ImageOrientationPatient;
                
                save([flowDir '/AcquisitionInfo_VelA-P'],'info');
                
            else
                disp('Magnitude');
                magDir = ([flowDir '/Mag/']);
                system(['mkdir ' magDir]); % make a new directory Mag
                
                cd([dicomDir char(Subjectdata{SubjectNum,1}) '/v' num2str(iSes) '/'  num2str(scanNo) '_4DFlow_WIP_p3_seg4_1mm']); % change to mag directory
                system('dcm2niix .'); % convert dicoms 4D nifti slab
                filename = dir(['__4DFlow_WIP_p3_seg4_1mm_' date '*.nii']);
                system(['fslsplit ' filename(1).name ' -t']) % then use fslsplit to split the slab into N timeframes
                system('gunzip *.nii.gz'); % convert slabs to niftis
                system(['cp ', '*.nii', ' ', magDir]);  % copy all the slabs into the appropriate directory in amorgan_PhD folder
                delete('vol*.nii'); % remove slabs from study folder
                
                cd(targetDir)
            end
        end
    end
end
