
dicomDir = '/home/s1005152/W/E182012_Brain_4D_Flow_Res/'; %point to 4D flow dicoms
targetDir = '/home/s1005152/W/amorgan_PhD/4DFlowProject/4Dv2Danalysis/'; % new directory for niftis separated into timeframes
addpath([targetDir]);

SubjectSpreadsheet ='/home/s1005152/W/amorgan_PhD/4DFlowProject/SubjectDatabase.xlsx';  % subject data
Subjectdata = readtable(SubjectSpreadsheet);

nTimeFrames = 25;



for iSubject=2:3
    
    SubjectNum=iSubject;
    
    flowDir = ([targetDir char(Subjectdata{SubjectNum,1}) '/FLOWIMAGES']);
    if ~exist(flowDir,'dir'); mkdir(flowDir); else delete([flowDir '/*']); end; % create directory for output or delete files in existing directory
    
    a=Subjectdata{SubjectNum,15}; b=Subjectdata{SubjectNum,16}; c=Subjectdata{SubjectNum,17}; d=Subjectdata{SubjectNum,18};
    orientation = char(Subjectdata{SubjectNum,23});
    nSlices = Subjectdata{SubjectNum,7};
    
    
    for scanNo = [a,b,c,d] % the numbers of the 4Dflow scans (velx, vely, velz, mag)
        if scanNo == a % vel-FH
            
            velDir = ([flowDir '/VelF-H/']);
            system(['mkdir ' velDir]); % make a new directory Vel-FH
            
            for TF = 1:nTimeFrames % for each timeframe
                firstDicom = (TF-1)*nSlices + 1;
                lastDicom = TF*nSlices;
                TFdir = ([velDir 'TF_' num2str(TF)]);
                if ~exist(TFdir,'dir'); mkdir(TFdir); else delete([TFdir '/*']); end; % make a directory for each TF
                
                cd([dicomDir char(Subjectdata{SubjectNum,1}) '/' num2str(scanNo) '_4DFlowWIP_p3_' orientation '_1mm_P']); % change to velocity directory
                for n = firstDicom:lastDicom
                    fname = [num2str(n) '.dcm'];
                    system(['cp ', fname, ' ', TFdir]);  % copy all the slices for this timeframe into the appropriate directory
                end
                cd(velDir)
                system(['dcm2niix TF_' num2str(TF)]);
                cd (['TF_' num2str(TF)])
                movefile('series07.nii', ['slab_F-H_' num2str(TF) '.nii'])
            end
            
            % retrieve the relevant header information for later use
            
            dicom = ([dicomDir char(Subjectdata{SubjectNum,1}) '/' num2str(scanNo) '_4DFlowWIP_p3_'  orientation '_1mm_P/1.dcm']);
            x=dicominfo(dicom);
            
            x=SiemensCsaParse(x);
            info.TR=x.RepetitionTime; info.TE=x.EchoTime; info.FA=x.FlipAngle; info.HR=round((1/(x.NominalInterval/1000))*60);
            info.TrigTime=x.(dicomlookup('0018','1060')); info.NomInterval=x.(dicomlookup('0018','1062'));
            info.PixelSpacing=x.(dicomlookup('0028','0030'));
            info.FlowVenc=x.csa.FlowVenc;
            info.orientation=x.ImageOrientationPatient;
            
            save([flowDir '/AcquisitionInfo_VelF-H'],'info');
           
            
        elseif scanNo == b
            velDir = ([flowDir '/VelL-R/']);
            system(['mkdir ' velDir]); % make a new directory Vel-LR
            
            for TF = 1:nTimeFrames % for each timeframe
                firstDicom = (TF-1)*nSlices + 1;
                lastDicom = TF*nSlices;
                TFdir = ([velDir 'TF_' num2str(TF)]);
                if ~exist(TFdir,'dir'); mkdir(TFdir); else delete([TFdir '/*']); end;
                
                cd([dicomDir char(Subjectdata{SubjectNum,1}) '/' num2str(scanNo) '_4DFlowWIP_p3_'  orientation '_1mm_P']); % change to velocity directory
                for n = firstDicom:lastDicom
                    fname = [num2str(n) '.dcm'];
                    system(['cp ', fname, ' ', TFdir]);  % copy all the slices for this timeframe into the appropriate directory
                end
                
                
                cd(velDir)
                system(['dcm2niix TF_' num2str(TF)]);
                cd (['TF_' num2str(TF)])
                movefile('series07.nii', ['slab_L-R_' num2str(TF) '.nii'])
            end
            
            
            % retrieve the relevant header information for later use
            dicom = ([dicomDir char(Subjectdata{SubjectNum,1}) '/' num2str(scanNo) '_4DFlowWIP_p3_'  orientation '_1mm_P/1.dcm']);
            x=dicominfo(dicom);
            
            x=SiemensCsaParse(x);
            info.TR=x.RepetitionTime; info.TE=x.EchoTime; info.FA=x.FlipAngle; info.HR=round((1/(x.NominalInterval/1000))*60);
            info.TrigTime=x.(dicomlookup('0018','1060')); info.NomInterval=x.(dicomlookup('0018','1062'));
            info.PixelSpacing=x.(dicomlookup('0028','0030'));
            info.FlowVenc=x.csa.FlowVenc;
            info.orientation=x.ImageOrientationPatient;
            
            save([flowDir '/AcquisitionInfo_VelL-R'],'info');
            
        elseif scanNo == c
            velDir = ([flowDir '/VelA-P/']);
            system(['mkdir ' velDir]); % make a new directory Vel-AP
            
            for TF = 1:nTimeFrames % for each timeframe
                firstDicom = (TF-1)*nSlices + 1;
                lastDicom = TF*nSlices;
                TFdir = ([velDir 'TF_' num2str(TF)]);
                if ~exist(TFdir,'dir'); mkdir(TFdir); else delete([TFdir '/*']); end; % make a directory for each TF
                
                cd([dicomDir char(Subjectdata{SubjectNum,1}) '/' num2str(scanNo) '_4DFlowWIP_p3_'  orientation '_1mm_P']); % change to velocity directory
                for n = firstDicom:lastDicom
                    fname = [num2str(n) '.dcm'];
                    system(['cp ', fname, ' ', TFdir]);  % copy all the slices for this timeframe into the appropriate directory
                end
                cd(velDir)
                system(['dcm2niix TF_' num2str(TF)]);
                cd (['TF_' num2str(TF)])
                movefile('series07.nii', ['slab_A-P_' num2str(TF) '.nii'])
            end
            
            
            % retrieve the relevant header information for later use
            dicom = ([dicomDir char(Subjectdata{SubjectNum,1}) '/' num2str(scanNo) '_4DFlowWIP_p3_'  orientation '_1mm_P/1.dcm']);
            x=dicominfo(dicom);
            
            x=SiemensCsaParse(x);
            info.TR=x.RepetitionTime; info.TE=x.EchoTime; info.FA=x.FlipAngle; info.HR=round((1/(x.NominalInterval/1000))*60);
            info.TrigTime=x.(dicomlookup('0018','1060')); info.NomInterval=x.(dicomlookup('0018','1062'));
            info.PixelSpacing=x.(dicomlookup('0028','0030'));
            info.FlowVenc=x.csa.FlowVenc;
            info.orientation=x.ImageOrientationPatient;
            
            save([flowDir '/AcquisitionInfo_VelA-P'],'info');
            
        else
            magDir = ([flowDir '/Mag/']);
            system(['mkdir ' magDir]); % make a new directory Mag
            
            for TF = 1:nTimeFrames % for each timeframe
                firstDicom = (TF-1)*nSlices + 1;
                lastDicom = TF*nSlices;
                TFdir = ([magDir 'TF_' num2str(TF)]);
                if ~exist(TFdir,'dir'); mkdir(TFdir); else delete([TFdir '/*']); end; % make a directory for each TF
                
                cd([dicomDir char(Subjectdata{SubjectNum,1}) '/' num2str(scanNo) '_4DFlowWIP_p3_'  orientation '_1mm']); % change to velocity directory
                for n = firstDicom:lastDicom
                    fname = [num2str(n) '.dcm'];
                    system(['cp ', fname, ' ', TFdir]);  % copy all the slices for this timeframe into the appropriate directory
                end
                cd(magDir)
                system(['dcm2niix TF_' num2str(TF)]);
                cd (['TF_' num2str(TF)])
                movefile('series07.nii', ['slab_Mag_' num2str(TF) '.nii'])
            end
        end
    end
    
    
    
    
end