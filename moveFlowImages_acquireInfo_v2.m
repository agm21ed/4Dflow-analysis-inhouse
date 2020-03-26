dicomDir = '/DSTORE/BRICIA/E182012_Brain_4D_Flow_Res/'; %point to 4D flow dicoms
targetDir = '/DSTORE/BRICIA/amorgan_PhD/4DFlowProject/4Dv2Danalysis/'; % new directory for niftis separated into timeframes
addpath([targetDir]);


SubjectSpreadsheet ='/DSTORE/BRICIA/amorgan_PhD/4DFlowProject/SubjectDatabase.xlsx';  % subject data
Subjectdata = readtable(SubjectSpreadsheet);

code = '120005_';


for iSubject=11
    
    SubjectNum=iSubject;
    disp(['Acquiring timeframe slabs for HV0' num2str(iSubject,'%02.f')]);
    flowDir = ([targetDir char(Subjectdata{SubjectNum,1}) '/FLOWIMAGES']);
        if ~exist(flowDir,'dir'); mkdir(flowDir); else delete([flowDir '/*']); end % create directory for output or delete files in existing directory
    
    a=Subjectdata{SubjectNum,15}; b=Subjectdata{SubjectNum,16}; c=Subjectdata{SubjectNum,17}; d=Subjectdata{SubjectNum,18};
    date = num2str(Subjectdata{SubjectNum,27});
    orientation = char(Subjectdata{SubjectNum,23});
    resolution = char(Subjectdata{SubjectNum,24});
    segments = char(Subjectdata{SubjectNum,25});
    grappa = char(Subjectdata{SubjectNum,26});
    nSlices = Subjectdata{SubjectNum,7};
    nTimeFrames = Subjectdata{SubjectNum,28};
    
    
    for scanNo = [a b c d] % b c d] % the numbers of the 4Dflow scans (velx, vely, velz, mag)
        if scanNo == a % vel-FH
            disp('Direction Foot-Head');
            velDir = ([flowDir '/VelF-H/']);
            system(['mkdir ' velDir]); % make a new directory Vel-FH
            
            cd([dicomDir char(Subjectdata{SubjectNum,1}) '/'  num2str(scanNo) '_4DFlowWIP_' grappa '_' orientation '_' resolution segments '_P']); % change to velocity directory
            system('dcm2niix .'); % convert dicoms 4D nifti slab
            system(['fslsplit __4DFlowWIP_' grappa '_' orientation '_' resolution segments '_' date code num2str(a) '_ph.nii -t']) % then use fslsplit to split the slab into N timeframes
            system('gunzip *.nii.gz'); % convert slabs to niftis
            system(['cp ', '*.nii', ' ', velDir]);  % copy all the slabs into the appropriate directory in amorgan_PhD folder
            delete('vol*.nii'); % remove slabs from study folder
            
            cd(targetDir)
%             
            %
            % retrieve the relevant header information for later use
            
            dicom = ([dicomDir char(Subjectdata{SubjectNum,1}) '/' num2str(scanNo) '_4DFlowWIP_' grappa '_' orientation '_' resolution segments '_P/1.dcm']);
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
            
            cd([dicomDir char(Subjectdata{SubjectNum,1}) '/'  num2str(scanNo) '_4DFlowWIP_' grappa '_' orientation '_' resolution segments '_P']); % change to velocity directory
            system('dcm2niix .'); % convert dicoms 4D nifti slab
            system(['fslsplit __4DFlowWIP_' grappa '_' orientation '_' resolution segments '_' date code num2str(b) '_ph.nii -t']) % then use fslsplit to split the slab into N timeframes
            system('gunzip *.nii.gz'); % convert slabs to niftis
            system(['cp ', '*.nii', ' ', velDir]);  % copy all the slabs into the appropriate directory in amorgan_PhD folder
            delete('vol*.nii'); % remove slabs from study folder
            
            cd(targetDir)
            %
            
            % retrieve the relevant header information for later use
            dicom = ([dicomDir char(Subjectdata{SubjectNum,1}) '/' num2str(scanNo) '_4DFlowWIP_' grappa '_' orientation '_' resolution segments '_P/1.dcm']);
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
            
            cd([dicomDir char(Subjectdata{SubjectNum,1}) '/'  num2str(scanNo) '_4DFlowWIP_' grappa '_' orientation '_' resolution segments '_P']); % change to velocity directory
            system('dcm2niix .'); % convert dicoms 4D nifti slab
            system(['fslsplit __4DFlowWIP_' grappa '_' orientation '_' resolution segments '_' date code num2str(c) '_ph.nii -t']) % then use fslsplit to split the slab into N timeframes
            system('gunzip *.nii.gz'); % convert slabs to niftis
            system(['cp ', '*.nii', ' ', velDir]);  % copy all the slabs into the appropriate directory in amorgan_PhD folder
            delete('vol*.nii'); % remove slabs from study folder
            
            cd(targetDir)
            
            
            %             retrieve the relevant header information for later use
            dicom = ([dicomDir char(Subjectdata{SubjectNum,1}) '/'  num2str(scanNo) '_4DFlowWIP_' grappa '_' orientation '_' resolution segments '_P/1.dcm']);
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
            
            cd([dicomDir char(Subjectdata{SubjectNum,1}) '/'  num2str(scanNo) '_4DFlowWIP_' grappa '_' orientation '_' resolution segments]); % change to mag directory
            system('dcm2niix .'); % convert dicoms 4D nifti slab
            system(['fslsplit __4DFlowWIP_' grappa '_' orientation '_' resolution segments '_' date code num2str(d) '.nii -t']) % then use fslsplit to split the slab into N timeframes
            system('gunzip *.nii.gz'); % convert slabs to niftis
            system(['cp ', '*.nii', ' ', magDir]);  % copy all the slabs into the appropriate directory in amorgan_PhD folder
            delete('vol*.nii'); % remove slabs from study folder
            
            cd(targetDir)
        end
    end
end