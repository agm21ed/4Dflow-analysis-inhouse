%this code is for analysing 2D PC-MRI images from the 4DFlow scans

%%first let's find our images
% procRoot='U:\Datastore\CMVM\scs\groups\BRICIA\amorgan_PhD\';
procRoot='/DSTORE/BRICIA/amorgan_PhD/4DFlowProject/';
% procRoot='U:\Datastore\CMVM\scs\groups\BRICIA\amorgan_PhD\4DFlowProject\';
addpath([procRoot(1:end-14) 'MATLAB/spm12']);
addpath([procRoot(1:end-14) 'MATLAB/dicom2nifti']);

SubjectSpreadsheet ='/DSTORE/BRICIA/amorgan_PhD/4DFlowProject/SubjectDatabase.xlsx';
Subjectdata = readtable(SubjectSpreadsheet);


% addpath('/PROJ/CVR_STUDIES/SOFTWARE/Flow1.3');
addpath('/DSTORE/BRICIA/resources/Temp_MS/Matlab/MJT_Matlab_Path/utilities');
% addpath('U:\Datastore\CMVM\scs\groups\BRICIA\resources/Temp_MS/Matlab/MJT_Matlab_Path/utilities');
% addpath( 'DSTORE/BRICIA/resources/Temp_MS/Matlab/MJT_Matlab_Path/random_downloaded_functions');
load([procRoot '4Dv2Danalysis/scanInfo.mat']);

allOpts=cell(scanInfo.N,13); %subject,flow scan - unsure of this part

directionNames = {'FH' 'LR' 'AP'};

maskNames= {'RMCA' 'LMCA' 'RACA' 'LACA' 'RPCA' 'LPCA' 'SSS' 'StS' 'RTS' 'LTS' 'RICA' 'LICA' 'BA'};

% flowSignCorrection_FH = [-1 -1 1 1 1 1 -1 -1 -1 -1 1 1 1];
% flowSignCorrection_RL = [1 -1 1 1 1 -1 1 1 1 1 -1 1 -1 1];
% flowSignCorrection_AP = [-1 -1 -1 -1 1 1 1 1 -1 -1 -1 1 1];

BGMaskNames= { 'RMCA_BG' 'LMCA_BG' 'RACA_BG' 'LACA_BG' 'RPCA_BG' 'LPCA_BG' 'SSS_BG' 'StS_BG' 'RTS_BG' 'LTS_BG' 'RICA_BG' 'LICA_BG' 'BA_BG'} ;

displayImages={ 'Mag/RMCA_Mag_Max' 'Mag/LMCA_Mag_Max' 'Mag/RACA_Mag_Max' 'Mag/LACA_Mag_Max' 'Mag/RPCA_Mag_Max' 'Mag/LPCA_Mag_Max'...
    'Mag/SSS_Mag_Max' 'Mag/StS_Mag_Max' 'Mag/RTS_Mag_Max' 'Mag/LTS_Mag_Max' 'Mag/RICA_Mag_Max' 'Mag/LICA_Mag_Max' 'Mag/BA_Mag_Max' }; %work out which display images are best for each scan


for iSubject=[1:5,7:scanInfo.N]
    
    for iFlowScan=1:13
        
        opts.HVNumberStr=scanInfo.HVNumberStr{iSubject}; %need scanInfo file for this study, HV = healthy volunteer

%         opts.scanName=scanNames{iFlowScan};
        opts.displayImage=displayImages{iFlowScan};
        opts.scanner='BRIC2';
        
        %data locations
        opts.FlowImagesDir=[procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject} '/FLOWIMAGES/' maskNames{iFlowScan}];
        opts.SubjectDir= [procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject} '/FLOWIMAGES/'];
        opts.masksDir=[procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject} '/ROIs/'];
        opts.SubjectDir2 = [procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject}];
        opts.dicomDir = ['/DSTORE/BRICIA/E182012_Brain_4D_Flow_Res/' scanInfo.HVNumberStr{iSubject} '/'];
        opts.noiseMasksDir = [procRoot '4Dv2Danalysis/noiseMasks/'];
        opts.noiseQuantDir = [procRoot 'Noise quantification/'];
        
        opts.maskNames = maskNames{iFlowScan};
        opts.BGMaskName = BGMaskNames{iFlowScan};
%         opts.flowSignCorrection_FH = flowSignCorrection_FH(iFlowScan);
%         opts.flowSignCorrection_RL = flowSignCorrection_RL(iFlowScan);
%         opts.flowSignCorrection_AP = flowSignCorrection_AP(iFlowScan);
        opts.aliasCorrection = 0;
        opts.flowDataDir1 =[procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject} '/FlowData_Basic/' maskNames{iFlowScan}];
        opts.flowDataDir2 =[procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject} '/FlowData_DotProduct/' maskNames{iFlowScan}];
        opts.flowDataDir3 =[procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject} '/FlowData_PartialVolumeCorrection/' maskNames{iFlowScan}];
        
        opts.overwriteImages=0;
        
        %% exceptions 
% reminder of vessel order: {1'RMCA' 2'LMCA' 3'RACA' 4'LACA' 5'RPCA' 6'LPCA' 7'SSS' 8'StS' 9'RTS' 10'LTS' 11'RICA' 12'LICA' 13'BA'}

%         switch opts.HVNumberStr
%             case 'HV002'
%                 flowSignCorrection_AP = [-1 -1 -1 -1 -1 1 -1 1 -1 -1 1 -1 -1];
%                 flowSignCorrection_FH = [1 -1 -1 1 1 1 -1 -1 -1 -1 1 1 1];
%                 flowSignCorrection_LR = [1 -1 1 1 1 -1 -1 -1 1 1 -1 -1 1 -1];end 
%             case 'HV004'
%                 flowSignCorrection_AP = [1 1 1 1 -1 1 1 1 1 1 1 1 -1 1] ;
%                 flowSignCorrection_FH = [-1 -1 -1 1 -1 1 1 1 1 1 1 1];
%                 flowSignCorrection_LR = [-1 -1 -1 -1 -1 -1 -1 1 -1 1 1 -1 1];end     

        
        allOpts{iSubject,iFlowScan}=opts;
    end
end

save('allFlowOptions','allOpts');