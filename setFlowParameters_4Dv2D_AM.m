%this code is for analysing 2D PC-MRI images from the 4DFlow scans

%%first let's find our images
% procRoot='U:\Datastore\CMVM\scs\groups\BRICIA\amorgan_PhD\';
procRoot='/home/s1005152/W/amorgan_PhD/4DFlowProject/';
addpath([procRoot(1:end-14) 'MATLAB/spm12']);
addpath([procRoot(1:end-14) 'MATLAB/dicom2nifti']);

% image_dir = [procRoot 'example_DICOMs\Vel_FH']; % here will be the location of the volunteer's DICOMS
%possibly need to load scanInfo.mat - need a getScanInfo file?

addpath('/ISIS/procB/CVR_STUDIES/SOFTWARE/Flow1.3');
addpath('/home/s1005152/W/resources/Temp_MS/Matlab/MJT_Matlab_Path/utilities');
addpath( '/ISIS/proc5/mthripp1/software/MJT_Matlab_Path/random_downloaded_functions');

load([procRoot '4Dv2Danalysis/scanInfo.mat']);

allOpts=cell(scanInfo.N,2,3); %subject,flow scan - unsure of this part

directionNames = {'FH' 'LR' 'AP'};

scanNames= {'ACAs' 'RMCA' 'Veins'}; % scan names

maskNames= { {'LACA' 'RACA' };... % vessel names
    { 'RMCA' } ;... 
    { 'SSS' 'StS' } };
flowSignCorrection_FH = { [1 1] ;... % ACAs flow upwards
    [-1];... % RMCA flows reverse plane
    [-1 -1] }; % SSS and StS flow reverse plane
flowSignCorrection_LR = { [1 1] ;... % ACAs have no real left or right flow
    [-1];... % RMCA flows reverse plane
    [-1 -1] }; % SSS and StS flow reverse plane
flowSignCorrection_AP = { [-1 -1] ;... % ACAs flow from back to front of head
    [-1];... % RMCA flows reverse plane
    [-1 -1] }; % SSS and StS flow reverse plane
BGMaskNames= { {'LACA_BG' 'RACA_BG'} ;...
    {'RMCA_BG' } ;...
    {'SSS_BG' 'StS_BG'} } ;
displayImages={'Mag/ACAs_Mag_Max' 'Mag/RMCA_Mag_Max' 'F-H/Veins_Vel_Mean' }; %work out which display images are best for each scan


for iSubject=1:scanInfo.N
        iSes=1;
    for iFlowScan=1:3
        
        opts.HVNumberStr=scanInfo.HVNumberStr{iSubject}; %need scanInfo file for this study, HV = healthy volunteer
        NMasks=size(maskNames{iFlowScan},2); %number of masks/ROIs for this scan
        
        %series-specific options
        switch iFlowScan
            case 1; opts.flowDicomDir=[scanInfo.ACADicomDir{iSubject,iSes}]; opts.SeriesNo=scanInfo.ACASeriesNo(iSubject,iSes); %retrieves DICOMS for each PC scan
            case 2; opts.flowDicomDir=[scanInfo.RMCADicomDir{iSubject,iSes}]; opts.SeriesNo=scanInfo.RMCASeriesNo(iSubject,iSes);
            case 3; opts.flowDicomDir=[scanInfo.VeDicomDir{iSubject,iSes}]; opts.SeriesNo=scanInfo.VeSeriesNo(iSubject,iSes);
        end
        if isnan(opts.SeriesNo); continue; end
        
        opts.scanName=scanNames{iFlowScan};
        opts.displayImage=displayImages{iFlowScan};
        opts.scanner='BRIC2';
        
        %data locations
        opts.FlowImagesDir=[procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject} '/FLOWIMAGES/' scanNames{iFlowScan}];
        opts.SubjectDir= [procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject} '/FLOWIMAGES/'];
        opts.masksDir=[procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject} '/ROIs/'];
        
        opts.maskNames = maskNames{iFlowScan};
        opts.BGMaskName = BGMaskNames{iFlowScan};
        opts.flowSignCorrection_FH = flowSignCorrection_FH{iFlowScan};
        opts.flowSignCorrection_LR = flowSignCorrection_LR{iFlowScan};
        opts.flowSignCorrection_AP = flowSignCorrection_AP{iFlowScan};
        opts.aliasCorrection = zeros(1,NMasks);
        opts.flowDataDir =[procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject} '/FlowData/' scanNames{iFlowScan}];
        
        opts.overwriteImages=0;
        
        %% exceptions
        switch opts.HVNumberStr
%             case 'HV_0041'
%                 if strcmp(scanName{iFlowScan},'Aqueduct'); opts.aliasCorrection=[zeros(25,1); 50*ones(1,1); 45*ones(2,1); 60*ones(2,1); 10*ones(3,1); 20*ones(2,1)];end
        end
        
        allOpts{iSubject,iSes,iFlowScan}=opts;
    end
end

save('allFlowOptions','allOpts');