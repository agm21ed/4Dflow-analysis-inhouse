%this code sets all the necessary parameters and directory names for later 4D flow analysis

%%first let's find our images

procRoot= % location of project directory
addpath(%location of spm matlab tools;
addpath(% location of nifti matlab tools

SubjectSpreadsheet =[% location of spreadsheet containing subject IDs and scan info
Subjectdata = readtable(SubjectSpreadsheet);


load(% location of scanInfo.mat);

allOpts=cell(scanInfo.N,13); %subject,flow scan 

directionNames = {'FH' 'LR' 'AP'};
maskNames= {'RMCA' 'LMCA' 'RACA' 'LACA' 'RPCA' 'LPCA' 'SSS' 'StS' 'RTS' 'LTS' 'RICA' 'LICA' 'BA'};
BGMaskNames= { 'RMCA_BG' 'LMCA_BG' 'RACA_BG' 'LACA_BG' 'RPCA_BG' 'LPCA_BG' 'SSS_BG' 'StS_BG' 'RTS_BG' 'LTS_BG' 'RICA_BG' 'LICA_BG' 'BA_BG'} ;
displayImages={ 'Mag/RMCA_Mag_Max' 'Mag/LMCA_Mag_Max' 'Mag/RACA_Mag_Max' 'Mag/LACA_Mag_Max' 'Mag/RPCA_Mag_Max' 'Mag/LPCA_Mag_Max'...
    'Mag/SSS_Mag_Max' 'Mag/StS_Mag_Max' 'Mag/RTS_Mag_Max' 'Mag/LTS_Mag_Max' 'Mag/RICA_Mag_Max' 'Mag/LICA_Mag_Max' 'Mag/BA_Mag_Max' }; %  display images that are best for each scan


for iSubject=[1:scanInfo.N]
    for iSes = 1:2
        
        for iFlowScan=1:13
            
            opts.HVNumberStr=scanInfo.HVNumberStr{iSubject}; %need scanInfo file for this study, HV = healthy volunteer
            
            %         opts.scanName=scanNames{iFlowScan};
            opts.displayImage=displayImages{iFlowScan};
            opts.scanner='BRIC2';
            
            %data locations
            opts.FlowImagesDir=[procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject} '/FLOWIMAGES/v' num2str(iSes) '/' maskNames{iFlowScan} '/'];
            opts.SubjectDir= [procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject} '/FLOWIMAGES/'];
            opts.masksDir1=[procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject} '/ROIs/'];
            opts.masksDir2=[opts.masksDir1 'v' num2str(iSes) '/'];            
            opts.SubjectDir2 = [procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject}];
            opts.FlowImagesDir2 = [opts.SubjectDir 'v' num2str(iSes) '/'];
            opts.dicomDir = ['/DSTORE/BRICIA/E202194_Brain_4D_Flow_2_RES/RIE' scanInfo.HVNumberStr{iSubject} '/'];

            
            opts.maskNames = maskNames{iFlowScan};
            opts.BGMaskName = BGMaskNames{iFlowScan};
            opts.aliasCorrection = 0;
            opts.flowDataDir2 =[procRoot '4Dv2Danalysis/' scanInfo.HVNumberStr{iSubject} '/FlowData_DotProduct/v' num2str(iSes) '/' maskNames{iFlowScan}];
                   
            opts.overwriteImages=0;
            
            
            allOpts{iSubject,iSes,iFlowScan}=opts;
        end
    end
end

save('allFlowOptions','allOpts');
