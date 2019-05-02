clear; %close all;

setFlowParameters_4Dv2D_AM; %run the set parameters matlab file
load('allFlowOptions');
load([procRoot '4Dv2Danalysis/scanInfo.mat']); %load in the scan info matlab file from wherever you put it


spm('defaults', 'FMRI'); spm_jobman('initcfg');

% addpath([procRoot '2Danalysis']) %point to where other matlab files are for processing

for iSubject=1%:7%:scanInfo.N
    
    for iFlowScan=1 % i.e. ACAs, RMCA, veins
                     
            opts=allOpts{iSubject,iSes,iFlowScan}; %load options
            if isempty(opts); continue; end
            opts.scanner='BRIC2';
            
            disp(['Processing subject ' opts.HVNumberStr ', ' opts.scanName]);
            
            %         pipeline_convertFlowDicoms_BRIC2_2DPC_AM(opts); %convert DICOMs to Niftis
            pipeline_extractFlow_4Dv2D_AM(opts); %extract flow data       
    end
    
end