clear; close all;
% clear classes

setFlowParameters_4Dv2D_AM; %run the set parameters matlab file
load('allFlowOptions');
load([procRoot '4Dv2Danalysis/scanInfo.mat']); %load in the scan info matlab file from wherever you put it


spm('defaults', 'FMRI'); spm_jobman('initcfg');

% addpath([procRoot '2Danalysis']) %point to where other matlab files are for processing

for iSubject=[11]
    
    
%% Comment/uncomment this section to use region-growing algorithm to draw masks
%     for iFlowScan=14 %  1.RMCA 2.LMCA 3.RACA 4.LACA 5.RPCA 6.LPCA 7.SSS 8.StS 9.RTS 10.LTS 11.RICA 12.LICA 13.BA
% %         regionGrowingFunction(iSubject,iFlowScan,0.6,0.4)     %iSubject, iFlowScan, lumen_threshold (higher), wall_threshold (lower)
% %         regionGrowingFunction_PROTOTYPE(iSubject,iFlowScan)   
%             regionGrowingFunction_INTERPOLATED(iSubject,iFlowScan)
%     end
    

%This section uses the ROI and BG masks to extract flow using different methods
    for iFlowScan=[9:10] %  1.RMCA 2.LMCA 3.RACA 4.LACA 5.RPCA 6.LPCA 7.SSS 8.StS 9.RTS 10.LTS 11.RICA 12.LICA 13.BA
%                      
            opts=allOpts{iSubject,iFlowScan}; %load options
            if isempty(opts); continue; end
            
            disp(['Processing subject ' opts.HVNumberStr ', ' opts.maskNames]);
%      
            pipeline_extractFlowBasic_4Dv2D(opts,iSubject,iFlowScan); % calculates flow by combining velocity components within each mask pixel
            pipeline_extractFlowDotProduct_4Dv2D(opts,iSubject,iFlowScan); % uses dot product of vector along vessel and velocity vector
%             pipeline_extractFlowPVC_4Dv2D(opts,iSubject,iFlowScan); % same as above but uses Buillot's PVC algorithm to calculate flow

        
            
% %    comment/uncomment this code to compare the three flow extraction methods above
%     
% %            methodComparisonPlot_function(opts);

%             
    end

%     %% this section is for tidying up and plotting the flow/pulsatility results
%     
%     pipeline_finalFlowCalculationsBasic(iSubject); % for sum of square velocity flow with accurate pixel areas
%     pipeline_finalFlowCalculationsDotProduct(iSubject); % for dot product flow with accurate pixel areas
%     pipeline_finalFlowCalculationsPVC(iSubject); % for partial volume corrected flow with accurate pixel areas
%     
%     %% noise quantification
%     opts=allOpts{iSubject};
%     noise_quantification_v2(opts,iSubject);
    
    
end