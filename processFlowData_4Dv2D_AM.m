% this is the main script where the 4D low data can be extracted from the previously-drawn ROIs

clear; close all;
% clear classes

setFlowParameters_4Dv2D_AM; %run the set parameters matlab file
load('allFlowOptions');
load(% load scanInfo.mat from where it is); 

spm('defaults', 'FMRI'); spm_jobman('initcfg');



for iSubject=[1]
    for iSes = 1
       
        
        %This section uses the ROI and BG masks to extract flow using different methods
        for iFlowScan=1%[1:6,11:13] % 1.RMCA 2.LMCA 3.RACA 4.LACA 5.RPCA 6.LPCA 7.SSS 8.StS 9.RTS 10.LTS 11.RICA 12.LICA 13.BA
            
            opts=allOpts{iSubject,iSes,iFlowScan}; %load options
            if isempty(opts); disp('error: opts is empty'); continue; end

            disp(['Processing subject ' opts.HVNumberStr ', visit ' num2str(iSes) ', '  opts.maskNames]);
%             %
                % FLOW EXTRACTION OF HAND-DRAWN MASKS
             pipeline_extractFlowDotProduct_4Dv2D(opts,iSubject,iSes,iFlowScan); % uses dot product of vector along vessel and velocity vector

        end
        
        %     %% this section is for tidying up and plotting the flow/pulsatility results
        %   uncomment when you want to run the final phase of processing which tidies up the data into nicer graphs and saves the flow/pulsatility variables into a spreadsheet
        % NOTE, best to run this for all subjects in a separate loop once all the flow extraction has been completed
        
%             pipeline_finalFlowCalculationsDotProduct(iSubject,iSes); % for dot product flow with accurate pixel areas



    end
    
end
