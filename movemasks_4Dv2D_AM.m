% simple code to copy the masks drawn around vessel lumens in FSLeyes to a masks directory named ROIs
% make sure to create the ROIs folder manually in the 4D analysis directory (or have the code create it just once)

setFlowParameters_4Dv2D_AM;
load('allFlowOptions');
load(%point to wherever scaninfo.mat is)

spm('defaults', 'FMRI'); spm_jobman('initcfg');

for iSubject=[11]
    for iFlowScan=7:10%1:13% 1.RMCA 2.LMCA 3.RACA 4.LACA 5.RPCA 6.LPCA 7.SSS 8.StS 9.RTS 10.LTS 11.RICA 12.LICA 13.BA
        
        opts=allOpts{iSubject,iFlowScan}; %load options
        
        if iFlowScan==1
            system(['mkdir ' opts.masksDir]);
        else
            %fprintf('Directory already exists\n')
        end
        
        if isempty(opts); continue; end
            system(['cp ' opts.FlowImagesDir '/*.nii.gz ' opts.masksDir]); %copy paste compressed nifti images to ROI directory

    end
    
    system(['for i in $(ls ' opts.masksDir '*.nii.gz);do fslchfiletype NIFTI ${i};done']); % convert compressed niftis to regular niftis
        gunzip([opts.masksDir '/*.nii.gz']);
    
end
