function pipeline_extractFlow_4Dv2D_AM(opts)

opengl('software');

close all force
set(0,'defaulttextinterpreter','none','DefaultAxesFontSize',10)

directionNames = {'FH' 'LR' 'AP'}; % named this way because variables can't have a '-' in
directionFileNames = {'F-H' 'L-R' 'A-P'}; % this exists to call the niftis with these names

%% create directory for output or delete files in existing directory
if ~exist(opts.flowDataDir,'dir'); mkdir(opts.flowDataDir); else delete([opts.flowDataDir '/*']); end

for iFlowDirection = 1:3
    %% Load 4D velocity image and acquisition parameters
    V_vel=spm_vol([opts.FlowImagesDir '/' directionFileNames{iFlowDirection} '/' opts.scanName '_' directionFileNames{iFlowDirection} '_Vel.nii']);
    [vel,xyz]=spm_read_vols(V_vel,0);
    load([opts.SubjectDir 'AcquisitionInfo_Vel' directionFileNames{iFlowDirection}]);
    venc=info.FlowVenc*10; % does this convert from cm to mm
    vmax=4096;
    
    %% Initialise variables
    NMasks=size(opts.maskNames,2);
    NFrames=size(vel,4);
    timeResSec=(60/(NFrames*info.HR)); %time resolution
    timeResNomSec=(60/(NFrames*60)); %time resolution (assuming HR=60)
    timeSec=((0:(NFrames-1))*timeResSec).'; %time of each frame (first frame=0)
    timeNomSec=((0:(NFrames-1))*timeResNomSec).'; %time of each frame (assuming HR=60)
    flow_mlps=nan(NFrames,NMasks); %flow (m/s)
    flowCorr_mlps=nan(NFrames,NMasks); %flow (corrected for background)
    BGVelocity_mmps=nan(NFrames,NMasks); %background velocity (mm/s)
    area_mm2=nan(1,NMasks); %area of ROI
        
    maskExists=zeros(1,NMasks);
    pixelVelRaw=cell(NFrames,NMasks); %velocity of individual pixels within the ROI (mm/s)
    pixelVel=cell(NFrames,NMasks); %velocity of individual pixels within the ROI after phase unwrapping (mm/s)
    pixelVelBG=cell(NFrames,NMasks); %velocity of individual pixels within the background ROI (mm/s)
    pixelVelCorr=cell(NFrames,NMasks); %velocity of individual pixels within the ROI corrected for background velocity (mm/s)
    
    %% Make figures showing image in background. Masks are added to images later.
    V_disp=spm_vol([opts.FlowImagesDir '/' opts.displayImage '.nii']);
    [displayImage,xyz]=spm_read_vols(V_disp);
    displayImage2=(rot90(displayImage) - min(displayImage(:))) / (max(displayImage(:)) - min(displayImage(:)));
    displayImageRGB=zeros([size(displayImage2) 3]); displayImageRGB(:,:,1)=displayImage2; displayImageRGB(:,:,2)=displayImage2; displayImageRGB(:,:,3)=displayImage2;
    figure(2)%, set(gcf,'Units','Centimeters','OuterPosition',[0 20 25 20],'PaperOrientation','Portrait','PaperType','A4','PaperPositionMode','Auto');
    image(displayImageRGB); axis image; hold on;
    
    %% Loop through masks to get velocities
    for iMask=1:NMasks
        
        if ~exist([opts.masksDir '/' opts.maskNames{iMask} '.nii'],'file'); continue; end %if mask doesn't exist, skip and continue
        maskExists(iMask)=1;
        [mask,temp]=spm_read_vols(spm_vol([opts.masksDir '/' opts.maskNames{iMask} '.nii']));
        if ~isempty(find(mask ~= 1 & mask ~= 0)); error('Mask contains voxels which are not 0/1!'); end; %check masks contain only 0 and 1
        mask(mask(:)==0)=nan; %replace zeros with nans
        area_mm2(iMask)=nansum(mask(:)) * info.PixelSpacing(1) * info.PixelSpacing(2);
        
        [BGMask,temp]=spm_read_vols(spm_vol([opts.masksDir '/' opts.BGMaskName{iMask} '.nii']));
        BGMask(BGMask(:)==0)=nan; %replace zeros with nans
        displayBGMask=rot90(BGMask); if ~isempty(find((BGMask ~= 1) & (~isnan(BGMask)))); error('Background mask contains voxels which are not 0/1!'); end;
        
        for iFrame=1:NFrames %calculate mean signal intensity
            temp2=((squeeze(vel(:,:,:,iFrame)).*BGMask)/vmax)*venc; %Final
            
            pixelVelBG{iFrame,iMask}=(eval(['opts.flowSignCorrection_' directionNames{iFlowDirection} '(' num2str(iMask) ')']) * temp2(~isnan(temp2)));
            BGVelocity_mmps(iFrame,iMask) = ( sum(pixelVelBG{iFrame,iMask}) / nansum(BGMask(:)));
            
            temp3=((squeeze(vel(:,:,:,iFrame)).*mask)/vmax)*venc; %Final
            
            pixelVelRaw{iFrame,iMask}=(eval(['opts.flowSignCorrection_' directionNames{iFlowDirection} '(' num2str(iMask) ')']) * temp3(~isnan(temp3)));
            pixelVel{iFrame,iMask}=(eval(['opts.flowSignCorrection_' directionNames{iFlowDirection} '(' num2str(iMask) ')']) * temp3(~isnan(temp3)));
            
            %% phase unwrapping
            if size(opts.aliasCorrection,1)>1
                centreVelocity = nanmedian(pixelVelRaw{iFrame,iMask}) + opts.aliasCorrection(iFrame,iMask); %allow correction to be different for each frame
            else
                centreVelocity = nanmedian(pixelVelRaw{iFrame,iMask}) + opts.aliasCorrection(iMask); %same correction for all frames
            end
            for iVoxel=1:size(pixelVelRaw{iFrame,iMask},1)
                if pixelVelRaw{iFrame,iMask}(iVoxel) < centreVelocity  - 1*venc; pixelVel{iFrame,iMask}(iVoxel)=pixelVelRaw{iFrame,iMask}(iVoxel) + 2*venc; end
                if pixelVelRaw{iFrame,iMask}(iVoxel) > centreVelocity  + 1*venc; pixelVel{iFrame,iMask}(iVoxel)=pixelVelRaw{iFrame,iMask}(iVoxel) - 2*venc; end
            end
            
            %% correct for background velocity
            pixelVelCorr{iFrame,iMask}=pixelVel{iFrame,iMask}-BGVelocity_mmps(iFrame,iMask);
            
            %% calculate flow
            flow_mlps(iFrame,iMask) = sum(pixelVel{iFrame,iMask}) * info.PixelSpacing(1) * info.PixelSpacing(2) * 0.001; %calculate total flow in ml/s
            flowCorr_mlps(iFrame,iMask) = sum(pixelVelCorr{iFrame,iMask}) * info.PixelSpacing(1) * info.PixelSpacing(2) * 0.001; %calculate total flow (corrected for background) in ml/s
            
            
            
        end;
        
        %% Add mask to image
        displayMask=rot90(mask);
        overlayRGB=zeros(size(displayImageRGB));
        for iRow=1:size(displayMask,1); for iCol=1:size(displayMask,2)
                if displayMask(iRow,iCol)==1 || displayBGMask(iRow,iCol)==1; overlayRGB(iRow,iCol,:)=[1 0 0]; end
            end; end
        figure(2),image(overlayRGB,'alphadata',0.3*overlayRGB(:,:,1));
        x=max(find(sum(displayMask==1))); y=max(find(sum(displayMask.'==1)));
        text(x,y,opts.maskNames{iMask},'Color','y','FontSize',10);
        x=max(find(sum(displayBGMask==1))); y=max(find(sum(displayBGMask.'==1)));
        text(x,y,opts.BGMaskName{iMask},'Color','y','FontSize',10);
        
        
    end
    set(gcf,'Renderer','OpenGL')
    saveas(gcf,[opts.flowDataDir '/flowROIs_' opts.HVNumberStr '_' opts.scanName '_' directionFileNames{iFlowDirection} ],'jpg'); %save image
    
    %% Save data
    save([opts.flowDataDir '/flowData_' directionNames{iFlowDirection}],'NMasks','NFrames','maskExists','timeResSec','timeResNomSec','timeSec','timeNomSec','flow_mlps','flowCorr_mlps','BGVelocity_mmps','pixelVel','pixelVelCorr','venc','area_mm2');
    
    %% Spreadsheet output
    NCols=3+NMasks*4; NRows=1+NFrames; xls=cell(NRows,NCols);
    xls(1,1:3)={'Frame' 'time (s)' 'nominal time (s)'};
    xls(2:NFrames+1,1)=num2cell(1:NFrames);
    xls(2:NFrames+1,2)=num2cell(timeSec);
    xls(2:NFrames+1,3)=num2cell(timeNomSec);
    for iMask=1:NMasks
        xls(1,(iMask*4):(iMask*4+3))={[opts.maskNames{iMask} ' area (mm2)'] [opts.maskNames{iMask} ' flow (ml/s)'] [opts.maskNames{iMask} ' corrected flow (ml/s)'] [opts.maskNames{iMask} ' BG velocity (mm2/s)']};
        xls(2:(NFrames+1),(iMask*4):(iMask*4+3))=num2cell( [ area_mm2(1,iMask)*ones(NFrames,1) flow_mlps(:,iMask) flowCorr_mlps(:,iMask) BGVelocity_mmps(:,iMask) ] );
    end
    if isfield(info,'NomInterval')
        xls2={ 'TR' info.TR; 'TE' info.TE; 'FA' info.FA; 'TrigTime' info.TrigTime; 'NomInterval' info.NomInterval; 'PixelSpacing' num2str(info.PixelSpacing.'); 'Orientation' num2str(info.orientation.'); 'Velocity Encoding' num2str(info.FlowVenc)};
    else
        info.NomInterval=0;
        xls2={ 'TR' info.TR; 'TE' info.TE; 'FA' info.FA; 'TrigTime' info.TrigTime; 'NomInterval' info.NomInterval; 'PixelSpacing' num2str(info.PixelSpacing.'); 'Orientation' num2str(info.orientation.'); 'Velocity Encoding' num2str(info.FlowVenc)};
    end

    
    fid=fopen([opts.flowDataDir '/flowData_' directionFileNames{iFlowDirection} '.csv'],'w');
    
    for m=1:size(xls,2)
        if m==size(xls,2)
            fprintf(fid,'%s \n',xls{1,m});
        else
            fprintf(fid,'%s,',xls{1,m});
        end
    end
    
    for n=1:size(xls,1)-1
        for m=1:size(xls,2)
            if m==1
                fprintf(fid,'%d,',xls{n+1,m});
            elseif m==size(xls,2)
                fprintf(fid,'%f\n',xls{n+1,m});
            else
                fprintf(fid,'%f,',xls{n+1,m});
            end
        end
    end
    fclose(fid);
    
    fid2=fopen([opts.flowDataDir '/flowData_parameters_' directionFileNames{iFlowDirection} '.csv'],'w');
    for n=1:size(xls2,1)
        if n==1 || n==2
            fprintf(fid2,'%s,%f \n',xls2{n,:});
        elseif n==3 || n==4 || n==5 || n==6
            fprintf(fid2,'%s,%d \n',xls2{n,:});
        else
            fprintf(fid2,'%s,%s \n',xls2{n,:});
        end
    end
    fclose(fid2);
    
    %% Plot flow data
    figure(1); set(gcf,'Units','Centimeters','OuterPosition',[10 0 15 25],'PaperOrientation','Portrait','PaperType','A4','PaperPositionMode','Auto');
    figure(4); set(gcf,'Units','Centimeters','OuterPosition',[10 0 15 25],'PaperOrientation','Portrait','PaperType','A4','PaperPositionMode','Auto');
    
    for iMask=1:NMasks
        if maskExists(iMask)==0; continue; end
        
        figure(1)% flow
        
        subplot(NMasks,2,(iMask*2)-1);
        %     x=timeNomSec;
        x=1:NFrames; y=flow_mlps(:,iMask); y2=flowCorr_mlps(:,iMask);
        plot(x,y,'b:',x,y2,'b-');
        title([opts.maskNames{iMask} ' HR=' num2str(info.HR)]);
        line([0 1],[0 0],'Color','k');
        xlim([0 32]); ylim([min([0 y.' y2.']) max([0 y.' y2.'])]);
        %     xlabel('time (nom) (/s)');
        xlabel('timeframe'); ylabel('flow (ml/s)');
        
        subplot(NMasks,2,iMask*2)
        %     x=timeNomSec;
        x=1:NFrames; y=BGVelocity_mmps(:,iMask);
        plot(x,y,'b-');
        title([opts.BGMaskName{iMask} ' HR=' num2str(info.HR)]);
        line([0 1],[0 0],'Color','k');
        xlim([0 32]); ylim([min([0 y.']) max([0 y.'])]);
        %     xlabel('time (nom) (/s)');
        xlabel('timeframe'); ylabel('mean BG velocity (mm/s)');
        
        figure(4); %pixel velocities
        
        subplot(NMasks,2,(iMask*2)-1);
        y=cell2mat(pixelVel(:,iMask).').';
        plot(x,y,'k.')
        %     xlabel('time (nom) (/s)');
        xlabel('timeframe'); ylabel('velocity (mm/s)');
        ylim([min([-venc min(y)]) max([venc max(y)])]);
        title([opts.maskNames{iMask}  ', venc=' num2str(venc) ' mm/s']);
        line([min(x) max(x)],venc*[1 1],'Color','r');
        line([min(x) max(x)],-venc*[1 1],'Color','r');
        
        subplot(NMasks,2,iMask*2)
        plot(x,cell2mat(pixelVelBG(:,iMask).').','k.'); title([opts.BGMaskName{iMask}]);
        %     xlabel('time (nom) (/s)');
        xlabel('timeframe'); ylabel('velocity (mm/s)'); ylim([-venc venc]);
    end
    
    print(1,'-djpeg','-r400',[opts.flowDataDir '/flowCurves_' opts.HVNumberStr '_'  opts.scanName '_' directionFileNames{iFlowDirection}]);
    print(4,'-djpeg','-r400',[opts.flowDataDir '/pixelVelocities_' opts.HVNumberStr '_' opts.scanName '_' directionFileNames{iFlowDirection}]);
    
end
end