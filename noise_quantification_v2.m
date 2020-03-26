% this function is designed to calculate spatial and temporal noise in 4D
% flow scans. This new version calculates temporal variance by looking at
% the standard deviations of all the pixels within the circular ROI and the
% spatial variance by looking at how the ROI varies within each timeframe
% separately, and calculating the mean standard deviation from all
% timeframes. SNR is also calculated as mean signal/spatial stdev

% REMEMBER TO OPEN UP THE MAG SLAB AND PLACE A 7MM CIRCULAR ROI IN THE
% CENTRE, SAVE AS SUBJECT ID AND CONVERT IT TO NIFTI USING CODE BELOW


function noise_quantification_v2(opts,iSubject)

opengl('software');

close all force
set(0,'defaulttextinterpreter','none','DefaultAxesFontSize',10)

%% read in subject spreadsheet
SubjectSpreadsheet ='/DSTORE/BRICIA/amorgan_PhD/4DFlowProject/SubjectDatabase.xlsx';
Subjectdata = readtable(SubjectSpreadsheet);
centreSlice = (Subjectdata{iSubject,7}/2)+1;

disp(['Noise processing for ' char(Subjectdata{iSubject,1})])

directions = {'velFH','velRL','velAP'};

%% create directory for output or delete files in existing directory
subjectNoiseDir = [opts.noiseQuantDir char(Subjectdata{iSubject,1})];
if ~exist(subjectNoiseDir,'dir'); mkdir(subjectNoiseDir); else delete([subjectNoiseDir '/*']); end;

%% Load 4D velocity image and acquisition parameters
mag_image=spm_vol([opts.SubjectDir 'Mag/combined.nii']);
[mag,~]=spm_read_vols(mag_image,0);
velFH_image=spm_vol([opts.SubjectDir 'VelF-H/combined.nii']);
[velFH,~]=spm_read_vols(velFH_image,0);
velRL_image=spm_vol([opts.SubjectDir 'VelR-L/combined.nii']);
[velRL,~]=spm_read_vols(velRL_image,0);
velAP_image=spm_vol([opts.SubjectDir 'VelA-P/combined.nii']);
[velAP,~]=spm_read_vols(velAP_image,0);

load([opts.SubjectDir '/AcquisitionInfo_VelA-P']); % chosen an arbitrary velocity direction to acquire scan info
%venc=ceil(max(vel(:))); %calculate venc (mm/s) as maximum velocity in image (this should be accurate as long as there is noise in the image)
% venc=info.FlowVenc*10; % does this convert from cm to mm
% vmax=4096;

%% Initialise variables
NMasks=size(opts.maskNames,2);
NFrames=size(mag,4);
BGVelocity_mmps=nan(NFrames,NMasks); %background velocity (mm/s)
area_mm2=nan(1,NMasks); %area of ROI

maskExists=zeros(1,NMasks);
pixelVel=cell(NFrames,NMasks); %velocity of individual pixels within the ROI after phase unwrapping (mm/s)

% use this command in case compressed niftis need converting
system(['for i in $(ls ' opts.noiseMasksDir '*.nii.gz);do fslchfiletype NIFTI ${i};done']); % convert compressed niftis to regular niftis
        gunzip([opts.noiseMasksDir '/*.nii.gz']);


%% Make figures showing image in background. Masks are added to images later.
V_disp=spm_vol([opts.SubjectDir '/Mag/vol0000.nii']);
[displayImage,xyz]=spm_read_vols(V_disp);
displayImage=displayImage(:,:,centreSlice); % select middle slice
displayImage2=(rot90(displayImage) - min(displayImage(:))) / (max(displayImage(:)) - min(displayImage(:)));
displayImageRGB=zeros([size(displayImage2) 3]);
displayImageRGB(:,:,1)=displayImage2;
displayImageRGB(:,:,2)=displayImage2;
displayImageRGB(:,:,3)=displayImage2;
figure(1)%, set(gcf,'Units','Centimeters','OuterPosition',[0 20 25 20],'PaperOrientation','Portrait','PaperType','A4','PaperPositionMode','Auto');
image(displayImageRGB); axis image; hold on;
title(['Centre slice: HV0' num2str(iSubject,'%02.f')]);

%% Loop through masks to get velocities


    if ~exist([opts.noiseMasksDir '/HV0' num2str(iSubject,'%02.f') '.nii'],'file'); return; end %if mask doesn't exist, skip and continue
    [mask,~]=spm_read_vols(spm_vol([opts.noiseMasksDir '/HV0' num2str(iSubject,'%02.f') '.nii'])); %load mask
    
    if ~isempty(find(mask ~= 1 & mask ~= 0)); error('Mask contains voxels which are not 0/1!'); end %check masks contain only 0 and 1
    mask(mask(:)==0)=nan; %replace zeros with nans
    mask = mask(:,:,centreSlice);
    
    pixelN = nansum(mask(:));
    area_mm2= pixelN * info.PixelSpacing(1) * info.PixelSpacing(2);
    VNR = nan(3,1);
    
    %% Magnitude image%%

pixels = nan(pixelN,NFrames);
temporalMeans = nan(pixelN,1); temporalSD = nan(pixelN,1);
               
    for iFrame=1:NFrames %calculate mean signal intensity
        
        temp2=squeeze(mag(:,:,centreSlice,iFrame)).*mask; % this finds the pixel intensities within the mask
        
        % find the mean value across all mask pixels at timeframe 1 (chosen arbitrarily)
        if iFrame == 1
            A = temp2(~isnan(temp2));
            meanPixelInt_space_mag = mean(A); stdevPixelInt_space_mag = std(A);
        end
        
        pixels(:,iFrame) = temp2(~isnan(temp2));

    end
    
    for i = 1:pixelN
        temporalMeans(i) = mean(pixels(i,:)); temporalSD(i) = std(pixels(i,:));
    end
    
%     spatialMeans = mean(pixels); spatialSD = std(pixels);
%     meanPixelInt_space_mag = mean(spatialMeans); stdevPixelInt_space_mag = mean(spatialSD);
    meanPixelInt_time_mag = mean(temporalMeans); stdevPixelInt_time_mag = mean(temporalSD);
    SNR = meanPixelInt_space_mag/stdevPixelInt_space_mag;
    
        %% Velocity images%%
        
 meanPixelInt_space_vels = nan(3,1); stdevPixelInt_space_vels = nan(3,1);
 meanPixelInt_time_vels = nan(3,1); stdevPixelInt_time_vels = nan(3,1);
        
for direction = 1:3
pixels = nan(pixelN,NFrames);
temporalMeans = nan(pixelN,1); temporalSD = nan(pixelN,1);
               
    for iFrame=1:NFrames %calculate mean signal intensity
        
        temp2=eval(['squeeze(' directions{direction} '(:,:,centreSlice,iFrame)).*mask']); % this finds the pixel intensities within the mask
        
        % find the mean value across all mask pixels at timeframe 1 (chosen arbitrarily)
        if iFrame == 1
            A = temp2(~isnan(temp2));
            meanPixelInt_space_vels(direction) = mean(A); stdevPixelInt_space_vels(direction) = std(A);
        end
        
        pixels(:,iFrame) = temp2(~isnan(temp2));

    end
    
    for i = 1:pixelN
        temporalMeans(i) = mean(pixels(i,:)); temporalSD(i) = std(pixels(i,:));
    end
    
%     spatialMeans = mean(pixels); spatialSD = std(pixels);
%     meanPixelInt_space_vels(direction) = mean(spatialMeans); stdevPixelInt_space_vels(direction) = mean(spatialSD);
    meanPixelInt_time_vels(direction) = mean(temporalMeans); stdevPixelInt_time_vels(direction) = mean(temporalSD);
    VNR(direction) = meanPixelInt_space_vels(direction)/stdevPixelInt_space_vels(direction);
end

    
    %% Add mask to image
    displayMask=rot90(mask);
    overlayRGB=zeros(size(displayImageRGB));
    for iRow=1:size(displayMask,1);
        for iCol=1:size(displayMask,2)
            if displayMask(iRow,iCol)==1; overlayRGB(iRow,iCol,:)=[1 0 0]; end % show ROI as red
        end
    end
    
%     overlayRGB(centreCol,centreRow,:)=[0 0 1]; % show single pixel as blue
    figure(1),image(overlayRGB,'alphadata',0.3*overlayRGB(:,:,1));
    x=max(find(sum(displayMask==1))); y=max(find(sum(displayMask.'==1)));
    text(x,y,'Noise ROI','Color','y','FontSize',10);
%     x=max(find(sum(displayBGMask==1))); y=max(find(sum(displayBGMask.'==1)));
%     text(x,y,opts.BGMaskName{iMask},'Color','y','FontSize',10);
    

set(gcf,'Renderer','OpenGL')
saveas(gcf,[subjectNoiseDir '/noiseROI_' opts.HVNumberStr],'jpg'); %save image

%% Save data
save([subjectNoiseDir '/noiseData'], 'meanPixelInt_time_mag', 'stdevPixelInt_time_mag', 'meanPixelInt_space_mag', 'stdevPixelInt_space_mag', 'meanPixelInt_time_vels', 'stdevPixelInt_time_vels', 'meanPixelInt_space_vels', 'stdevPixelInt_space_vels','VNR','area_mm2');

%% Plot bar charts 
figure(2); %%%% magnitude 
hold on
title(['Magnitude noise quantification: HV0' num2str(iSubject,'%02.f')]);
y = [meanPixelInt_space_mag meanPixelInt_time_mag];
y2 = [stdevPixelInt_space_mag stdevPixelInt_time_mag];
ylabel('Mean pixel signal');

% here plot the mean and stdevs for spatial and temporal noise
 bar(y);
 errorbar(y, y2,'.')
xticks([1 2]);
xticklabels({'Spatial','Temporal'});
gtext(['SNR: ' num2str(round(SNR,2))]);
disp('Click desired location of SNR label')
 

print(2,'-djpeg','-r400',[subjectNoiseDir '/magnitudeNoiseChart_' opts.HVNumberStr]);

figure(3); %%%% velocity
hold on
title(['Velocity noise quantification: HV0' num2str(iSubject,'%02.f')]);
% y = [meanPixelInt_space_vels;meanPixelInt_time_vels];
y2 = [stdevPixelInt_space_vels'; stdevPixelInt_time_vels'];
ylabel('Pixel velocity standard deviation');

% here plot the mean and stdevs for spatial and temporal noise
 bar(y2);
%  errorbar(y, y2,'.')
xticks([1:2]);
xticklabels({'Spatial','Temporal'});
legend(directions{1},directions{2},directions{3})
% gtext(['F-H VNR: ' num2str(round(VNR(1),2)); 'R-L VNR: ' num2str(round(VNR(2),2)); 'A-P VNR: ' num2str(round(VNR(3),2))]);
% disp('Click desired locations of VNR labels')

print(3,'-djpeg','-r400',[subjectNoiseDir '/velocitiesNoiseChart_' opts.HVNumberStr]);

end