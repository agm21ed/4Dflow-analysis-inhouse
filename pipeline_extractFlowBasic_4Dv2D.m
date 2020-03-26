% this function extracts the velocity values from the vessel ROIs, and
% combines the x,y and z direction components into throughplane velocity.
% If the vector normal to the ROI does not contain a component, the
% relevant velocities are not included

% this code is also my attempt at calculating a more accurate pixel area
% for the resliced 2D image - rather than just using the voxel dimensions
% from the scan itself

function pipeline_extractFlowBasic_4Dv2D(opts,iSubject,iFlowScan)

if ~exist([opts.masksDir '/' opts.maskNames '.nii'],'file'); return; end %if mask doesn't exist, skip and continue
opengl('software');

close all force
set(0,'defaulttextinterpreter','none','DefaultAxesFontSize',10)

directionNames = {'FH' 'RL' 'AP'}; % named this way because variables can't have a '-' in
directionFileNames = {'F-H' 'R-L' 'A-P'}; % this exists to call the niftis with these names

% call in normal vector components
VectorSpreadsheet = '/DSTORE/BRICIA/amorgan_PhD/4DFlowProject/4Dv2Danalysis/Vessel_centres+normals.xlsx';
VesselLocations = readtable(VectorSpreadsheet); % want to be able to call in vector components from spreadsheet

% centre = [str2double(VesselLocations{iSubject+1, iFlowScan*7 -4}) 180-(str2double(VesselLocations{iSubject+1, iFlowScan*7 -5}))  str2double(VesselLocations{iSubject+1, iFlowScan*7 -3})];
% dist = [str2double(VesselLocations{iSubject+1, iFlowScan*7 -1}) 180-(str2double(VesselLocations{iSubject+1, iFlowScan*7 -2})) str2double(VesselLocations{iSubject+1, iFlowScan*7})];

centre = [str2double(VesselLocations{iSubject+1, iFlowScan*7 -4}) (str2double(VesselLocations{iSubject+1, iFlowScan*7 -5}))  str2double(VesselLocations{iSubject+1, iFlowScan*7 -3})];
dist = [str2double(VesselLocations{iSubject+1, iFlowScan*7 -1}) (str2double(VesselLocations{iSubject+1, iFlowScan*7 -2})) str2double(VesselLocations{iSubject+1, iFlowScan*7})];
vec = dist - centre; normx = -1*vec(1); normy = vec(2); normz = vec(3);

%% create directory for output or delete files in existing directory
if ~exist(opts.flowDataDir1,'dir'); mkdir(opts.flowDataDir1);
else
    delete([opts.flowDataDir1 '/*']);
end

%% Initialise variables
NFrames=25; %size(vel,4);

directionVelocities = cell(NFrames,4); %a matrix with one column for each direction plus one for combined velocity
directionVelocitiesCorr = cell(NFrames,4);
BGvectorMatrix = cell(3,NFrames);
flow_mlps=nan(NFrames); %flow (m/s)
flowCorr_mlps=nan(NFrames); %flow (corrected for background)

for iFlowDirection = 1:3
    %% Load 4D velocity image and acquisition parameters
    V_vel=spm_vol([opts.FlowImagesDir '/' directionFileNames{iFlowDirection} '/' opts.maskNames '_' directionFileNames{iFlowDirection} '_Vel.nii']);
    [vel,~]=spm_read_vols(V_vel,0);
    load([opts.SubjectDir 'AcquisitionInfo_Vel' directionFileNames{iFlowDirection}]);
    venc=info.FlowVenc*10; % does this convert from cm to mm
    vmax=4095;
    
    %% Initialise variables
    timeResSec=(60/(NFrames*info.HR)); %time resolution
    timeResNomSec=(60/(NFrames*60)); %time resolution (assuming HR=60)
    timeSec=((0:(NFrames-1))*timeResSec).'; %time of each frame (first frame=0)
    timeNomSec=((0:(NFrames-1))*timeResNomSec).'; %time of each frame (assuming HR=60)
    flow_mlps=nan(NFrames,1); %flow (m/s)
    flowCorr_mlps=nan(NFrames,1); %flow (corrected for background)
    BGVelocity_mmps=nan(NFrames,1); %background velocity (mm/s)
    
    area_mm2=nan; %area of ROI
    
    maskExists=0;
    pixelVelRaw=cell(NFrames,1); %velocity of individual pixels within the ROI (mm/s)
    pixelVel=cell(NFrames,1); %velocity of individual pixels within the ROI after phase unwrapping (mm/s)
    pixelVelBG=cell(NFrames,1); %velocity of individual pixels within the background ROI (mm/s)
    pixelVelCorr=cell(NFrames,1); %velocity of individual pixels within the ROI corrected for background velocity (mm/s)
    
    % load in background phase model
%     load(['/DSTORE/BRICIA/amorgan_PhD/4DFlowProject/Background_phase_fitting/HV0' num2str(iSubject,'%02.f')  '/FittedValues_' directionFileNames{iFlowDirection} '.mat'],'fitted_value'); % load in fitted BG model

    
    
    %% Make figures showing image in background. Masks are added to images later.
    if iFlowDirection == 1 % only need one ROI image
        V_disp=spm_vol([opts.FlowImagesDir '/' opts.displayImage '.nii']);
        [displayImage,xyz]=spm_read_vols(V_disp);
%         displayImage2=(rot90(displayImage) - min(displayImage(:))) / (max(displayImage(:)) - min(displayImage(:)));
        displayImage2=((displayImage) - min(displayImage(:))) / (max(displayImage(:)) - min(displayImage(:)));
        displayImageRGB=zeros([size(displayImage2) 3]); displayImageRGB(:,:,1)=displayImage2; displayImageRGB(:,:,2)=displayImage2; displayImageRGB(:,:,3)=displayImage2;
        figure(2)%, set(gcf,'Units','Centimeters','OuterPosition',[0 20 25 20],'PaperOrientation','Portrait','PaperType','A4','PaperPositionMode','Auto');
        image(displayImageRGB); axis image; hold on; title(['Basic flow extraction: ' opts.maskNames]);
    end
    
    
    %% calculate new pixel area to more accurately calculate flow - needs to be done in 2 sections I think
    if iFlowDirection == 1 % only calculate new pixel dimensions on first loop
        
        [requiresTrigH, requiresTrigW, originalDimensions] = normalVectorFunction(normx,normy,normz); % this function tells us which (if any) dimensions need to be calculated
        info.PixelSpacing
        
        % section 1: x and y angle to tell us the pixel width
        
        if requiresTrigW ==1
            theta1 = atand(normy/normx); % calculate angle of vector, convert to degrees
            %         sliceAngle1 = 90 - theta1; % angle of slice is perpedicular to angle of normal
            if theta1 >= 0 && theta1 <= 90 % if the angle is positive, subtract it from 90 to get perpendicular angle
                sliceAngle1 = 90 - theta1;
            elseif theta1 < 0 && theta1 >= -90 % if the angle is negative, add 90 to it to get perpendicular angle
                sliceAngle1 = theta1 + 90;
            else
                print('error! angle too large');
            end
            % maximum pixel height is when sliced at 45 degrees, above this and the
            % height reduces back to the original isotropic value - the following
            % code allows for this mirroring effect
            if sliceAngle1 > 45 && sliceAngle1 <= 90
                sliceAngle1 = 45 - (sliceAngle1-45);
            elseif sliceAngle1 > 90
                error('error! angle too large');
            else % do nothing
            end
            % now use sliceAngle to calculate new dimension of 2D plane pixels
            newDimension1 = info.PixelSpacing(2)/(cosd(sliceAngle1));
        else newDimension1 = info.PixelSpacing(1);
        end
        
        
        % section 2: x/y and z angle to tell us pixel height
        
        if requiresTrigH == 1
            if abs(normx) > 0 && abs(normz) > 0 % if x and z components exist, do this method
                theta2 = atand(normz/normx);
            else % otherwise use y and z to calculate pixel height
                theta2 = atand(normz/normy);
            end
            
            %         sliceAngle2 = 90 - theta2; % angle of slice is perpedicular to angle of normal
            if theta2 >= 0 && theta2 <= 90 % if the angle is positive, subtract it from 90 to get perpendicular angle
                sliceAngle2 = 90 - theta2;
            elseif theta2 < 0 && theta2 >= -90 % if the angle is negative, add 90 to it to get perpendicular angle
                sliceAngle2 = theta2 + 90;
            else
                error('error! angle too large');
            end
            % maximum pixel height is when sliced at 45 degrees, above this and the
            % height reduces back to the original isotropic value - the following
            % code allows for this mirroring effect
            if sliceAngle2 > 45 && sliceAngle2 <= 90
                sliceAngle2 = 45 - (sliceAngle2-45);
            elseif sliceAngle2 > 90
                error('error! angle too large');
            else % do nothing
            end
            
            % now use sliceAngle to calculate new dimension of 2D plane pixels
            % - remember pixelspacing(1) and (2) are width and height
            % respectively
            newDimension2 = info.PixelSpacing(1)/cosd(sliceAngle2);
        else
            newDimension2 = info.PixelSpacing(2); % if either x or y components are 0, keep original size
        end
        
        if originalDimensions == 1
            newDimension1 = info.PixelSpacing(1); newDimension2 = info.PixelSpacing(2);
        end
    end
    
    %% Loop through masks to get velocities
    
    if ~exist([opts.masksDir '/' opts.maskNames '.nii'],'file'); return; end %if mask doesn't exist, skip and continue
    maskExists=1;
    [mask,temp]=spm_read_vols(spm_vol([opts.masksDir '/' opts.maskNames '.nii']));
    if ~isempty(find(mask ~= 1 & mask ~= 0)); error('Mask contains voxels which are not 0/1!'); end; %check masks contain only 0 and 1
    mask(mask(:)==0)=nan; %replace zeros with nans
    area_mm2=nansum(mask(:)) * newDimension1 * newDimension2;
    
    [BGMask,temp]=spm_read_vols(spm_vol([opts.masksDir '/' opts.BGMaskName '.nii']));
    BGMask(BGMask(:)==0)=nan; %replace zeros with nans
%     displayBGMask=rot90(BGMask); if ~isempty(find((BGMask ~= 1) & (~isnan(BGMask)))); error('Background mask contains voxels which are not 0/1!'); end;
    displayBGMask=BGMask;  if ~isempty(find((BGMask ~= 1) & (~isnan(BGMask)))); error('Background mask contains voxels which are not 0/1!'); end;
    
    for iFrame=1:NFrames %calculate mean signal intensity
%         min = min(min(squeeze(vel(:,:,:,iFrame))))
%         max = max(max(squeeze(vel(:,:,:,iFrame))))
        temp2=((squeeze(vel(:,:,:,iFrame)).*BGMask)/vmax)*venc; %Final
        
%         pixelVelBG{iFrame}=(eval(['opts.flowSignCorrection_' directionNames{iFlowDirection}]) * temp2(~isnan(temp2)));
        pixelVelBG{iFrame}=temp2(~isnan(temp2));
        BGVelocity_mmps(iFrame) = ( sum(pixelVelBG{iFrame}) / nansum(BGMask(:)));
        
        temp3=((squeeze(vel(:,:,:,iFrame)).*mask)/vmax)*venc; %Final
        
%         pixelVelRaw{iFrame}=(eval(['opts.flowSignCorrection_' directionNames{iFlowDirection}]) * temp3(~isnan(temp3)));
          pixelVelRaw{iFrame}=temp3(~isnan(temp3));
%         pixelVel{iFrame}=(eval(['opts.flowSignCorrection_' directionNames{iFlowDirection}]) * temp3(~isnan(temp3)));
          pixelVel{iFrame}= temp3(~isnan(temp3));
        
        %% phase unwrapping
        if size(opts.aliasCorrection,1)>1
            centreVelocity = nanmedian(pixelVelRaw{iFrame}) + opts.aliasCorrection(iFrame); %allow correction to be different for each frame
        else
            centreVelocity = nanmedian(pixelVelRaw{iFrame}) + opts.aliasCorrection; %same correction for all frames
        end
        for iVoxel=1:size(pixelVelRaw{iFrame},1)
            if pixelVelRaw{iFrame}(iVoxel) < centreVelocity  - 1*venc; pixelVel{iFrame}(iVoxel)=pixelVelRaw{iFrame}(iVoxel) + 2*venc; end
            if pixelVelRaw{iFrame}(iVoxel) > centreVelocity  + 1*venc; pixelVel{iFrame}(iVoxel)=pixelVelRaw{iFrame}(iVoxel) - 2*venc; end
        end
        
        %% correct for background velocity
        pixelVelCorr{iFrame}=pixelVel{iFrame}-BGVelocity_mmps(iFrame);
        BGvectorMatrix{iFlowDirection,iFrame} = pixelVelBG{iFrame}; % for plotting BG velocities later
        
%         pixelVelCorrBGModel{iFrame}=pixelVel{iFrame}-fitted_value(centre(2),centre(1),centre(3))*10; % x10 turns cm to mm to match rest of code
   
        %% enter velocities into direction matrix
        directionVelocities{iFrame,iFlowDirection} = pixelVel{iFrame}; % uncorrected
        directionVelocitiesCorr{iFrame,iFlowDirection} = pixelVelCorr{iFrame}; % corrected
%         directionVelocitiesCorrBGModel{iFrame,iFlowDirection} = pixelVelCorrBGModel{iFrame}; % corrected    
        
    end
    
    %% correct for background velocity using model
%     load(['/DSTORE/BRICIA/amorgan_PhD/4DFlowProject/Background_phase_fitting/HV001/FittedValues_' directionFileNames{iFlowDirection} '.mat'],'fitted_value'); % load in fitted BG model
%     pixelVelCorr=pixelVel{iFrame}-fitted_value(centre(1),centre(2),centre(3));     % Use vessel coordinates to find correct BG phase in model
%     
%     directionVelocities{:,iFlowDirection} = pixelVel; % uncorrected
%     directionVelocitiesCorr{:,iFlowDirection} = pixelVelCorr; % corrected
    
end
%% calculate combined velocity across x, y and z directions and then calculate flow from this
% pixelNo = size(pixelVel{1,1}); pixelNo = pixelNo(1); % calculate number of pixels, then specify row size only
% if normx == 0 % remove velocity components if the normal vector doesn't include it e.g. [0 1 4]
%     for iFrame = 1:NFrames
%         directionVelocities{iFrame,1} = zeros(pixelNo,1); directionVelocitiesCorr{iFrame,1} = zeros(pixelNo,1);
%     end
% end
% if normy == 0
%     for iFrame = 1:NFrames
%         directionVelocities{iFrame,2} = zeros(pixelNo,1); directionVelocitiesCorr{iFrame,2} = zeros(pixelNo,1);
%     end
% end
% if normz == 0
%     for iFrame = 1:NFrames
%         directionVelocities{iFrame,3} = zeros(pixelNo,1); directionVelocitiesCorr{iFrame,3} = zeros(pixelNo,1);
%     end
% end

for iFrame = 1:NFrames % this needs to go through every ROI pixel at a given timeframe and calculate combined velocity
    directionVelocities{iFrame,4} = sqrt(directionVelocities{iFrame,1}.^2 + directionVelocities{iFrame,2}.^2 + directionVelocities{iFrame,3}.^2); % uncorrected
    directionVelocitiesCorr{iFrame,4} = sqrt(directionVelocitiesCorr{iFrame,1}.^2 + directionVelocitiesCorr{iFrame,2}.^2 + directionVelocitiesCorr{iFrame,3}.^2); % corrected
%     directionVelocitiesCorrBGModel{iFrame,4} = sqrt(directionVelocitiesCorrBGModel{iFrame,1}.^2 + directionVelocitiesCorrBGModel{iFrame,2}.^2 + directionVelocitiesCorrBGModel{iFrame,3}.^2); % corrected
    flow_mlps(iFrame) = sum(directionVelocities{iFrame,4}) * newDimension1 * newDimension2 * 0.001; %calculate total flow in ml/s
    flowCorr_mlps(iFrame) = sum(directionVelocitiesCorr{iFrame,4})* newDimension1 * newDimension2 * 0.001; %calculate total flow (corrected for background) in ml/s
%     flowCorrBGModel_mlps(iFrame) = sum(directionVelocitiesCorrBGModel{iFrame,4})* newDimension1 * newDimension2 * 0.001; %calculate total flow (corrected for background) in ml/s
end



%% Add mask to image
% displayMask=rot90(mask);
displayMask=mask;
overlayRGB=zeros(size(displayImageRGB));
for iRow=1:size(displayMask,1); 
    for iCol=1:size(displayMask,2)
        if displayMask(iRow,iCol)==1 || displayBGMask(iRow,iCol)==1; overlayRGB(iRow,iCol,:)=[1 0 0]; end
    end
end
figure(2),image(overlayRGB,'alphadata',0.3*overlayRGB(:,:,1));
x=max(find(sum(displayMask==1))); y=max(find(sum(displayMask.'==1)));
text(x,y,opts.maskNames,'Color','y','FontSize',10);
x=max(find(sum(displayBGMask==1))); y=max(find(sum(displayBGMask.'==1)));
text(x,y,opts.BGMaskName,'Color','y','FontSize',10);


set(gcf,'Renderer','OpenGL')
saveas(gcf,[opts.flowDataDir1 '/flowROIs_' opts.HVNumberStr '_' opts.maskNames],'jpg'); %save image

%% Save data
save([opts.flowDataDir1 '/flowData'],'NFrames','maskExists','timeResSec','timeResNomSec','timeSec','timeNomSec','flow_mlps','flowCorr_mlps','BGVelocity_mmps','directionVelocities','directionVelocitiesCorr','venc','area_mm2');

%% Spreadsheet output
NCols=6; NRows=1+NFrames; xls=cell(NRows,NCols);
xls(1,1:3)={'Frame' 'time (s)' 'nominal time (s)'}; % column headers
xls(2:NFrames+1,1)=num2cell(1:NFrames); % timeframes
xls(2:NFrames+1,2)=num2cell(timeSec); % time (s)
xls(2:NFrames+1,3)=num2cell(timeNomSec); % time (nominal)
xls(1,4:6)={[opts.maskNames ' area (mm2)'] [opts.maskNames ' flow (ml/s)'] [opts.maskNames ' corrected flow (ml/s)']}; % column headers
xls(2:(NFrames+1),4:6)=num2cell( [ area_mm2*ones(NFrames,1) flow_mlps flowCorr_mlps  ] );

if isfield(info,'NomInterval')
    xls2={ 'TR' info.TR; 'TE' info.TE; 'FA' info.FA; 'TrigTime' info.TrigTime; 'NomInterval' info.NomInterval; 'PixelSpacing' num2str(info.PixelSpacing.'); 'Orientation' num2str(info.orientation.'); 'Velocity Encoding' num2str(info.FlowVenc)};
else
    info.NomInterval=0;
    xls2={ 'TR' info.TR; 'TE' info.TE; 'FA' info.FA; 'TrigTime' info.TrigTime; 'NomInterval' info.NomInterval; 'PixelSpacing' num2str(info.PixelSpacing.'); 'Orientation' num2str(info.orientation.'); 'Velocity Encoding' num2str(info.FlowVenc)};
end


fid=fopen([opts.flowDataDir1 '/flowData_' directionFileNames{iFlowDirection} '.csv'],'w');

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

fid2=fopen([opts.flowDataDir1 '/flowData_parameters.csv'],'w');
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
figure(1);% set(gcf,'Units','Centimeters','OuterPosition',[10 0 15 25],'PaperOrientation','Portrait','PaperType','A4','PaperPositionMode','Auto');

if maskExists==0; return; end

figure(1)% flow

%     x=timeNomSec;
%WITH BG MODEL CORRECTION
% x=1:NFrames; y=flow_mlps(:); y2=flowCorr_mlps(:); y3=flowCorrBGModel_mlps(:);
% plot(x,y,'b--',x,y2,'b-',x,y3,'r-'); % raw is dotted line, corrected is solid line
% title(['Basic flow extraction: ' opts.maskNames ' (HR=' num2str(info.HR) ')']);
% legend({'Uncorrected', 'BG corrected','BG model corrected'}, 'Location', 'bestoutside');

% WITHOUT BG MODEL CORRECTION
x=1:NFrames; y=flow_mlps(:); y2=flowCorr_mlps(:);
plot(x,y,'b--',x,y2,'b-'); % raw is dotted line, corrected is solid line
title(['Basic flow extraction: ' opts.maskNames ' (HR=' num2str(info.HR) ')']);
legend({'Uncorrected', 'BG corrected'}, 'Location', 'bestoutside');

% line([0 1],[0 0],'Color','k');
xlim([1 NFrames]); ylim([-inf inf]); %ylim([min([0 y.' y2.']) max([0 y.' y2.'])]);
%     xlabel('time (nom) (/s)');
xlabel('timeframe'); ylabel('flow (ml/s)');


print(1,'-djpeg','-r400',[opts.flowDataDir1 '/flowCurves_' opts.HVNumberStr '_'  opts.maskNames]);

%% plot pixel velocities 

    % plot pixel velocities in each direction, each set
    % having a different colour
        
    figure(3); %pixel velocities     
    subplot(1,2,1);
    sgtitle(['Basic Flow Extraction: ' opts.maskNames ' (venc=' num2str(venc/10) ' cm/s)']);
    
    NPixels = size(pixelVelRaw{iFrame},1);
    if iSubject == 3
        for iFlowDirection = 1:3
            for i = 1:NPixels
                for iFrame = 1:NFrames
                    y = directionVelocities{iFrame,iFlowDirection}(i);
                    if iFlowDirection == 1; plot(iFrame,y,'bx'); hold on; end
                    if iFlowDirection == 2; plot(iFrame,y,'gd'); hold on; end
                    if iFlowDirection == 3; plot(iFrame,y,'ro'); hold on; end
                end
            end
        end
    else
        for iFlowDirection = 1:3
            for i = 1:NPixels
                for iFrame = 1:NFrames
                    y = directionVelocities{iFrame,iFlowDirection}(i);
                    if iFlowDirection == 1; plot(iFrame,y,'ro'); hold on; end
                    if iFlowDirection == 2; plot(iFrame,y,'bx'); hold on; end
                    if iFlowDirection == 3; plot(iFrame,y,'gd'); hold on; end
                end
            end
        end
    end
        
        % need this code to get correct legend
        h = zeros(3, 1);
        h(1) = plot(NaN,NaN,'ro');
        h(2) = plot(NaN,NaN,'bx');
        h(3) = plot(NaN,NaN,'gd');
        legend(h,directionFileNames);
    
    xlabel('timeframe'); ylabel('velocity (mm/s)');
    xlim([0 25]); ylim([min([-venc min(y)]) max([venc max(y)])]);
    title('Pixel velocities'); 
%     line([min(x) max(x)],venc*[1 1],'Color','r'); % lines depecting venc
%     line([min(x) max(x)],-venc*[1 1],'Color','r');



               
       
    subplot(1,2,2) % background velocity
    
    NPixels = size(pixelVelBG{iFrame},1);
    for iFlowDirection = 1:3
        for i = 1:NPixels
            for iFrame = 1:NFrames
                y = BGvectorMatrix{iFlowDirection,iFrame}(i);
                if iFlowDirection == 1; plot(iFrame,y,'ro'); hold on; end
                if iFlowDirection == 2; plot(iFrame,y,'bx'); hold on; end
                if iFlowDirection == 3; plot(iFrame,y,'gd'); hold on; end
            end
        end
    end
    
    h(1) = plot(NaN,NaN,'ro');
    h(2) = plot(NaN,NaN,'bx');
    h(3) = plot(NaN,NaN,'gd');
    legend(h,directionFileNames);
    title('Background velocities');
    xlabel('timeframe'); ylabel('velocity (mm/s)'); ylim([-venc venc]); xlim([0 25]);

    
    print(3,'-djpeg','-r400',[opts.flowDataDir1 '/pixelVelocities_' opts.HVNumberStr '_'  opts.maskNames]);

end
