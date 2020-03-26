% a function that loads the flow waveforms for a vessel across the 3 main
% methods and plots them
function methodComparisonPlot_function(opts)%,iSubject,iFlowScan)

NFrames = 25;
%% load in data
% basic flow calculated using root of sum squares
basicData = load([opts.flowDataDir1 '/flowData']);
% dot product flow
dpData = load([opts.flowDataDir2 '/flowData']);
% PVC flow
pvcData = load([opts.flowDataDir3 '/flowData']);

% flow values
flowBasic = basicData.flowCorr_mlps;
flowDotProduct = dpData.flowCorr_mlps;
flowPVC = pvcData.flowPVC_mlps;

[~, index_min] = min(flowBasic);
flowBasic = circshift(flowBasic,-(index_min-1));  % flow values shifted so min is last
flowDotProduct = circshift(flowDotProduct,-(index_min-1));
flowPVC = circshift(flowPVC,-(index_min-1));

% ROI areas
areaHand = basicData.area_mm2;
areaRGtotal = pvcData.area_total_mm2;
areaRGlumen = pvcData.area_lumen_mm2;
areaRGwall = areaRGtotal - areaRGlumen;


%% Plot flow data
figure(1); set(gcf,'Units','Centimeters','OuterPosition',[10 0 15 25],'PaperOrientation','Portrait','PaperType','A4','PaperPositionMode','Auto');

x=1:NFrames;  y1=flowBasic; y2=flowDotProduct; y3=flowPVC;
plot(x',y1,'b-',x,y2,'g-',x,y3,'r-'); % raw is dashed line, corrected is solid line, partial volume and BG corrected is red line
title(['Flow extraction comparison: ' opts.maskNames]);
line([0 1],[0 0],'Color','k');
xlim([1 NFrames]); ylim([-inf inf]); %ylim([min([0 y.' y2.']) max([0 y.' y2.'])]);
%     xlabel('time (nom) (/s)');
xlabel('Timeframe'); ylabel('Flow (ml/s)');
legend('Basic', 'Dot product', 'Partial volume corrected');

print(1,'-djpeg','-r400',[opts.SubjectDir2 '/flowMethodsComparison_' opts.HVNumberStr '_'  opts.maskNames]);

%% Mask size bar chart
figure(2); set(gcf,'Units','Centimeters','OuterPosition',[10 0 15 25],'PaperOrientation','Portrait','PaperType','A4','PaperPositionMode','Auto');

%  y = [areaHand; areaRGtotal, areaRGlumen];
y = [areaHand 0 0; 0 areaRGlumen areaRGwall];
bar(y,'stacked');
title (['ROI area comparison: ' opts.maskNames]);
% ylim([-inf inf]);
ylabel('Area (mm2)');
legend('Hand drawn', 'Region-growing (lumen)', 'Region-growing (wall)', 'Location', 'Northwest');
set(gca,'XTickLabel',[]);

print(2,'-djpeg','-r400',[opts.SubjectDir2 '/ROIDrawingComparison_' opts.HVNumberStr '_'  opts.maskNames]);


end

