function pipeline_finalFlowCalculationsDotProduct(iSubject,iSes)


SubjectSpreadsheet =% location of spreadsheet containing subject IDs, scan info etc;
Subjectdata = readtable(SubjectSpreadsheet);


vesselNames = {'RMCA' 'LMCA' 'RACA' 'LACA' 'RPCA' 'LPCA' 'SSS' 'StS' 'RTS' 'LTS' 'RICA' 'LICA' 'BA'};


Ntimeframes = 25; % number of 4D flow timeframes
[r,c] = size(vesselNames);
variableTable = zeros(c,4); % a matrix for iSubject, with a row for each vessel and a column for each flow variable
flowTable = zeros(c,Ntimeframes); % a table for flow at each timeframe for each vessel
PI = zeros(c,1); RI = zeros(c,1);peakVel = nan(c,1); peakVel2 = nan(c,1); FVP = zeros(c,1);
meanVel = nan(Ntimeframes,1);

ss_data= {'Subject' 'Heart rate (bpm)' 'Temporal res (s)'...
    'RMCA total flow (ml)' 'RMCA mean flow (ml/s)' 'RMCA PI' 'RMCA RI' 'RMCA FVP' 'RMCA max pixel velocity (cm/s)'...
    'LMCA total flow (ml)' 'LMCA mean flow (ml/s)' 'LMCA PI' 'LMCA RI' 'LMCA FVP' 'LMCA max pixel velocity (cm/s)' ...
    'RACA total flow (ml)' 'RACA mean flow (ml/s)' 'RACA PI' 'RACA RI' 'RACA FVP' 'RACA max pixel velocity (cm/s)' ...
    'LACA total flow (ml)' 'LACA mean flow (ml/s)' 'LACA PI' 'LACA RI' 'LACA FVP' 'LACA max pixel velocity (cm/s)' ...
    'RPCA total flow (ml)' 'RPCA mean flow (ml/s)' 'RPCA PI' 'RPCA RI' 'RPCA FVP' 'RPCA max pixel velocity (cm/s)' ...
    'LPCA total flow (ml)' 'LPCA mean flow (ml/s)' 'LPCA PI' 'LPCA RI' 'LPCA FVP' 'LPCA max pixel velocity (cm/s)' ...
    'SSS total flow (ml)' 'SSS mean flow (ml/s)' 'SSS PI' 'SSS RI' 'SSS FVP' 'SSS max pixel velocity (cm/s)' ...
    'StS total flow (ml)' 'StS mean flow (ml/s)' 'StS PI' 'StS RI' 'StS FVP' 'StS max pixel velocity (cm/s)' ...
    'RTS total flow (ml)' 'RTS mean flow (ml/s)' 'RTS PI' 'RTS RI' 'RTS FVP' 'RTS max pixel velocity (cm/s)' ...
    'LTS total flow (ml)' 'LTS mean flow (ml/s)' 'LTS PI' 'LTS RI' 'LTS FVP' 'LTS max pixel velocity (cm/s)' ...
    'RICA total flow (ml)' 'RICA mean flow (ml/s)' 'RICA PI' 'RICA RI' 'RICA FVP' 'RICA max pixel velocity (cm/s)' ...
    'LICA total flow (ml)' 'LICA mean flow (ml/s)' 'LICA PI' 'LICA RI' 'LICA FVP' 'LICA max pixel velocity (cm/s)' ...
    'BA total flow (ml)' 'BA mean flow (ml/s)' 'BA PI' 'BA RI' 'BA FVP' 'BA max pixel velocity (cm/s)'
    };

[~, Nvariables] = size(ss_data);


resultsdir1  = % directory location for this subject's results
resultsdir2 = ([resultsdir1 '/v' num2str(iSes)]); % directory for this visit/session
system(['mkdir ' resultsdir1]); system(['mkdir ' resultsdir2]);
results_spreadsheet = '/flow_results.csv';
if exist(results_spreadsheet); delete results_spreadsheet; end

fvp_func = @(q) max(cumtrapz(q-mean(q)))-min(cumtrapz(q-mean(q))); % this is the single-line function for calculating flow volume pulsatility - is used for each vessel
maxFlowPlot = 0; % this is just so I can keep a record of which vessel has the highest peak flow and scale the plots accordingly

for iFlowScan = 1:13
    %% load data
    flowDataDir  = % directory location of the previously acquired flow data for this patient, visit and vessel
    if ~exist([flowDataDir '/flowData.mat']); continue; end % if the data hasn't been collected (e.g. due to no mask) then move onto next vessel
    
    disp(['Processing ' char(Subjectdata{iSubject,1}) ', visit ' num2str(iSes) ', ' vesselNames{iFlowScan}]);
    
    FD1=load([flowDataDir '/flowData.mat']);
    timeRes = FD1.timeResSec;
    RR_interval = Ntimeframes * FD1.timeResSec;  % time resolution in carotid scans
    HR = 60/RR_interval;
    
    %% flow
    % Flows in mL/s
    FD1.flowCorr_mlps(isnan(FD1.flowCorr_mlps))=0;  % guard against missing values (NaN)
    q = FD1.flowCorr_mlps(:,1);
    
    
    %the next bit of code circshifts all the flow values so that
    %the peak comes at the start - this is for the individual flow graphs
    [min_flow, index_min] = min(q);
    [max_flow, index_max] = max(q);
    q1 = circshift(q,-(index_min-1));
    q_mean = mean(q1); variableTable(iFlowScan,2) = q_mean;
    q_sum = q_mean*RR_interval;  variableTable(iFlowScan,1) = q_sum;     %total flow in a cardiac cycle
    if max_flow > maxFlowPlot
        maxFlowPlot = max_flow;
    end
    
    %this bit of code circshifts all the flow values so that
    % they align with each other - this is for the overall flow graph
    if iFlowScan == 1 % i.e. RMCA
        [~, index_RMCAmin] = min(q);
    end
    
    q2 = circshift(q,-(index_RMCAmin-1)); % flow values for given vessel now shifted to align with other vessels
    flowTable(iFlowScan,:) = q2;
    
    
    %% Calculate pulsatility and resistance index of waveforms
    
    PI(iFlowScan) = (max(q)-min(q))/mean(q); variableTable(iFlowScan,3) = PI(iFlowScan);
    RI(iFlowScan) = (max(q)-min(q))/max(q); variableTable(iFlowScan,4) = RI(iFlowScan);
    FVP(iFlowScan) = fvp_func(q); variableTable(iFlowScan,5) = FVP(iFlowScan);
    
    
    variableTable(iFlowScan,6) = peakVel2(iFlowScan); % max pixel velocity is blank until the code for calculating it is added
    
    %% plot figures
    
    figure(1)
    subplot(2,2,1);% plots the flow values for this vessel
    sgtitle(['Subject' Subjectdata{iSubject,1}]);
    plot(q1,'r','LineWidth',2);
    axis([1 25 -inf inf]);
    xlabel('Timeframe'); ylabel('Flow (ml/s)');
    lgd1 = legend(vesselNames{iFlowScan});
    set(lgd1,'Location','best') % places legend so that it doesn't cover plot
    
    % this bit may be obsolete
    subplot (2,2,2) % bar chart of this vessel's pulsatility
    x = PI(iFlowScan);
    y = vesselNames{iFlowScan};%, "LMCA", "SSS", "StS"];
    bar(x,'facecolor','red')
    set(gca, 'xticklabel',y);
    xlabel('Vessel'); ylabel('Pulsatility Index');
    
    figure1 = ['print -dtiff ' resultsdir2 '/flow_waveforms_' char(Subjectdata{iSubject,1}) '_' vesselNames{iFlowScan} '_v' num2str(iSes) '.tif'];
    eval(figure1);
    saveas(1,[resultsdir2 '/flow_waveforms_' char(Subjectdata{iSubject,1}) '_' vesselNames{iFlowScan} '_v' num2str(iSes)]);
end
%% add variables to matrix

variableTable(variableTable(iFlowScan,:)==0) = NaN; % replace zeros with NaNs

ss_data(iSubject+1,1:3)=num2cell([iSubject HR timeRes]);
for iFlowScan = 1:c
    ss_data(iSubject+1,iFlowScan*6-2:iFlowScan*6+3)=num2cell(variableTable(iFlowScan,1:6));
end

save('alldata.mat','ss_data')

%% save to spreadsheet


fid=fopen([resultsdir2 '/' results_spreadsheet],'w');

for n=1:size(ss_data,2)
    if n==size(ss_data,2)
        fprintf(fid,'%s \n',ss_data{1,n});
    else
        fprintf(fid,'%s,',ss_data{1,n});
    end
end

for n=2:size(ss_data,1)
    for m=1:size(ss_data,2)
        if m==1
            fprintf(fid,'%s,',char(Subjectdata{ss_data{n,m},1})); % first column in string form
        elseif m==Nvariables
            fprintf(fid,'%f\n',ss_data{n,m}); % if final column, new line
        else
            fprintf(fid,'%f,',ss_data{n,m}); % everything else is a floating number
        end
    end
end



% fillmissing(results_spreadsheet,'constant','NaN');

fclose(fid);

figure (2) % flow curves all vessels
subplot(2,1,1)
plot(flowTable(1,:),'LineWidth',2); hold on; % RMCA
plot(flowTable(2,:),'LineWidth',2); hold on; % LMCA
plot(flowTable(3,:),'LineWidth',2); hold on; % RACA
plot(flowTable(4,:),'LineWidth',2); hold on; % LACA
plot(flowTable(5,:),'LineWidth',2); hold on; % RPCA
plot(flowTable(6,:),'LineWidth',2); hold on; % LPCA
plot(flowTable(11,:),'LineWidth',2); hold on; % RICA
plot(flowTable(12,:),'k','LineWidth',2); hold on; % LICA
plot(flowTable(13,:),'r-','LineWidth',2); hold on; % BA
axis([1 25 0 maxFlowPlot]);
xlabel('Timeframe'); ylabel('Flow (ml/s)');
legend(vesselNames{[1:6 11:13]}, 'Location', 'bestoutside');
title('Arteries')

subplot(2,1,2)
plot(flowTable(7,:),'LineWidth',2); hold on;
plot(flowTable(8,:),'LineWidth',2); hold on;
plot(flowTable(9,:),'LineWidth',2); hold on;
plot(flowTable(10,:),'LineWidth',2); hold on;
axis([1 25 -inf maxFlowPlot]);
xlabel('Timeframe'); ylabel('Flow (ml/s)');
legend(vesselNames{7:10}, 'Location', 'bestoutside');
title('Veins')

sgtitle(['Subject' Subjectdata{iSubject,1}]);


figure2 = ['print -dtiff ' resultsdir2 '/flow_' char(Subjectdata{iSubject,1})  '_v' num2str(iSes) '.tif'];
eval(figure2);
saveas(2,[resultsdir2 '/flow_' char(Subjectdata{iSubject,1}) '_v' num2str(iSes)]);

figure(3) % PI and RI all vessels
subplot (2,1,1) % bar chart of this vessel's pulsatility
x = PI;
y = categorical(vesselNames); y = reordercats(y,{'RICA', 'LICA', 'BA', 'RMCA', 'LMCA', 'RACA', 'LACA', 'RPCA', 'LPCA', 'StS', 'SSS', 'RTS', 'LTS'});
N = numel(x);
for i = 1:N
    h = bar(y(i),x(i));
    if i == 1, hold on, end
    if i > 6 && i <=10, col = 'b';
    else col = 'r';
    end
    set (h, 'FaceColor', col)
end
xtickangle(45)
xlabel('Vessel'); ylabel('Pulsatility Index');

subplot (2,1,2)
x = RI;
y = categorical(vesselNames); y = reordercats(y,{'RICA', 'LICA', 'BA', 'RMCA', 'LMCA', 'RACA', 'LACA', 'RPCA', 'LPCA', 'StS', 'SSS', 'RTS', 'LTS'});
N = numel(x);
for i = 1:N
    h = bar(y(i),x(i));
    if i == 1, hold on, end
    if i > 6 && i <=10, col = 'b';
    else col = 'r';
    end
    set (h, 'FaceColor', col)
end

xtickangle(45)
xlabel('Vessel'); ylabel('Resistance Index');


figure3 = ['print -dtiff ' resultsdir2 '/puls_res_' char(Subjectdata{iSubject,1})  '_v' num2str(iSes) '.tif'];
eval(figure3);
saveas(3,[resultsdir2 '/puls_res_' char(Subjectdata{iSubject,1})  '_v' num2str(iSes)]);

close all
end
