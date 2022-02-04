% a simple piece of code to run before 4D analysis which saves the subject IDs to 'scanInfo'
% this can be done in a much better way but I just never got round to making a more sophisticated method
clear; 


%% determine number and code of subjects

scanInfo.HVNumberStr = { ...

% subject IDs
'HV001' 'HV002' 'HV003' 'HV004' 'HV005' 'HV006'  'HV007' 'HV008' 'HV009' 'HV010' 'HV011'};


N=size(scanInfo.HVNumberStr,1);

%% initialise variables
scanInfo.N=N;
scanInfo.dateStr=cell(N,2);
scanInfo.examNoStr=cell(N,2);


scanInfo.N=size(scanInfo.HVNumberStr,2);

procRoot= %  location of 4D analysis directory 

cd(procRoot)
save scanInfo scanInfo % make sure in correct directory

