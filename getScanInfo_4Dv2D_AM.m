clear; 

[DB_MR_NUM,DB_MR_TXT,DB_MR]=xlsread('/DSTORE/BRICIA/amorgan_PhD/4DFlowProject/SubjectDatabase.xlsx');

%% determine number and code of subjects

scanInfo.HVNumberStr = { ...


% flow variable
'HV001' 'HV002' 'HV003' 'HV004' 'HV005' 'HV006'  'HV007' 'HV008' 'HV009' 'HV010' 'HV011'};


N=size(scanInfo.HVNumberStr,1);

%% initialise variables
scanInfo.N=N;
scanInfo.dateStr=cell(N,2);
scanInfo.examNoStr=cell(N,2);

scanNames={'TOF' 'TOF';'ACA' 'ACA';'RMCA' 'RMCA';'Ve' 'Veins'}; %each row gives field name for scan and corresponding name used in database
for iScan=1:size(scanNames,1) %define fields for series numbers and dicom dirs of all scans
    scanInfo.([scanNames{iScan,1} 'SeriesNo'])=nan(N,2);
    scanInfo.([scanNames{iScan,1} 'DicomDir'])=cell(N,2);
end

scanInfo.N=size(scanInfo.HVNumberStr,2);

procRoot='/DSTORE/BRICIA/amorgan_PhD/4DFlowProject/4Dv2Danalysis';

%% extract data
for iSubject=[1:5,7:scanInfo.N]
    
    
    iSubject
    row = find(strcmp(scanInfo.HVNumberStr(iSubject) , DB_MR(2:end,strcmp(DB_MR(1,:),'Subject')))) + 1 %find correct row in MR databse
    row=row(1);
    
    %for iSes=1:2
        iSes=1
    
    yyyy = DB_MR{row,strcmp(DB_MR(1,:),['Date YYYY'])}; 
    mm = DB_MR{row,strcmp(DB_MR(1,:),['Date MM'])};
    dd = DB_MR{row,strcmp(DB_MR(1,:),['Date DD'])};
    scanInfo.dateStr{iSubject,iSes}=[num2str(yyyy,'%04g') num2str(mm,'%02g') num2str(dd,'%02g')];

% dicom location
    for iScan=1:size(scanNames,1)
        iScan
            scanInfo.([scanNames{iScan,1} 'SeriesNo'])(iSubject,iSes)=DB_MR{row,strcmp(DB_MR(1,:),scanNames{iScan,2})};
            temp_name=char(DB_MR_TXT(row,1));
            direc=dir(['/DSTORE/BRICIA/E182012_Brain_4D_Flow_RO/' temp_name(1:end) '/' scanInfo.dateStr{iSubject} '*']);
            size(scanNames)
            [scanNames{iScan,1} 'DicomDir']
            scanInfo.([scanNames{iScan,1} 'DicomDir']){iSubject,iSes}=['/DSTORE/BRICIA/E182012_Brain_4D_Flow_RO/' ...
            temp_name(1:end) '/' direc(1).name '/' num2str(scanInfo.([scanNames{iScan,1} 'SeriesNo'])(iSubject,iSes)) ];
    end
    

end

cd(scanInfo.HVNumberStr{iSubject})
save scanInfo scanInfo % make sure in correct directory

