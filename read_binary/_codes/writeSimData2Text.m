% write SWMM simulation results to TXT data (4Alex)

% run ReadSwmmOutTest.m before 

% itemlist = {'node 138a_Mess-Sch. depth'; 'node 58_Spezialsch._ depth'; 'node RK59_MW depth';...
%     'link Link_2 flow'; 'node RB_59 depth'; 'node RÜB+PW80_Industrie depth'; 'node RK59_Deckel01 depth'}; 
% ids = {'RKB Usterstr. (B60)'; 'Überlaufschacht 58 (B5F)'; 'Schacht Zulauf ARA (n257)';...
%     'Zulauf ARA Ls-1 (Qin)'; 'RKB ARA (n256)'; 'RKB Industrie (W769)'; 'RK59_Deckel01 (SK ARA)'};

%% intitialise
StartTimeOffset = datenum('30-Dec-1899')+d.starttime;
formatSpec = '%s\t %5.1f \n';

%% B60 - h in mm @ RKB Usterstr
wData = [res(:,1)+StartTimeOffset res(:,2).*1000];

for rowcount = 1:length(wData)
    celldata(rowcount,:) = {datestr(wData(rowcount,1),'dd-mmm-yyyy HH:MM:SS'), wData(rowcount,2)};
end

fileID = fopen('B60sim_h_mm.txt','w');
% fprintf(fileID,'%s \n',strcat('h [mm] @ ',ids{1}));

[nrows,ncols] = size(celldata);

for row = 1:nrows
    fprintf(fileID,formatSpec,celldata{row,:});
end
fclose(fileID);
clear celldata

%% B5F - h in mm @ Schieber Zürichstr
wData = [res(:,1)+StartTimeOffset res(:,3).*1000];

for rowcount = 1:length(wData)
    celldata(rowcount,:) = {datestr(wData(rowcount,1),'dd-mmm-yyyy HH:MM:SS'), wData(rowcount,2)};
end
fileID = fopen('B5Fsim_h_mm.txt','w');
[nrows,ncols] = size(celldata);
for row = 1:nrows
    fprintf(fileID,formatSpec,celldata{row,:});
end
fclose(fileID);
clear celldata

%% n257 - h in mm @ Zulaufsammler ARA
wData = [res(:,1)+StartTimeOffset res(:,4).*1000];

for rowcount = 1:length(wData)
    celldata(rowcount,:) = {datestr(wData(rowcount,1),'dd-mmm-yyyy HH:MM:SS'), wData(rowcount,2)};
end
fileID = fopen('n257sim_h_mm.txt','w');
[nrows,ncols] = size(celldata);
for row = 1:nrows
    fprintf(fileID,formatSpec,celldata{row,:});
end
fclose(fileID);
clear celldata

%% n256 / n258 - h in mm @ RKB ARA
wData = [res(:,1)+StartTimeOffset res(:,6).*1000];

for rowcount = 1:length(wData)
    celldata(rowcount,:) = {datestr(wData(rowcount,1),'dd-mmm-yyyy HH:MM:SS'), wData(rowcount,2)};
end
fileID = fopen('n256sim_h_mm.txt','w');
[nrows,ncols] = size(celldata);
for row = 1:nrows
    fprintf(fileID,formatSpec,celldata{row,:});
end
fclose(fileID);
clear celldata

%% W769 - h in mm @ PW Industrie
wData = [res(:,1)+StartTimeOffset res(:,7).*1000];

for rowcount = 1:length(wData)
    celldata(rowcount,:) = {datestr(wData(rowcount,1),'dd-mmm-yyyy HH:MM:SS'), wData(rowcount,2)};
end
fileID = fopen('W769sim_h_mm.txt','w');
[nrows,ncols] = size(celldata);
for row = 1:nrows
    fprintf(fileID,formatSpec,celldata{row,:});
end
fclose(fileID);
clear celldata
