function [node, nodeNames, flag] = readSWMMvar
% based on:
% reads large TXT files, i.e. SWMM output files [node nodenames] = read_swmm;
% #fbl, Nov 2005, rev. May 2006, rev. Apr 2008 (skr,fbl)

flag = false;
wannaGetOrgResults = true;
[fileName, pathName] = uigetfile({'*.txt','SWMM result file (*.txt)'},'Choose File to Read');
if not(isequal(pathName,0) || isequal(fileName,0))
    if not(exist([pathName, fileName],'file'))
        error(['File "', pathName, fileName, '" not found (',mfilename,').']);
    end%if
else
    return;
end%if

filename = [pathName, fileName];
options.Resize='on';
instr = '0    % no quality';

answer = inputdlg('Enter number of quality parameters:','quality parameters', 1, {instr}, options);
if isempty(answer)
    return;
end%if
numberOfQualiParam = eval(answer{1});

%set SWMM typical conventions
%format of the TXT file is defined
format = '%s %f %f %f %f %f %f %f';
if numberOfQualiParam == 1
    format = '%s %f %f %f %f %f %f %f %f';
elseif numberOfQualiParam == 2
    format = '%s %f %f %f %f %f %f %f %f %f';
elseif numberOfQualiParam == 3
    format = '%s %f %f %f %f %f %f %f %f %f %f';
end%if

segcount = 200;

%opens the large TXT file for reading the #of nodes
file_id = fopen(filename, 'r');

%skip the 5 headerlines
nodenumbers = textscan(file_id, '%f','headerLines',numberOfQualiParam+5);
no_nodes = nodenumbers{1};
numberOfAllNodes = no_nodes;
segarray = textscan(file_id,'%s',no_nodes,'headerLines',1);
%closes the large TXT file
fclose(file_id);
nodeNames = segarray{1};
% Completes pending drawing events
drawnow;
pause(0.1);

% lists all nodes (nodenames) in the SWMM result file (*.txt)
[selection, ok] = listdlg('ListString', nodeNames, 'PromptString', 'Select results:', 'Name', 'Import SWMM data');
if not(ok)
    return;
end%if
tic;
no_nodes = length(selection);
nodeNames = nodeNames(selection);

%re-opens the large TXT file for reading
file_id = fopen(filename, 'r');
%reads 1st part of the txt-file skipping a '(no_nodes + 7)' headerlines
%reading starts @ no_nodes + 7 (depending on the number of headerlines)
%segarray = textscan(file_id,format,(no_nodes*1500),'headerlines',(no_nodes + 7));
segarray = textscan(file_id,format,(numberOfAllNodes),'headerlines',(numberOfQualiParam+numberOfAllNodes + 7));
k = 1;
m = 1;
for i = 1:no_nodes
    year = segarray{2}(numberOfAllNodes*(m-1)+selection(i));
    mon = segarray{3}(numberOfAllNodes*(m-1)+selection(i));
    day = segarray{4}(numberOfAllNodes*(m-1)+selection(i));
    hr = segarray{5}(numberOfAllNodes*(m-1)+selection(i));
    min = segarray{6}(numberOfAllNodes*(m-1)+selection(i));
    sec = segarray{7}(numberOfAllNodes*(m-1)+selection(i));
    seg{k}.flow(m,i) = segarray{8}(numberOfAllNodes*(m-1)+selection(i));
    %calculation of the corresponding numeric dates
    seg{k}.time(m,i) = datenum(year, mon, day, hr, min, sec);
end
if numberOfQualiParam > 0
    for i = 1:no_nodes
        seg{k}.cod(m,i) = segarray{9}(numberOfAllNodes*(m-1)+selection(i));
    end%for
    if numberOfQualiParam > 1
        for i = 1:no_nodes
            seg{k}.tkn(m,i) = segarray{10}(numberOfAllNodes*(m-1)+selection(i));
        end%for
    end%if
end%if

%splits large TXT file into 'segcount' segments à 12000 rows
for k = 2:segcount
    segarray = textscan(file_id,format,(numberOfAllNodes*1500));
    for m = 1:length(segarray{1,1})/numberOfAllNodes
        for i = 1:no_nodes
            year = segarray{2}(numberOfAllNodes*(m-1)+selection(i));
            mon = segarray{3}(numberOfAllNodes*(m-1)+selection(i));
            day = segarray{4}(numberOfAllNodes*(m-1)+selection(i));
            hr = segarray{5}(numberOfAllNodes*(m-1)+selection(i));
            min = segarray{6}(numberOfAllNodes*(m-1)+selection(i));
            sec = segarray{7}(numberOfAllNodes*(m-1)+selection(i));
            seg{k}.flow(m,i) = segarray{8}(numberOfAllNodes*(m-1)+selection(i));
            %calculation of the corresponding numeric dates
            seg{k}.time(m,i) = datenum(year, mon, day, hr, min, sec);
        end
    end
    if numberOfQualiParam > 0
        for m = 1:length(segarray{1,1})/numberOfAllNodes
            for i = 1:no_nodes
                seg{k}.cod(m,i) = segarray{9}(numberOfAllNodes*(m-1)+selection(i));
            end%for
        end%for
        if numberOfQualiParam > 1
            for m = 1:length(segarray{1,1})/numberOfAllNodes
                for i = 1:no_nodes
                    seg{k}.tkn(m,i) = segarray{10}(numberOfAllNodes*(m-1)+selection(i));
                end%for
            end%for
        end%if
    end%if
end
fclose(file_id);

% concatenate flow and time data (collected in vertically; generate 'nodes'=[time flow]
for i = 1:no_nodes
    node{i}.flow = [];
    node{i}.time = [];
    for k = 1:length(seg)
        node{i}.flow = vertcat(node{i}.flow, seg{k}.flow(:,i));
        node{i}.time = vertcat(node{i}.time, seg{k}.time(:,i));
    end
end
if numberOfQualiParam > 0
    for i = 1:no_nodes
        node{i}.cod = [];
        for k = 1:length(seg)
            node{i}.cod = vertcat(node{i}.cod, seg{k}.cod(:,i));
        end%for
    end%for
    if numberOfQualiParam > 1
        for i = 1:no_nodes
            node{i}.tkn = [];
            for k = 1:length(seg)
                node{i}.tkn = vertcat(node{i}.tkn, seg{k}.tkn(:,i));
            end%for
        end%for
    end%if
end%if

if wannaGetOrgResults
    %PLOTTING of 'node'
    figure
    for i = 1:no_nodes
        plot(node{1,i}.time,node{1,i}.flow,'Color',[1*(1/no_nodes*i) 0 1*(1/no_nodes*i)]);
        hold on;
    end

    ylabel('\bf Flow [Ls-1]','fontsize',14);
    grid on;
    legend(nodeNames);
    h = title(fileName);
    set(h,'Fontsize',14);
    set(gca,'Fontsize',14);
    dtick;
    save(filename(1:end-4), 'node', 'nodeNames');
else
    data = node{1,i}.time;
    for i = 1:no_nodes
        data(:,i+1) = node{1,i}.flow;
    end%for
    node = data;
    nodeNames = [{'Time'}, nodeNames'];
end%if
flag = true;
toc;