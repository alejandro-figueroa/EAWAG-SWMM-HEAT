% script generates a MAT variable 'itemList' in  the format of {element ids
% paramter} from a list of SWMM element IDs 'ids' | 
% The attribute 'element' can be 'link' or 'node'; 'ids' is the
% ID of the corresponding element; 'parameter' is the SWMM state
% variable for the corresponding element. | 'itemlist' is used in readswmmout.m and
% readOutItemlist.m; INPUT: 'ids'; OUTPUT: 'itemList.mat' fbl, 12/02/2012,
% rev: fbl, 26/04/2016

clear all;
% get list of SWMM 'elements' from XLS data (copy 'n paste from *.inp
% to XLS)
[fileName, pathName] = uigetfile({'*.xlsx'}, 'Choose XLSX file');
cd(pathName);
[num,ids,raw] = xlsread(fileName,'Sheet3');
clear num raw;

id=[ids(1:end) ];
type='link ';
item='flow';
itemlist = cell(length(id),1);

for nr=1:length(id)
    itemlist(nr,:)=strcat(type, {' '}, id(nr), {' '}, item);
end

% write MAT file with itemList
evalstr = strcat(['save ',[strcat('itemlist',num2str(nr),item)],' itemlist ids']);
eval(evalstr);
clearvars -except itemlist ids