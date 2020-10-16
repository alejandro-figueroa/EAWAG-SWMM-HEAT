%% USER INPUT
% path were the SWMM-TEMP files are
clear all
clc

%%
cd('C:\Users\figueral\Documents\SWMM-TEMP\PacketEAWAG\PacketEAWAG\test_case_pulse\new_swmm_temp\')
simpath = 'C:\Users\figueral\Documents\SWMM-TEMP\PacketEAWAG\PacketEAWAG\test_case_pulse\new_swmm_temp\';
resultpath = 'C:\Users\figueral\Documents\SWMM-TEMP\PacketEAWAG\PacketEAWAG\test_case_pulse\new_swmm_temp\';
runpath = 'C:\Users\figueral\Documents\SWMM-TEMP\PacketEAWAG\PacketEAWAG\test_case_pulse\new_swmm_temp\_codes';
%%addpath('C:\Users\figueral\Documents\SWMM-TEMP\PacketEAWAG\PacketEAWAG\test_case_pulse\new_swmm_temp\_readSwmmOutShell\_codes');


%%
a=readtable('Temp_Sens.txt');
b=readtable('Flow_Sens.txt');

a.Var2(:)=2;
a.Var2(5000:7000)=2;


b.Var2(:)=20;
b.Var2(5000:7000)=20;

writetable(a,'Flow_Sens','FileType','text','WriteVariableNames',0,'Delimiter','tab')
writetable(b,'Temp_Sens','FileType','text','WriteVariableNames',0,'Delimiter','tab')

%% Launch SWMM from commandline
name='Winter';
filename_in = [simpath,sprintf('%s.inp',name)];
filename_out = [simpath,sprintf('%s.out',name)];
filename_report = [simpath,sprintf('%s.rpt',name)];
[status,cmdout] = system([runpath,'swmm5_temp.exe ',filename_in,' ',filename_report,' ',filename_out,' ']);

%%
[d] =  readswmmout3('open',filename_out);
% pre-define period for which data are relevant
fromidx = 1;
toidx = d.Nperiods; %d.Nperiods is the maxiumum number of periods that can be read
% get all data
% [err,res] =  readswmmout('get',d,[fromidx,toidx]);
% load OR generate 'itemlist'
itemlist = {'node Out inflow';'link Out_link air_velocity'; 'link 591_592 air_velocity'; 'node Out wtemperature'}; ids = {'Out'};
%itemlist = {'node 592 inflow';'node Out inflow';'node 592 temperature';'node Out temperature'}; ids = {'592';'Out'};
% read data for elements and attributes as specified in 'itemlist'
flowunit = 'lps';
[err,res,unitlist]=  readswmmout3('getitems',d,[fromidx,toidx],itemlist,flowunit); 
% create startTimeOffset due to deviation in readswmmout3.m
StartTimeOffset = datenum('30-Dec-1899')+d.starttime;
simRes = [res(:,1)+StartTimeOffset res(:,2:end)];

    %%
    simFlow = table(simRes(:,1),simRes(:,2));
    simFlow.Properties.VariableNames = {'time','value'};
    simFlow.time=datetime(simFlow.time,'ConvertFrom','datenum','Format','MMM-dd-yyyy HH:mm');
    simFlow=table2timetable(simFlow);
    %a=~isnan(simFlow.value);
    %simFlow=simFlow(a,:);
    plot(simFlow.time,simFlow.value)
    hold on
    simVel1 = table(simRes(:,1),simRes(:,3));
    simVel1.Properties.VariableNames = {'time','value'};
    simVel1.time=datetime(simVel1.time,'ConvertFrom','datenum','Format','MMM-dd-yyyy HH:mm');
    simVel1=table2timetable(simVel1);
    %b=~isnan(simTemp.value);
    %simTemp=simTemp(b,:);
    plot(simVel1.time,simVel1.value)
    hold on
    simVel2 = table(simRes(:,1),simRes(:,4));
    simVel2.Properties.VariableNames = {'time','value'};
    simVel2.time=datetime(simVel2.time,'ConvertFrom','datenum','Format','MMM-dd-yyyy HH:mm');
    simVel2=table2timetable(simVel2);
    %b=~isnan(simTemp.value);
    %simTemp=simTemp(b,:);
    plot(simVel2.time,simVel2.value)
    hold on
    simTemp = table(simRes(:,1),simRes(:,5));
    simTemp.Properties.VariableNames = {'time','value'};
    simTemp.time=datetime(simTemp.time,'ConvertFrom','datenum','Format','MMM-dd-yyyy HH:mm');
    simTemp=table2timetable(simTemp);
    %b=~isnan(simTemp.value);
    %simTemp=simTemp(b,:);
    plot(simTemp.time,simTemp.value)



%%
cd('C:\Users\joshipra\polybox\Master thesis\Codes')
[NSE_SWMM_TEMP,RMSE_SWMM_TEMP,PBIAS_SWMM_TEMP]=gof(simTemp.value(2:end),obsTemp.value(2:end));

%%
figure('units','normalized','outerposition',[0 0 1 1])
plot(simTemp.time,simTemp.value,'DisplayName','Simulated')
hold on
plot(obsTemp.time,obsTemp.value,'DisplayName','Observed')
hold off
legend('FontSize',20)
ylim([10 18])
str={['NSE (SWMM_TEMP) = ',num2str(NSE_SWMM_TEMP,'%.2f')];['RMSE (SWMM_TEMP) = ',num2str(RMSE_SWMM_TEMP,'%.2f')]};
annotation('textbox',[.14 .82 .1 .1],'string',str,'FitBoxToText','on','FontSize',20,'Interpreter','none','EdgeColor','none')
ylabel('Wastewater temperature [^o C]','FontSize',20)
set(gca,'FontSize',20)

cd('C:\Users\joshipra\polybox\Master thesis\Comparison')
saveas(gcf,'FehraltorfSc1_SWMM_NoSoilConductance', 'png')

%%
% store data
upTemp = table(simRes(:,1),simRes(:,4));
upTemp.Properties.VariableNames = {'time','value'};
sewerTemp = table(simRes(:,1),simRes(:,5));
sewerTemp.Properties.VariableNames = {'time','value'};

upFlow = table(simRes(:,1),simRes(:,2));
upFlow.Properties.VariableNames = {'time','value'};
sewerFlow = table(simRes(:,1),simRes(:,3));
sewerFlow.Properties.VariableNames = {'time','value'};

%%
figure('units','normalized','outerposition',[0 0 1 1])
yyaxis left
plot(datetime(upFlow.time,'ConvertFrom','datenum'),upFlow.value,'k-','LineWidth',2,'DisplayName','Upstream flow')
hold on
plot(datetime(sewerFlow.time,'ConvertFrom','datenum'),sewerFlow.value,'b--','LineWidth',2,'DisplayName','Downstream flow')
hold off
legend('FontSize',20)
ylim([20 60])
ylabel('Wastewater flow rate [L/s]','FontSize',20)
set(gca,'FontSize',20);
xlabel('FontSize',20)


%figure('units','normalized','outerposition',[0 0 1 1])
yyaxis right
plot(datetime(upTemp.time,'ConvertFrom','datenum'),upTemp.value,'r-','LineWidth',2,'DisplayName','Upstream temp')
hold on
plot(datetime(sewerTemp.time,'ConvertFrom','datenum'),sewerTemp.value,'g--','LineWidth',2,'DisplayName','Downstream temp')
hold off
legend('FontSize',20)
ylim([8 55])
ylabel('Wastewater temperature [^o C]','FontSize',20)
%set(gca,'FontSize',20);
%xlabel('FontSize',20)

Tdiff=mean(sewerTemp.value)-mean(upTemp.value);
Tdiff_pc=(mean(sewerTemp.value)-mean(upTemp.value))/mean(upTemp.value)*100;

%%
figure('units','normalized','outerposition',[0 0 1 1])
yyaxis left
plot(datetime(upFlow.time(2:end),'ConvertFrom','datenum'),upFlow.value(2:end),'LineWidth',2,'DisplayName','Upstream flow')
hold on
plot(datetime(sewerFlow.time(2:end),'ConvertFrom','datenum'),sewerFlow.value(2:end),'LineWidth',2,'DisplayName','Downstream flow')
hold off
legend('FontSize',30,'Location','Northwest')
ylim([20 60])
ylabel('Wastewater flow rate [L/s]','FontSize',30)
set(gca,'FontSize',30);
%xlabel('FontSize',20)


%figure('units','normalized','outerposition',[0 0 1 1])
yyaxis right
plot(datetime(upTemp.time(2:end),'ConvertFrom','datenum'),upTemp.value(2:end),'LineWidth',2,'DisplayName','Upstream temp')
hold on
plot(datetime(sewerTemp.time(2:end),'ConvertFrom','datenum'),sewerTemp.value(2:end),'LineWidth',2,'DisplayName','Downstream temp')
hold off
legend('FontSize',30,'Location','Northwest')
ylim([min(upTemp.value)-2 max(upTemp.value)+2])
ylabel('Wastewater temperature [^oC]','FontSize',30)
%set(gca,'FontSize',20);

%%
%xlabel('FontSize',20)
saveas(gcf, 'Noded_Sc4_ZeroAirVelocity', 'png')

%%
% read the observed and input values
obsTemp_Cal=readtable('Cal_TempDown_Obs_SWMM.txt');
obsFlow_Cal=readtable('Cal_FlowDown_Obs_SWMM.txt');

obsTemp_Val=readtable('Summer_TempDown_Obs_SWMM.txt');
obsFlow_Val=readtable('Summer_FlowDown_Obs_SWMM.txt');

Temp_up=readtable('Summer_TempUp_Obs_SWMM.txt');
Flow_up=readtable('Summer_FlowUp_Obs_SWMM.txt');

obsTemp_Cal.Properties.VariableNames={'time','value'};
obsFlow_Cal.Properties.VariableNames={'time','value'};
obsFlow_Val.Properties.VariableNames={'time','value'};
obsTemp_Val.Properties.VariableNames={'time','value'};
Temp_up.Properties.VariableNames={'time','value'};
Flow_up.Properties.VariableNames={'time','value'};

%%
%evalStart=datenum('08-Apr-2019 00:00:00');
evalEnd=datenum('11-Apr-2019 11:18:20');

evalStart=datenum('24-Jun-2019 00:00:00');
evalEnd=datenum('27-Jun-2019 11:18:20');


i = find(evalStart<=datenum(sewerTemp.time) & evalEnd>=datenum(sewerTemp.time));
j = find(evalStart<=datenum(obsTemp_Val.time) & evalEnd>=datenum(obsTemp_Val.time));
k = find(evalStart<=datenum(Temp_up.time) & evalEnd>=datenum(Temp_up.time));

l = find(evalStart<=datenum(sewerFlow.time) & evalEnd>=datenum(sewerFlow.time));
m = find(evalStart<=datenum(obsFlow_Val.time) & evalEnd>=datenum(obsFlow_Val.time));
n = find(evalStart<=datenum(Flow_up.time) & evalEnd>=datenum(Flow_up.time));


obsTemp=obsTemp_Val(j,:);
simTemp=sewerTemp(i,:);
Temp_up=Temp_up(k,:);

obsFlow=obsFlow_Val(m,:);
simFlow=sewerFlow(l,:);
Flow_up=Flow_up(n,:);

formatOut='dd-MMM-yyyy HH:mm:ss';
simTemp.time=datetime(simTemp.time,'ConvertFrom','datenum','Format',formatOut);
obsTemp.time=datetime(obsTemp.time,'ConvertFrom','datenum','Format',formatOut);
Temp_up.time=datetime(Temp_up.time,'ConvertFrom','datenum','Format',formatOut);

simFlow.time=datetime(simFlow.time,'ConvertFrom','datenum','Format',formatOut);
obsFlow.time=datetime(obsFlow.time,'ConvertFrom','datenum','Format',formatOut);
Flow_up.time=datetime(Flow_up.time,'ConvertFrom','datenum','Format',formatOut);

%%
cd('C:\Users\joshipra\polybox\Master thesis\Codes')

[NSE_temp,RMSE_temp,PBIAS_temp]=gof(obsTemp.value,simTemp.value);
[NSE_flow,RMSE_flow,PBIAS_flow]=gof(obsFlow.value/1000,simFlow.value/1000);


%%
a=innerjoin(Temp_up,innerjoin(simTemp,obsTemp,'Keys','time'),'Keys','time');
figure('units','normalized','outerposition',[0 0 1 1])
plot(a.time,a.value,'ko','MarkerSize',3,'DisplayName','Temp_{up,obs}')
hold on
plot(a.time,a.value_simTemp,'b.-','DisplayName','Temp_{down,sim}')
plot(a.time,a.value_obsTemp,'ro','MarkerSize',3,'DisplayName','Temp_{down,obs}')
legend()
ylabel('Water temperature [^oC]')
ylim([8 25])
grid on
str={['NSE = ',num2str(NSE_temp)];['RMSE = ',num2str(RMSE_temp)]};
annotation('textbox',[.14 .80 .1 .1],'string',str,'FitBoxToText','on','FontSize',30)
set(gcf, 'Color', 'w');
set(gca,'FontSize',30);
grid on

%%
cd('C:\Users\joshipra\polybox\Master thesis\Pictures')
saveas(gcf, 'Val_temp_SWMMtemp.fig');
export_fig Val_temp_SWMMtemp.png;

%%
cd('C:\Users\joshipra\polybox\Master thesis\Comparison')
name2='Summer';
writetimetable(simTemp,sprintf('%sT_SWMM_new.txt',name2),'delimiter','\t','WriteVariableNames',0)
writetimetable(simFlow,sprintf('%sF_SWMM_new.txt',name2),'delimiter','\t','WriteVariableNames',0)


%%
b=innerjoin(Flow_up,innerjoin(simFlow,obsFlow,'Keys','time'),'Keys','time');
figure('units','normalized','outerposition',[0 0 1 1])
plot(b.time,b.value/1000,'ko','MarkerSize',3,'DisplayName','Flow_{up,obs}')
hold on
plot(b.time,b.value_simFlow/1000,'b.-','DisplayName','Flow_{down,sim}')
plot(b.time,b.value_obsFlow/1000,'ro','MarkerSize',3,'DisplayName','Flow_{down,obs}')
legend()
ylabel('Flow [m^3/s]')
%ylim([0 0.10])
grid on
str={['NSE = ',num2str(NSE_flow)];['RMSE = ',num2str(RMSE_flow)]};
annotation('textbox',[.14 .80 .1 .1],'string',str,'FitBoxToText','on','FontSize',30)
set(gcf, 'Color', 'w');
set(gca,'FontSize',30);

%%
cd('C:\Users\joshipra\polybox\Master thesis\Pictures')
saveas(gcf, 'Cal_flow_SWMMtemp.fig');
export_fig Cal_flow_SWMMtemp.png;


%%
figure('name','raw simulation results');
subplot(4,1,1)
hold on
title('Downstream')
plot(sewerTemp.time,sewerTemp.value,'-','Color',[0.2 0.2 0.2],'MarkerSize',1)
%legend({'555 data','555 sim'});
%ylim([10 25])
dtick; grid on;
ylabel(strcat('Temperature [°C]'),'Fontsize',12);
set(gca,'Fontsize',12);


%%
% [NSE_flow,RMSE_flow,PBIAS_flow]=gof(obs_flow(1:end)./1000,flow_5min);
[NSE_temp,RMSE_temp,PBIAS_temp]=gof(obsFlow_Cal(1:4999),sewerFlow.value);
peakT_obs=max(obs_temp.Tw(1:4999));
peakT_sim=max(sim.Tw);

[NSE_flow,RMSE_flow,PBIAS_flow]=gof(obs_flow.Qw(1:4999)./1000,sim.Qw);
peakQ_obs=max(obs_flow.Qw(1:4999)./1000);
peakQ_sim=max(sim.Qw);

%%
plot(Temp_up.value,'ko')
hold on
plot(simTemp.value,'b')
plot(obsTemp.value,'ro')