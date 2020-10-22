% Creates efficiency -time- graphs for 1D-structure, 
% bothe Brownian (mu=4) and Levy steps (mu=2.5) 

clear all; clf

dist=500; %%fixed distance
simu=3000; %number of cells
%nsteps=3000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pareto-Levy
%alph=1.5;
%k_gp=1/alph;
%xmin=2.1;
%sig_gp=xmin/alph; 


mean_R=2.1;  % mean step length
alphR=1.5;%3.5;%alpha=4.5;  1
alphR_LW=alphR;%3.5;%alpha=4.5;  1
k_gpR=1/alphR;
xminR=mean_R*((alphR-1)/alphR);
sig_gpR=xminR/alphR; 
mean_PaR=alphR*xminR/(alphR-1);
sig_GR=mean_PaR*sqrt(pi()/2);
mean_GuR=sig_GR*sqrt(2/pi());
max_stepR=300;


for i=1:simu %number of simulated cells
d=0;timeTOd=0;

while d<=dist                                                           
%st_stepln_Pareto1=gprnd(k_gp,sig_gp,sig_gp/k_gp,1,1);
%st_stepln_Pareto=min(st_stepln_Pareto1,300);


%brlen1R=gprnd(k_gpR,sig_gpR,sig_gpR/k_gpR,1);
%step=min(brlen1R,max_stepR);

st_stepln_Pareto12=gprnd(k_gpR,sig_gpR,sig_gpR/k_gpR,1);
st_stepln_Pareto2=min(st_stepln_Pareto12,max_stepR);


d=d+st_stepln_Pareto2;
timeTOd=timeTOd+1;
end
Par_st_timeTOd(i)=timeTOd;
end                                     
mean_LW=mean(Par_st_timeTOd)
median_LW=median(Par_st_timeTOd)
var_LW=var(Par_st_timeTOd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Gaussian
%mean_Pa=alph*xmin/(alph-1)
%sig_G=mean_Pa*sqrt(pi()/2);

mean_R=2.1;
alphR=4.5;
alphR_BW=alphR;%3.5;%alpha=4.5;  1
k_gpR=1/alphR;
xminR=mean_R*((alphR-1)/alphR);
sig_gpR=xminR/alphR; 
mean_PaR=alphR*xminR/(alphR-1);
sig_GR=mean_PaR*sqrt(pi()/2);
mean_GuR=sig_GR*sqrt(2/pi());
max_stepR=300;

for i=1:simu %number of simulated cells
    %i
d=0;timeTOd=0;
while d<=dist
%st_stepln_Gaussian=random('HalfNormal',0,sig_G, 1,1); 

%brlen1R=gprnd(k_gpR,sig_gpR,sig_gpR/k_gpR,1);
%step=min(brlen1R,max_stepR);

st_stepln_Pareto1=gprnd(k_gpR,sig_gpR,sig_gpR/k_gpR,1);
st_stepln_Pareto=min(st_stepln_Pareto1,max_stepR);


d=d+st_stepln_Pareto;
timeTOd=timeTOd+1;
end
Gau_st_timeTOd(i)=timeTOd;
end
mean_Gu=sig_GR*sqrt(2/pi());
mean_CBW=mean(Gau_st_timeTOd)
median_CBW=median(Gau_st_timeTOd)
var_CBW=var(Gau_st_timeTOd)


save Simulations_data_1D-LW-25-BW-45


% exporting data into csv
%dist=500; %%fixed distance
%simu=3000; %number of cells
exportdata_1D=transpose([Gau_st_timeTOd;Par_st_timeTOd]);
csvwrite("simulation_data-1D-BW-45-LW-25-3000.csv",exportdata_1D)


% loading simulation results (n=3000 cells)
clear all
load Simulations_data_1D-LW-25-BW-45
%s

% plotting the figure
dt=0.03;
figure(); 
hold on
h=histogram(log10(Gau_st_timeTOd),'Normalization','probability','BinWidth',dt);
%h.EdgeColor=[0 0 1]
h.EdgeColor='k'
h.FaceColor=[1 1 1]
h.LineWidth=2
h.LineStyle='-.'
%    'EdgeColor',[0.4660 0.6740 0.1880], 'FaceColor', [0.4940 0.1840 0.5560], 'LineWidth', 2);
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'log')
set(gca,'FontSize',14)
xlim([1 3]);ylim([0.001 1])

h=histogram(log10(Par_st_timeTOd),'Normalization','probability','BinWidth',dt);%,...
%    'EdgeColor',[0.4660 0.6740 0.1880], 'FaceColor', [0.4940 0.1840 0.5560], 'LineWidth', 2);
h.EdgeColor=[1 0 0]
h.FaceColor=[1 1 1]
h.EdgeAlpha=0.5
h.LineWidth=2
h.LineStyle='-'

legend('BW (\mu=4.5): T='+compose("%3.0f",round(mean_CBW)),'LW (\mu=2.5) T='+compose("%3.0f",round(mean_LW)),'Location','northwest')
xlabel('Time to reach a target at a fixed distance (log_{10} steps)');
ylabel('Probability')
hold off



% calculating time to target for different number of T cells searching
% for the target


ncell=100;
min(randsample(Par_st_timeTOd,ncell,true)) % example

LW_T0d_cells=zeros(simu,1);
BW_T0d_cells=zeros(simu,1);

rng(1)   % setting seed for reproducing results later
for i = 1:simu
    LW_T0d_cells(i,1)= min(randsample(Par_st_timeTOd,ncell,true));
    BW_T0d_cells(i,1)= min(randsample(Gau_st_timeTOd,ncell,true));
end

mean(Par_st_timeTOd)  % search is done by 1 cell with LW
mean(Gau_st_timeTOd)  % search is done by 1 cell with BW
mean(LW_T0d_cells)  % search is done by ncell cells with LW
mean(BW_T0d_cells)  % search is done by ncell cells with BW


% varying the number of cells searching for infection
LW_T0d_cells=zeros(simu,1);
BW_T0d_cells=zeros(simu,1);

ncells = [1 2 5 10 30 100 300 1000]
LW_meanT=zeros(length(ncells),3)
BW_meanT=zeros(length(ncells),3)

rng(2)   % setting seed for reproducing results later
for n=1:length(ncells)
    ncell=ncells(n)
    for i = 1:simu
        LW_T0d_cells(i,1)= min(randsample(Par_st_timeTOd,ncell,true));
        BW_T0d_cells(i,1)= min(randsample(Gau_st_timeTOd,ncell,true));
    end
    LW_meanT(n,1) = mean(LW_T0d_cells);  % search is done by ncell cells with LW
    LW_meanT(n,2:3) = quantile(LW_T0d_cells,[0.025 0.975]);  % 95%CIs
    BW_meanT(n) = mean(BW_T0d_cells);  % search is done by ncell cells with BW
    BW_meanT(n,2:3) = quantile(BW_T0d_cells,[0.025 0.975]);  % 95%CIs
end

% Plot fit with data: linear scale for time
figure()
hold on
%errorbar(x,y,yneg,ypos,xneg,xpos,'o')
errorbar(ncells,BW_meanT(:,1), BW_meanT(:,1)-BW_meanT(:,2), BW_meanT(:,3)-BW_meanT(:,1),'-ok','LineWidth',1.5);
%plot(ncells,BW_meanT(:,1), '-ok','LineWidth',1.5);
%plot(ncells,LW_meanT(:,1), '--^r','LineWidth',1.5);
errorbar(ncells,LW_meanT(:,1), LW_meanT(:,1)-LW_meanT(:,2), LW_meanT(:,3)-LW_meanT(:,1),'--^r','LineWidth',1.5);

set(gca,'XScale','log')
set(gca,'YScale','linear')
set(gca,'FontSize',14)
pos = get(gca, 'Position');
pos(2) = 0.05;
pos(3) = 0.05;
%set(gca, 'Position', pos)
axis([0.9,1000,0,350])
legend('BW (\mu=4.5)', 'LW (\mu=2.5)','LineWidth',1.5);
% Label axes
xlabel('Number of cells searching for a single target')
ylabel('Mean time to target in the liver (steps)')

% end of plot


% Plot fit with data: log scale for time

BW_meanT2=log10(BW_meanT);
LW_meanT2=log10(LW_meanT);

figure()
hold on
errorbar(ncells,BW_meanT2(:,1), BW_meanT2(:,1)-BW_meanT2(:,2), BW_meanT2(:,3)-BW_meanT2(:,1),'-ok','LineWidth',1.5);
errorbar(ncells,LW_meanT2(:,1), LW_meanT2(:,1)-LW_meanT2(:,2), LW_meanT2(:,3)-LW_meanT2(:,1),'--^r','LineWidth',1.5);

set(gca,'XScale','log')
set(gca,'YScale','linear')
set(gca,'FontSize',18)
pos = get(gca, 'Position');
pos(2) = 0.05;
pos(3) = 0.05;
%set(gca, 'Position', pos)
%axis([0.9,1000,0,4])
axis([0.9,1000,0,4])
legend('BW (\mu=4.5)', 'LW (\mu=2.5)','LineWidth',1.5);
% Label axes
xlabel('Number of cells searching for a single target')
ylabel('Mean time to target in liver (log_{10} steps)')

% end of plot



