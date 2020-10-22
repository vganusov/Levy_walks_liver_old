clear all; clf
%rng shuffle
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Note that 
%LfaWT; KO; w1k1; 1k2; AllneW1326: 27.5; 27.5;  6.51; 3.25; 26; secs
%P14; OT1; 1 min
%LN 20 sec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data (open the relevant data links): 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OT1

dataset1 = 'OT1_Dis';
    
%B) %%%%%%%%%%%% Liver Data P14 and OT1  
num=xlsread('All_datasets_for_curate6.xlsx',dataset1);%'AllData_Wk1');

Table=array2table(num, 'VariableNames',{...
    'ID'	'time' 'Fixed_time2' 'x' 'y' 'z'	'NewID' 'Dis'});
    Table.Fixed_time2=Table.Fixed_time2.*60;


%%%%%%%%%%%%%%%%%%%%%%%%%%
SPACE_UNITS = 'µm';
TIME_UNITS = 'sec';
UA=unique(Table.NewID);
[co xx]= size(UA);
N_PARTICLES=co; 
N_DIM = 3;

dT=Table.Fixed_time2(1);

SIZE = 3; % µm
tracks = cell(N_PARTICLES, 1);
for i = 1 : N_PARTICLES
    TT=(Table.NewID==UA(i));
    x=Table.x(TT);y=Table.y(TT); 
    z=Table.z(TT);
    X=[x y z];
    time=cumsum(Table.Fixed_time2(TT));
    tracks{i} = [time X];
end

maOT1 = msdanalyzer(3, SPACE_UNITS, TIME_UNITS);
maOT1 = maOT1.addAll(tracks);
maOT1 = maOT1.computeMSD;

% do log-log regression and estimate slope
fit1=maOT1.fitMeanMSDloglog(0.1);  % taking only 10% of the data points
p1=maOT1.plotMeanMSD(gca, false);
p1.Marker = '>';
p1.Color = 'k';
ax = gca; % current axes
ax.FontSize = 12;
%ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0 4.5];
ax.XLim = [0.5 4];
rl=refline([fit1.p1 fit1.p2])
rl.Color = 'k';
rl.LineWidth=1.5
rl.LineStyle = '--';
hold on

dataset2 = 'PyTCR_Dis';
    
%B) %%%%%%%%%%%% Liver Data P14 and OT1  
num=xlsread('All_datasets_for_curate6.xlsx',dataset2);%'AllData_Wk1');

Table=array2table(num, 'VariableNames',{...
    'ID'	'time' 'Fixed_time2' 'x' 'y' 'z'	'NewID' 'Dis'});
    Table.Fixed_time2=Table.Fixed_time2.*60;


%%%%%%%%%%%%%%%%%%%%%%%%%%
SPACE_UNITS = 'µm';
TIME_UNITS = 'sec';
UA=unique(Table.NewID);
[co xx]= size(UA);
N_PARTICLES=co; 
N_DIM = 3;

dT=Table.Fixed_time2(1);

SIZE = 3; % µm
tracks = cell(N_PARTICLES, 1);
for i = 1 : N_PARTICLES
    TT=(Table.NewID==UA(i));
    x=Table.x(TT);y=Table.y(TT); 
    z=Table.z(TT);
    X=[x y z];
    time=cumsum(Table.Fixed_time2(TT));
    tracks{i} = [time X];
end


maP14 = msdanalyzer(3, SPACE_UNITS, TIME_UNITS);
maP14 = maP14.addAll(tracks);
maP14 = maP14.computeMSD;

% do log-log regression and estimate slope
fit2=maP14.fitMeanMSDloglog(0.1);  % taking only 10% of the data points
p2=maP14.plotMeanMSD(gca, false);
p2.Marker = '*';
p2.Color = 'b';
rl=refline([fit2.p1 fit2.p2])
rl.Color = 'b';
rl.LineWidth=1.5
rl.LineStyle = '--';
hold on

legend([p1 p2], {"OT1"+ " (\gamma="+sprintf('%.2f',fit1.p1)+")", ...
    "P14"+ " (\gamma="+sprintf('%.2f',fit2.p1)+")"}, ... 
    'FontSize',14,'Location','northwest')

%set(leg,'Interpreter', 'none')
%text(0.5,3.5,['\gamma=' fit.p1],'FontSize',14)



% just making the graph
clf
fit1=maOT1.fitMeanMSDloglog(0.1);  % taking only 10% of the data points
p1=maOT1.plotMeanMSD(gca, false);
p1.Marker = '>';
p1.Color = 'k';
ax = gca; % current axes
ax.FontSize = 14;
%ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0 4.5];
ax.XLim = [0.5 4];
rl=refline([fit1.p1 fit1.p2])
rl.Color = 'k';
rl.LineWidth=1
rl.LineStyle = '--';
hold on

% do log-log regression and estimate slope
fit2=maP14.fitMeanMSDloglog(0.1);  % taking only 10% of the data points
p2=maP14.plotMeanMSD(gca, false);
p2.Marker = '*';
p2.Color = 'b';
rl=refline([fit2.p1 fit2.p2])
rl.Color = 'b';
rl.LineWidth=1
rl.LineStyle = '--';
hold on

xtick=0:1:5;
xticklab = cellstr(num2str(xtick(:), '10^{%d}'));
set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex')
ytick=0:1:5;
yticklab = cellstr(num2str(ytick(:), '10^{%d}'));
set(gca,'YTick',ytick,'YTickLabel',yticklab,'TickLabelInterpreter','tex')
xlabel('Delay (sec)')
ylabel('MSD (\mum^2)')

legend([p1 p2], {"OT1"+ " (\gamma="+sprintf('%.2f',fit1.p1)+")", ...
    "P14"+ " (\gamma="+sprintf('%.2f',fit2.p1)+")"}, ... 
    'FontSize',14,'Location','northwest')



% yt = get(gca,'ytick');
% for j=1:length(yt)
%     % With log plots, MATLAB defaulted to exponential format, that is difficult for lay
%     % readerst to understand. In my case, I wanted integer format.
%     YTL{1,j} = num2str(yt(j),'%d');
% end
% yticklabels(YTL);


