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

datasets = ["New_AllDis", "Lfa1_WTDis", "Lfa1_KODis", "Wk1_Dis", ... 
    "Wk4_Dis", "Newd_13", "Newd_26"];

dataset = 'Wk1_Dis';

%num=xlsread('All_datasets_for_curate6.xlsx','Lfa1_WTDis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Lfa1_KODis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Wk1_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Wk4_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Newd_13');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Newd_26');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','PyTCR_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','OT1_Dis');%'AllData_Wk1');
num=xlsread('All_datasets_for_curate6.xlsx',dataset);%'AllData_Wk1');

Table=array2table(num, 'VariableNames',{...
    'Mouse'	'x'	'y' 'z'	'Unit'	'Category'...
    'Collection'	'TimeSlice'	'TrackID'...
    'ID'	'Time' 	'Dis'	' Timediff'	'Speed'	'CelID'	'Type'	'Fixed_time1' 'Fixed_time3' 'Fixed_time2' 'addT' 'NewID'});
  
%B) %%%%%%%%%%%% Liver Data P14 and OT1  
%num=xlsread('All_datasets_for_curate6.xlsx','PyTCR_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','OT1_Dis');%'AllData_Wk1');

%Table=array2table(num, 'VariableNames',{...
%    'ID'	'time' 'Fixed_time2' 'x' 'y' 'z'	'NewID' 'Dis'});
%    Table.Fixed_time2=Table.Fixed_time2.*60;


%C) %%%%%%%%%%%%% LN data
%num=xlsread('0409_Doc5_P14',' Position(2)');%'AllData_Wk1');
%num=xlsread('0409_Doc6_P14','Position(2)');%'AllData_Wk1');
%num=xlsread('0409_Doc7_P14','Position(2)');%'AllData_Wk1');
%num=xlsread('0209_Doc6_P14','Position(2)');%'AllData_Wk1');
%num=xlsread('0209_Doc7_P14','Position(2)');%'AllData_Wk1');
%num=xlsread('0209_Doc8_P14','Position(2)');%'AllData_Wk1');

%Table=array2table(num, 'VariableNames',{...    
%'x' 'y' 'z'	'time' 'Fixed_time2' 'NewID2x'...
%'TrackID'	'NewID'});
%    Table.Fixed_time2=Table.Fixed_time2*20;


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


clear i X dX time X0
ma = msdanalyzer(3, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);
ma = ma.computeMSD;

% VVG addition
% do log-log regression and estimate slope
fit=ma.fitMeanMSDloglog(0.1);  % taking only 10% of the data points


ma.plotMeanMSD(gca, false);
% adding extra to the graph
rl=refline([fit.p1 fit.p2])
rl.Color = 'b';
rl.LineWidth=1.5
rl.LineStyle = '--';
ax = gca; % current axes
ax.FontSize = 12;
%ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0 4.5];
ax.XLim = [0.5 4];
legend([erase(dataset,["_" "Dis"]) "\gamma="+sprintf('%.2f',fit.p1)], ...
    'FontSize',14,'Location','northwest')

%set(leg,'Interpreter', 'none')
%text(0.5,3.5,['\gamma=' fit.p1],'FontSize',14)


