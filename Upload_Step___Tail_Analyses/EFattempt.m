clear all; 

%WK1

dt=0.1;
x=0.01:dt:100;

% A) %%%%%%%% Liver Data  
%num=xlsread('All_datasets_for_curate6.xlsx','Lfa1_WTDis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Lfa1_KODis');%'AllData_Wk1');
num=xlsread('All_datasets_for_curate6.xlsx','Wk1_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Wk4_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','New_AllDisT');%'AllData_Wk1');

Table=array2table(num, 'VariableNames',{...
    'Mouse'	'x'	'y' 'z'	'Unit'	'Category'...
    'Collection'	'TimeSlice'	'TrackID'...
    'ID'	'Time' 	'Dis'	' Timediff'	'Speed'	'CelID'	'Type'	'Fixed_time1' 'Fixed_time3' 'Fixed_time2' 'addT' 'NewID'});
  Table.Fixed_time2=Table.Fixed_time2;  
 
UA=unique(Table.NewID);
[co xx]= size(UA);
N_PARTICLES=co; 

jj=1;
for i = 1 : N_PARTICLES
    % Store
    TT=(Table.NewID==UA(i));
    xx=Table.x(TT);
    yy=Table.y(TT); 
    zz=Table.z(TT);
    X=[xx yy zz];
    time=cumsum(Table.Fixed_time2(TT));

dis(jj:jj+length(xx)-2)=sqrt((xx(2:length(xx))-xx(1:length(xx)-1)).^2+...
    (yy(2:length(xx))-yy(1:length(xx)-1)).^2+...
    (zz(2:length(xx))-zz(1:length(xx)-1)).^2);
jj=jj+length(xx)-1;

end

%Scaled
yd=dis./mean(dis);
yd_Wk1=yd;

%Lfa1w

dt=0.1;
x=0.01:dt:100;

% A) %%%%%%%% Liver Data  
num=xlsread('All_datasets_for_curate6.xlsx','Lfa1_WTDis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Lfa1_KODis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Wk1_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Wk4_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','New_AllDisT');%'AllData_Wk1');

Table=array2table(num, 'VariableNames',{...
    'Mouse'	'x'	'y' 'z'	'Unit'	'Category'...
    'Collection'	'TimeSlice'	'TrackID'...
    'ID'	'Time' 	'Dis'	' Timediff'	'Speed'	'CelID'	'Type'	'Fixed_time1' 'Fixed_time3' 'Fixed_time2' 'addT' 'NewID'});
  Table.Fixed_time2=Table.Fixed_time2;  
 
UA=unique(Table.NewID);
[co xx]= size(UA);
N_PARTICLES=co; 

jj=1;
for i = 1 : N_PARTICLES
    % Store
    TT=(Table.NewID==UA(i));
    xx=Table.x(TT);
    yy=Table.y(TT); 
    zz=Table.z(TT);
    X=[xx yy zz];
    time=cumsum(Table.Fixed_time2(TT));

dis(jj:jj+length(xx)-2)=sqrt((xx(2:length(xx))-xx(1:length(xx)-1)).^2+...
    (yy(2:length(xx))-yy(1:length(xx)-1)).^2+...
    (zz(2:length(xx))-zz(1:length(xx)-1)).^2);
jj=jj+length(xx)-1;

end

%Scaled
yd=dis./mean(dis);
yd_lfa1w=yd;


%Wk4

dt=0.1;
x=0.01:dt:100;

% A) %%%%%%%% Liver Data  
%num=xlsread('All_datasets_for_curate6.xlsx','Lfa1_WTDis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Lfa1_KODis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Wk1_Dis');%'AllData_Wk1');
num=xlsread('All_datasets_for_curate6.xlsx','Wk4_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','New_AllDisT');%'AllData_Wk1');

Table=array2table(num, 'VariableNames',{...
    'Mouse'	'x'	'y' 'z'	'Unit'	'Category'...
    'Collection'	'TimeSlice'	'TrackID'...
    'ID'	'Time' 	'Dis'	' Timediff'	'Speed'	'CelID'	'Type'	'Fixed_time1' 'Fixed_time3' 'Fixed_time2' 'addT' 'NewID'});
  Table.Fixed_time2=Table.Fixed_time2;  
 


UA=unique(Table.NewID);
[co xx]= size(UA);
N_PARTICLES=co; 

jj=1;
for i = 1 : N_PARTICLES
    % Store
    TT=(Table.NewID==UA(i));
    xx=Table.x(TT);
    yy=Table.y(TT); 
    zz=Table.z(TT);
    X=[xx yy zz];
    time=cumsum(Table.Fixed_time2(TT));

dis(jj:jj+length(xx)-2)=sqrt((xx(2:length(xx))-xx(1:length(xx)-1)).^2+...
    (yy(2:length(xx))-yy(1:length(xx)-1)).^2+...
    (zz(2:length(xx))-zz(1:length(xx)-1)).^2);
jj=jj+length(xx)-1;

end

%Scaled
yd=dis./mean(dis);
yd_Wk4=yd;

%lfa1k

dt=0.1;
x=0.01:dt:100;

% A) %%%%%%%% Liver Data  
%num=xlsread('All_datasets_for_curate6.xlsx','Lfa1_WTDis');%'AllData_Wk1');
num=xlsread('All_datasets_for_curate6.xlsx','Lfa1_KODis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Wk1_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Wk4_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','New_AllDisT');%'AllData_Wk1');

Table=array2table(num, 'VariableNames',{...
    'Mouse'	'x'	'y' 'z'	'Unit'	'Category'...
    'Collection'	'TimeSlice'	'TrackID'...
    'ID'	'Time' 	'Dis'	' Timediff'	'Speed'	'CelID'	'Type'	'Fixed_time1' 'Fixed_time3' 'Fixed_time2' 'addT' 'NewID'});
  Table.Fixed_time2=Table.Fixed_time2;  
 
UA=unique(Table.NewID);
[co xx]= size(UA);
N_PARTICLES=co; 

jj=1;
for i = 1 : N_PARTICLES
    % Store
    TT=(Table.NewID==UA(i));
    xx=Table.x(TT);
    yy=Table.y(TT); 
    zz=Table.z(TT);
    X=[xx yy zz];
    time=cumsum(Table.Fixed_time2(TT));

dis(jj:jj+length(xx)-2)=sqrt((xx(2:length(xx))-xx(1:length(xx)-1)).^2+...
    (yy(2:length(xx))-yy(1:length(xx)-1)).^2+...
    (zz(2:length(xx))-zz(1:length(xx)-1)).^2);
jj=jj+length(xx)-1;

end

%Scaled
yd=dis./mean(dis);
yd_lfa1k=yd;

%1326

dt=0.1;
x=0.01:dt:100;

% A) %%%%%%%% Liver Data  
%num=xlsread('All_datasets_for_curate6.xlsx','Lfa1_WTDis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Lfa1_KODis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Wk1_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Wk4_Dis');%'AllData_Wk1');
num=xlsread('All_datasets_for_curate6.xlsx','New_AllDis');%'AllData_Wk1');

Table=array2table(num, 'VariableNames',{...
    'Mouse'	'x'	'y' 'z'	'Unit'	'Category'...
    'Collection'	'TimeSlice'	'TrackID'...
    'ID'	'Time' 	'Dis'	' Timediff'	'Speed'	'CelID'	'Type'	'Fixed_time1' 'Fixed_time3' 'Fixed_time2' 'addT' 'NewID'});
  Table.Fixed_time2=Table.Fixed_time2;  
 
UA=unique(Table.NewID);
[co xx]= size(UA);
N_PARTICLES=co; 

jj=1;
for i = 1 : N_PARTICLES
    % Store
    TT=(Table.NewID==UA(i));
    xx=Table.x(TT);
    yy=Table.y(TT); 
    zz=Table.z(TT);
    X=[xx yy zz];
    time=cumsum(Table.Fixed_time2(TT));

dis(jj:jj+length(xx)-2)=sqrt((xx(2:length(xx))-xx(1:length(xx)-1)).^2+...
    (yy(2:length(xx))-yy(1:length(xx)-1)).^2+...
    (zz(2:length(xx))-zz(1:length(xx)-1)).^2);
jj=jj+length(xx)-1;

end

%Scaled
yd=dis./mean(dis);
yd_1326=yd;



%[sims xxx]=size(yd);
%data=yd+0.001; 
%delta=Table.Fixed_time2(1);
%global data vel delta xmin1 xmin2 

clf
figure(1); 
[mu, xmin, L]=tailfit(yd_Wk1+0.1);
mu_tail_Wk1=mu
r_min_Wk1=xmin
p1=tailplot_Wk1(yd_Wk1+0.1, r_min_Wk1, mu_tail_Wk1);

hold on

[mu, xmin, L]=tailfit(yd_Wk4+0.1);
mu_tail_Wk4=mu
r_min_Wk4=xmin
p2=tailplot_Wk4(yd_Wk4+0.1, r_min_Wk4, mu_tail_Wk4);

hold on

[mu, xmin, L]=tailfit(yd_lfa1k+0.1);
mu_tail_lfa1k=mu
r_min_lfa1k=xmin
p3=tailplot_Lfa1KO(yd_lfa1k+0.1, r_min_lfa1k, mu_tail_lfa1k);

hold on

[mu, xmin, L]=tailfit(yd_lfa1w+0.1);
mu_tail_lfa1w=mu
r_min_lfa1w=xmin
p4=tailplot_LFa1WT(yd_lfa1w+0.1, r_min_lfa1w, mu_tail_lfa1w);

hold on

[mu, xmin, L]=tailfit(yd_1326+0.1);
mu_tail_1326=mu
r_min_1326=xmin
p5=tailplot_New1326(yd_1326+0.1, r_min_1326, mu_tail_1326);



xlabel('Scaled displacement (\rho)')
ylabel('Pr (P\geq\rho)')
ax = gca
ax.FontSize = 14;

legend([p1 p2 p3 p4 p5],{"Wk1"+ " (\mu="+sprintf('%.2f',mu_tail_Wk1)+")", ...
    "Lfa1k"+ " (\mu"+sprintf('%.2f',mu_tail_lfa1k)+")", ...
    "Wk4"+ " (\mu="+sprintf('%.2f',mu_tail_Wk4)+")", ...
    "Lfa1w"+ " (\mu"+sprintf('%.2f',mu_tail_lfa1w)+")", ...
    "1326"+ " (\mu"+sprintf('%.2f',mu_tail_1326)+")"}, ...
    'FontSize',14,'Location','northeast')

