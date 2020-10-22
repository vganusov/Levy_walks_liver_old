clear all; 
dt=0.1;
x=0.01:dt:100;


%OT1 data

%B) %%%%%%%% Liver Data P14 and OT1  
%num=xlsread('All_datasets_for_curate6.xlsx','PyTCR_Dis');%'AllData_Wk1');
num=xlsread('All_datasets_for_curate6.xlsx','OT1_Dis');%'AllData_Wk1');

Table=array2table(num, 'VariableNames',{...
    'ID'	'time' 'Fixed_time2' 'x' 'y' 'z'	'NewID' 'Dis'});
    Table.Fixed_time2=Table.Fixed_time2.*60;

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


% open the following for scaled vs non-scaled;

% open for Scaled 
yd=dis./mean(dis);
yd_OT1=yd;

%P14 data

dt=0.1;
x=0.01:dt:100;
%B) %%%%%%%% Liver Data P14 and OT1  
num=xlsread('All_datasets_for_curate6.xlsx','PyTCR_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','OT1_Dis');%'AllData_Wk1');

Table=array2table(num, 'VariableNames',{...
    'ID'	'time' 'Fixed_time2' 'x' 'y' 'z'	'NewID' 'Dis'});
    Table.Fixed_time2=Table.Fixed_time2.*60;


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


% open the following for scaled vs non-scaled;

% open for Scaled 
yd=dis./mean(dis);
yd_P14=yd;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tail plot
clf
figure(1); 
[mu, xmin, L]=tailfit(yd_OT1+0.1);
mu_tail_OT1=mu
r_min_OT1=xmin
p1=tailplot_step_OT1(yd_OT1+0.1, r_min_OT1, mu_tail_OT1);

hold on

[mu, xmin, L]=tailfit(yd_P14+0.1);
mu_tail_P14=mu
r_min_P14=xmin
p2=tailplot_step_P14(yd_P14+0.1, r_min_P14, mu_tail_P14);
xlabel('Scaled displacement (\rho)')
ylabel('Pr (P\geq\rho)')
ax = gca
ax.FontSize = 14;


% legend({"OT1"+ " (\mu="+sprintf('%.2f',mu_tail_OT1)+")", ...
%     "P14"+ " (\mu="+sprintf('%.2f',mu_tail_P14)+")"}, ... 
%     'FontSize',14,'Location','northwest')
% 

legend([p1 p2],{"OT1"+ " (\mu="+sprintf('%.2f',mu_tail_OT1)+")", ...
    "P14"+ " (\mu="+sprintf('%.2f',mu_tail_P14)+")"}, ... 
    'FontSize',14,'Location','northeast')


