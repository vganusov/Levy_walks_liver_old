clear all; 
dt=0.1;
x=0.01:dt:100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Notes 
%LfaWT; KO; w1k1; 1k2; AllneW1326: 27.5; 27.5;  6.51; 3.25; 26; secs
%P14; OT1; 1 min
%LN 20 sec;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data (open the relevant data links): 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A) %%%%%%%% Liver Data  
%num=xlsread('All_datasets_for_curate6.xlsx','Lfa1_WTDis');%'AllData_Wk1');
num=xlsread('All_datasets_for_curate6.xlsx','Lfa1_KODis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Wk1_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','Wk4_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','New_26DisT');%'AllData_Wk1');

Table=array2table(num, 'VariableNames',{...
    'Mouse'	'x'	'y' 'z'	'Unit'	'Category'...
    'Collection'	'TimeSlice'	'TrackID'...
    'ID'	'Time' 	'Dis'	' Timediff'	'Speed'	'CelID'	'Type'	'Fixed_time1' 'Fixed_time3' 'Fixed_time2' 'addT' 'NewID'});
  Table.Fixed_time2=Table.Fixed_time2;  
 

%B) %%%%%%%% Liver Data P14 and OT1  
%num=xlsread('All_datasets_for_curate6.xlsx','PyTCR_Dis');%'AllData_Wk1');
%num=xlsread('All_datasets_for_curate6.xlsx','OT1_Dis');%'AllData_Wk1');

%Table=array2table(num, 'VariableNames',{...
%    'ID'	'time' 'Fixed_time2' 'x' 'y' 'z'	'NewID' 'Dis'});
%    Table.Fixed_time2=Table.Fixed_time2.*60;


%C) %%%%%%%%%%%%% LN data
%num=xlsread('0409_Doc5_P14','Position(2)');%'AllData_Wk1');
%num=xlsread('0409_Doc6_P14','Position(2)');%'AllData_Wk1');
%num=xlsread('0409_Doc7_P14','Position(2)');%'AllData_Wk1');
%num=xlsread('0209_Doc6_P14','Position(2)');%'AllData_Wk1');
%num=xlsread('0209_Doc7_P14','Position(2)');%'AllData_Wk1');
%num=xlsread('0209_Doc8_P14','Position(2)');%'AllData_Wk1');

%Table=array2table(num, 'VariableNames',{...    
%'x' 'y' 'z'	'time' 'Fixed_time2' 'NewID2x'...
%'TrackID'	'NewID'});
%    Table.Fixed_time2=Table.Fixed_time2*20;


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

% open for non-Scaled 
%yd=dis;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sims xxx]=size(yd);
data=yd+0.001; 
%delta=Table.Fixed_time2(1);
global data vel delta xmin1 xmin2 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tail plot
figure(); clf
[mu, xmin, L]=tailfit(yd+0.1);
tailplot_step(yd+0.1, xmin, mu);
mu_tail=mu
r_min=xmin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(); clf
%subplot(1,4,1);;hold on
dt2=0.3;
dt=0.3;
h=histogram(data,'Normalization','probability','BinWidth',dt); %
options = optimset('Display','iter');
x=0.01:dt2:1000;
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%Stable Caushy
sig_022=1.2;
[x_3,fval_3]=fminsearch('objectivefn_Stable_Caushy',sig_022);
pd_c= makedist('Stable','alpha',1,'beta',0,'gam',x_3(1)^2,'delta',0);
y_est_3 = pdf(pd_c,x);
z_est_3=y_est_3;
plot(x,z_est_3.*dt2, '--','color', 'b',...
    'LineWidth', 4)
nLL_Caushy=fval_3
para_Caushy=x_3(1)^2
[aic_Cau,bicCau] = aicbic(-nLL_Caushy,1,length(yd))

hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%Stable Levy
sig_0=0.02;
[x_4L,fval_4L]=fminsearch('objectivefn_Stable_Levy',sig_0);
pd_lL= makedist('Stable','alpha',0.5,'beta',1,'gam',x_4L(1)^2,'delta',0);
y_est_4L = pdf(pd_lL,x);
z_est_4=y_est_4L;
plot(x,z_est_4.*dt2, ':','color',  'r',...
    'LineWidth', 4)
nLL_Levy=fval_4L
para_Levy=x_4L^2
[aic_Lv,bicLv] = aicbic(-nLL_Levy,1,length(yd))


%%%%%%%%%%%%%%%%%%%%%%%%%%Gen_Perato_FULL 5 2 parameter, location=0;
pd_gp=fitdist(data','GeneralizedPareto');
fval_gp=negloglik(pd_gp);
y_est_8gp = gppdf(x,pd_gp.k,pd_gp.sigma,pd_gp.theta);
z_est_8gp=y_est_8gp;
plot(x,z_est_8gp.*dt2, '-','color', 'r',......
    'LineWidth', 4, 'MarkerSize', 4)
nLL_GP=fval_gp
para_GP=[pd_gp.k,pd_gp.sigma]
[aic_GP,bicGP] = aicbic(-nLL_GP,2,length(yd))

%%%%%%%%%%%%%%%%%%%%%%%%%%Exp 
pd_ep=fitdist(data','Exp');
fval_6ep=negloglik(pd_ep);
y_est_6ep = pdf(pd_ep,x);
z_est_6ep=y_est_6ep./sum(y_est_6ep);
nLL_Exp=fval_6ep
para_Exp=pd_ep.mu
[aic_Exp,bicExp] = aicbic(-nLL_Exp,1,length(yd))


%%%%%%%%%%%%%%%%%%%%%%%%%%Half Normal
pd_nhn=fitdist(data','HalfNormal');
fval_HfN=negloglik(pd_nhn);
y_est_5hn = pdf(pd_nhn,x);;
z_est_5hn=y_est_5hn;
plot(x,z_est_5hn.*dt2, 'square','color',  'g',...
'LineWidth', 2)
nLL_HfN=fval_HfN
para_HfN=pd_nhn.sigma
[aic_HfN,bicHfN] = aicbic(-nLL_HfN,1,length(yd))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot
legend('Data','Caushy 1p','Levy 1p',...
   'Generalized Perato 2p',...
   'Half Normal 1p')
xlabel('Step lengths(\mum)');ylabel('Prb')

xlim([0.1 100]);ylim([0.00001 0.2]); %non scaled

grid on;
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

% printing parameters
% alternative models fitted to displacement data
[para_Levy para_Caushy para_GP para_Exp para_HfN]

% tailfit analysis
[r_min mu_tail]  % mu and rmin
[r_min mu_tail-1];  % alpha=mu-1 and rmin

% AIC
aic=[aic_Lv aic_Cau aic_GP aic_Exp aic_HfN];
uint64(aic)
% weights
w=exp((min(aic)-aic)/2)/sum(exp((min(aic)-aic)/2))


