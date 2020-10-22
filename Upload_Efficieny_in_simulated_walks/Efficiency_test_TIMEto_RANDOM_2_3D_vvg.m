% Creates efficiency -time- graphs for no-structure, 
% both Brownian (mu=4.5) and Levy steps (mu=2.5) 

clear all;

%mu=20;
%sigma=2;
%number of branches encountering in a walk 
%sim=1000; % number of branches & cells
nstep=100; %number of steps

km=1000;%1000; %numbr of cells;

st_x=nan(km,nstep+1);
st_y=nan(km,nstep+1);
st_z=nan(km,nstep+1);
st_time=nan(km,nstep+1);
%data_an=angle_dat;


sense=40; %sensing distance
mod_radius=100; %modeling shphere radius  # one set of simulations
max_step=100; %max step size
celn=km;% number of cells (k,j)
targ_cells=1; % number of target cells (t)
a=50; b=mod_radius; %between 20 and 100 um

for t=1:targ_cells  
%%%% Target placement
len_rnd_point =90;% (b-a).*rand()+a;
  rn = randn(3,1); % Use a large n
  poi= bsxfun(@rdivide,rn,sqrt(sum(rn.^2)))*len_rnd_point;
  xtg=poi(1);
  ytg=poi(2);
  ztg=poi(3);
  
%      
  
alph=1.5;
k_gp=1/alph;
xmin=2.1;
sig_gp=xmin/alph; 
mean_Pa=alph*xmin/(alph-1)
sig_G=mean_Pa*sqrt(pi()/2);
mean_Gu=sig_G*sqrt(2/pi())
%brlen1=random('HalfNormal',0,sig_G, 1,1); 
%brlen=min(brlen1,max_step);
%brlen1=gprnd(k_gp,sig_gp,sig_gp/k_gp,1,1);
%brlen=min(brlen1,max_step);

mean_RG=mean_Gu;
alphRG=4.5;%alpha=4.5;  1
k_gpRG=1/alphRG;
xminRG=mean_RG*((alphRG-1)/alphRG);
sig_gpRG=xminRG/alphRG; 
mean_PaRG=alphRG*xminRG/(alphRG-1);
sig_GRG=mean_PaRG*sqrt(pi()/2);
mean_GuRG=sig_GRG*sqrt(2/pi());
max_stepR=max_step;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BW

for k=1:celn

    clear st_brlen st_x st_y st_z x y z
        x(1)=0; y(1)=0; z(1)=0; %initial points

    t,k
flag_targ=0;
l=0;
while flag_targ==0
 % for l=1:kmax; %cells
    l=l+1;

%flag=0;
%while flag==0 

%flag=0;
%while flag==0 
    
%brlen1=gprnd(k_gp,sig_gp,sig_gp/k_gp,1,1);
%brlen1=random('HalfNormal',0,sig_G, 1,1); 
%brlen=min(brlen1,max_step);


brlen1RG=gprnd(k_gpRG,sig_gpRG,sig_gpRG/k_gpRG,1,1);
brlen=min(brlen1RG,max_step);

%brlen=abs(normrnd(0,1));
%brlen=exprnd(2);
% brlen=gprnd(0.52,18.2,0);
%%%%%%%%%%%%%%Levy
 %  pd3 = makedist('Stable','alpha',0.5,'beta',1,'gam',0.15,'delta',0);
 %  step=icdf(pd3,rand());
% brlen=gprnd(1.10,1.2,0);
 %   if step<100 & step>0.1
 %       fll=0;
 %   end
  % brlen=abs(step);
    %%%%%%%%%%%%%%%%%%%%%%
%brlen=abs(levy(1,1,0.5)); %(no.steps,no. of dimensions, beta-powerlaw exponent)   
%if brlen>0.1 %minimum branch length
%    flag=1;
%end
%if brlen<max_step %maxmum branch length
%    flag=1;
%end
%end
st_brlen(l)=brlen ;


flag_area=0;  
while flag_area==0;  
%branch end-point cordinates
%  x(i+1)=poi(1);y(i+1)=poi(2);z(i+1)=poi(3);  
  rn = randn(3,1); % Use a large n
  poi= bsxfun(@rdivide,rn,sqrt(sum(rn.^2)))*brlen;
%if 
%end
%branch end-point cordinate
  x(l+1)=x(l)+poi(1);y(l+1)=y(l)+poi(2);z(l+1)=z(l)+poi(3);  
  position_radius=sqrt(x(l+1)^2 +y(l+1)^2 +z(l+1)^2);
  if position_radius<mod_radius
  flag_area=1;
  end
end
 

st_x(l,1:length(x))=x;
st_y(l,1:length(x))=y;
st_z(l,1:length(x))=z;
st_time(l,1:length(x))=1;

%length to target
len_tar=sqrt((x(l+1)-xtg).^2+(y(l+1)-ytg).^2+(z(l+1)-ztg).^2);
if len_tar<=sense
    flag_targ=1;
end

end
%end

st_time_to_target_Gau(k)=l;
dista_travel_Gau(k)=sum(st_brlen);
end

%figure(); 
%subplot(3,1,1);
%plot3(x,y,z,'o')
%hold on
%plot3(x,y,z)
%hold on 
%plot3(xtg,ytg,ztg,'o')



%%%%%%%%%%%%%%%%%%%%%%%%%% LW
for j=1:celn
    t,j
        clear st_brlen st_x st_y st_z x y z
    x(1)=0; y(1)=0; z(1)=0; %initial points

flag_targ=0;
l=0;
while flag_targ==0
 % for l=1:kmax; %cells
    l=l+1;

%flag=0;
%while flag==0 

%flag=0;
%while flag==0 
    
brlen1=gprnd(k_gp,sig_gp,sig_gp/k_gp,1,1);
%brlen1=random('HalfNormal',0,sig_G, 1,1); 
brlen=min(brlen1,max_step);


%brlen=abs(normrnd(0,1));
%brlen=exprnd(2);
% brlen=gprnd(0.52,18.2,0);
%%%%%%%%%%%%%%Levy
 %  pd3 = makedist('Stable','alpha',0.5,'beta',1,'gam',0.15,'delta',0);
 %  step=icdf(pd3,rand());
% brlen=gprnd(1.10,1.2,0);
 %   if step<100 & step>0.1
 %       fll=0;
 %   end
  % brlen=abs(step);
    %%%%%%%%%%%%%%%%%%%%%%
%brlen=abs(levy(1,1,0.5)); %(no.steps,no. of dimensions, beta-powerlaw exponent)   
%if brlen>0.1 %minimum branch length
%    flag=1;
%end
%if brlen<max_step %maxmum branch length
%    flag=1;
%end
%end
st_brlen(l)=brlen ;


flag_area=0;  
while flag_area==0;  
%branch end-point cordinates
%  x(i+1)=poi(1);y(i+1)=poi(2);z(i+1)=poi(3);  
  rn = randn(3,1); % Use a large n
  poi= bsxfun(@rdivide,rn,sqrt(sum(rn.^2)))*brlen;
%if 
%end
%branch end-point cordinate
  x(l+1)=x(l)+poi(1);y(l+1)=y(l)+poi(2);z(l+1)=z(l)+poi(3);  
  position_radius=sqrt(x(l+1)^2 +y(l+1)^2 +z(l+1)^2);
  if position_radius<mod_radius
  flag_area=1;
  end
end
 

st_x(l,1:length(x))=x;
st_y(l,1:length(x))=y;
st_z(l,1:length(x))=z;
st_time(l,1: length(x))=1;

%length to target
len_tar=sqrt((x(l+1)-xtg).^2+(y(l+1)-ytg).^2+(z(l+1)-ztg).^2);
if len_tar<=sense
    flag_targ=1;
end

end
%end

st_time_to_target_Par(j)=l;
dista_travel_Par(j)=sum(st_brlen);
end

st_time_Gau(t,:)=st_time_to_target_Gau;
st_dist_Gau(t,:)=dista_travel_Gau;
st_time_Par(t,:)=st_time_to_target_Par;
st_dist_Par(t,:)=dista_travel_Par;
end
%figure(); 
%subplot(3,1,1);
%plot3(x,y,z,'o')
%hold on
%plot3(x,y,z)
%hold on 
%plot3(xtg,ytg,ztg,'o')


dt=5;
hold on
figure(1);hold on
h=histogram(st_time_to_target_Par,'Normalization','probability','BinWidth',dt);%,...
%    'EdgeColor',[0.4660 0.6740 0.1880], 'FaceColor', [0.4940 0.1840 0.5560], 'LineWidth', 2);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
%xlim([0.1 400]);ylim([0.00001 1])

h=histogram(st_time_to_target_Gau,'Normalization','probability','BinWidth',dt);%,...
%    'EdgeColor',[0.4660 0.6740 0.1880], 'FaceColor', [0.4940 0.1840 0.5560], 'LineWidth', 2);

legend('LW (\mu=2.5)','BW (\mu=4)')

%linear diances of branch points from the initial point
%brlen_cum= cumsum(st_brlen);

meantime_LW=mean(st_time_to_target_Par)
mediantime_LW=median(st_time_to_target_Par)
sigmatime_LW=sqrt(var(st_time_to_target_Par))

meantime_BW=mean(st_time_to_target_Gau)
mediantime_BW=median(st_time_to_target_Gau)
sigmatime_BW=sqrt(var(st_time_to_target_Gau))

%figure
%pdSix = fitdist(st_time_to_target_Gau','Kernel','BandWidth',4);
%pdSix = fitdist(st_time_to_target_Gau','Gamma');
%pdSix = fitdist(st_time_to_target_Gau','GeneralizedExtremeValue');
%x = 0:10:100000;
%ySix = pdf(pdSix,x);
%plot(x,ySix*20,'k-','LineWidth',2)

%save Simulations_data_3D-LW-25-BW-45
save Simulations_data_3D-LW-25-BW-45



% loading simulation results (n=1000 cells)
clear all;
load Simulations_data_3D-LW-25-BW-45
%st_time_to_target_Par
%st_time_to_target_Gau


% exporting data
exportdata = transpose([st_time_to_target_Gau; st_time_to_target_Par])
csvwrite('simulation_data-3D-BW-45-LW-25-1000.csv',exportdata)

% plotting the figure
dt=0.2;
figure(1); 
hold on
h=histogram(log10(st_time_to_target_Gau),'Normalization','probability','BinWidth',dt);
%h.EdgeColor=[0 0 1]
h.EdgeColor='k'
h.FaceColor=[1 1 1]
h.LineWidth=2
h.LineStyle='-.'
%    'EdgeColor',[0.4660 0.6740 0.1880], 'FaceColor', [0.4940 0.1840 0.5560], 'LineWidth', 2);
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'log')
xlim([0 4]);ylim([0.001 1])

h=histogram(log10(st_time_to_target_Par),'Normalization','probability','BinWidth',dt);%,...
%    'EdgeColor',[0.4660 0.6740 0.1880], 'FaceColor', [0.4940 0.1840 0.5560], 'LineWidth', 2);
h.EdgeColor=[1 0 0]
h.FaceColor=[1 1 1]
h.EdgeAlpha=0.5
h.LineWidth=2
h.LineStyle='-'
set(gca,'FontSize',14)
legend('BW (\mu=4.5): T='+compose("%3.0f",round(meantime_BW)),'LW (\mu=2.5) T='+compose("%3.0f",round(meantime_LW)))
xlabel('Time to reach a target at a fixed distance (log_{10} steps)');
ylabel('Probability')
hold off

% calculating time to target for different number of T cells searching
% for the target

simu=1000
ncell=100;
min(randsample(st_time_to_target_Par,ncell,true)) % example

LW_T0d_cells=zeros(simu,1);
BW_T0d_cells=zeros(simu,1);

rng(1)   % setting seed for reproducing results later
for i = 1:simu
    LW_T0d_cells(i,1)= min(randsample(st_time_to_target_Par,ncell,true));
    BW_T0d_cells(i,1)= min(randsample(st_time_to_target_Gau,ncell,true));
end

mean(st_time_to_target_Par)  % search is done by 1 cell with LW
mean(st_time_to_target_Gau)  % search is done by 1 cell with BW
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
        LW_T0d_cells(i,1)= min(randsample(st_time_to_target_Par,ncell,true));
        BW_T0d_cells(i,1)= min(randsample(st_time_to_target_Gau,ncell,true));
    end
    LW_meanT(n,1) = mean(LW_T0d_cells);  % search is done by ncell cells with LW
    LW_meanT(n,2:3) = quantile(LW_T0d_cells,[0.025 0.975]);  % 95%CIs
    BW_meanT(n) = mean(BW_T0d_cells);  % search is done by ncell cells with BW
    BW_meanT(n,2:3) = quantile(BW_T0d_cells,[0.025 0.975]);  % 95%CIs
end

BW_meanT2=log10(BW_meanT);
LW_meanT2=log10(LW_meanT);



% plot time to target vs. # T cells - log scale

hold on
figure()
hold on
%errorbar(x,y,yneg,ypos,xneg,xpos,'o')
%plot(ncells,BW_meanT(:,1), '-ok','LineWidth',1.5);
%plot(ncells,LW_meanT(:,1), '--^r','LineWidth',1.5);
errorbar(ncells,BW_meanT2(:,1), BW_meanT2(:,1)-BW_meanT2(:,2), BW_meanT2(:,3)-BW_meanT2(:,1),'-ok','LineWidth',1.5);
errorbar(ncells,LW_meanT2(:,1), LW_meanT2(:,1)-LW_meanT2(:,2), LW_meanT2(:,3)-LW_meanT2(:,1),'--^r','LineWidth',1.5);

set(gca,'XScale','log')
set(gca,'YScale','linear')
set(gca,'FontSize',18)
pos = get(gca, 'Position');
pos(2) = 0.05;
pos(3) = 0.05;
%set(gca, 'Position', pos)
axis([0.9,1000,0,4])
legend('BW (\mu=4.5)', 'LW (\mu=2.5)','LineWidth',1.5);
% Label axes
xlabel('Number of cells searching for a single target')
ylabel('Mean time to target in 3D (log_{10} steps)')

% end of plot



