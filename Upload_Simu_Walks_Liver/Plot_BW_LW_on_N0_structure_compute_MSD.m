% Plot BW & LW on N0_structure & compute MSD

clear all; hold on;  
rng shuffle
warning off;
%clf;

figure()
%alpha=4.6 

r0=0.1;
nstpes=50; %no. of steps
ncells=200;%200; %no. of cells
nrep=1; %no. of respaples


mean_R=2.1;
alphR=1.5;%alpha=1.5; 3.5;  
k_gpR=1/alphR;
%xminR=0.1;
xminR=mean_R*((alphR-1)/alphR);
sig_gpR=xminR/alphR; 
mean_PaR=alphR*xminR/(alphR-1);
sig_GR=mean_PaR*sqrt(pi()/2);
mean_GuR=sig_GR*sqrt(2/pi());
max_stepR=300;

for j=1:ncells
    j

for i=1:nrep;
    clear poi
    poi(1,1)=rand()*512; poi(2,1)=rand()*512; poi(3,1)=rand()*44; %initial points
    

    brlen1R=gprnd(k_gpR,sig_gpR,sig_gpR/k_gpR,1,nstpes);
    brlenR=min(brlen1R,max_stepR); 
    rn = randn(3,nstpes);
    poi(:,2:nstpes+1)=[bsxfun(@rdivide,rn,sqrt(sum(rn.^2))).*brlenR];

    
    x=cumsum(poi')';
    plot3(x(1,:),x(2,:),x(3,:));
    hold on

    sq_dis(i,:)=(x(1,:).^2+x(2,:).^2+x(3,:).^2);
end    
    msd=log(mean(sq_dis+0.001));
    time=log([1:nstpes]);
end

grid on
view(150,75)

axis([1 512 1 512 1 44]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Levy

clear all; hold on;  
rng shuffle
warning off;
%clf;

figure()


r0=0.1;
nstpes=50; %no. of steps
ncells=200;%200; %no. of cells
nrep=1; %no. of respaples

mean_R=2.1;
alphR=1.5;%alpha=1.5; 3.5;  
k_gpR=1/alphR;
xminR=mean_R*((alphR-1)/alphR);
sig_gpR=xminR/alphR; 
mean_PaR=alphR*xminR/(alphR-1);
sig_GR=mean_PaR*sqrt(pi()/2);
mean_GuR=sig_GR*sqrt(2/pi());
max_stepR=300;

for j=1:ncells
    j

for i=1:nrep;
    clear poi
    poi(1,1)=rand()*512; poi(2,1)=rand()*512; poi(3,1)=rand()*44; %initial points
    

    brlen1R=gprnd(k_gpR,sig_gpR,sig_gpR/k_gpR,1,nstpes);
    brlenR=min(brlen1R,max_stepR); 
    rn = randn(3,nstpes);
    poi(:,2:nstpes+1)=[bsxfun(@rdivide,rn,sqrt(sum(rn.^2))).*brlenR];

    
    x=cumsum(poi')';
    plot3(x(1,:),x(2,:),x(3,:));
    hold on

    sq_dis(i,:)=(x(1,:).^2+x(2,:).^2+x(3,:).^2);
end    
    msd=log(mean(sq_dis+0.001));
    time=log([1:nstpes]);
end

grid on
view(150,75)

axis([1 512 1 512 1 44]);






