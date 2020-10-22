% Computes branch length distribution -Liver 
% Fits Pareto tail to branch length data,
% Computes and plot histogram of Angle distribution
% Computes and plot histogram of stratightness distribution

clear all; close all; clf 
load('IMG_OPENED_3D_Evans_old.mat');
skel = Skeleton3D(IMG_OPENED_3D_Evans_old); 
skel_2=skel;
clear skel
skel=skel_2(1:512,1:512,:);
w = size(skel,1);
l = size(skel,2);
h = size(skel,3);
[~,node,link] = Skel2Graph3D(skel,0);
wl = sum(cellfun('length',{node.links}));
skel2 = Graph2Skel3D(node,link,w,l,h);
[~,node2,link2] = Skel2Graph3D(skel2,0);
wl_new = sum(cellfun('length',{node2.links}));
while(wl_new~=wl)
    wl = wl_new;   
     skel2 = Graph2Skel3D(node2,link2,w,l,h);
     [A2,node2,link2] = Skel2Graph3D(skel2,0);
     wl_new = sum(cellfun('length',{node2.links}));
end;
hold on;
coun=1;
for i=1:length(node2)
    x1 = node2(i).comx;
    y1 = node2(i).comy;
    z1 = node2(i).comz;
    if(node2(i).ep==1)
        ncol = 'c';
    else
        ncol = 'y';
    end;
    for j=1:length(node2(i).links)   
        if(node2(node2(i).conn(j)).ep==1)
            col='k';
        else
            col='r'; 
        end;
        if(node2(i).ep==1)
            col='k';
        end;
        for k=1:length(link2(node2(i).links(j)).point)-1            
           [x3,y3,z3]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k));
            [x2,y2,z2]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k+1));
            %%%%%%%%%%%%%%%%%%%%%%%%%
            sto_x2(i,j,k)=x2;
            sto_y2(i,j,k)=y2;
            sto_z2(i,j,k)=z2;
            sto_x3(i,j,k)=x3;
            sto_y3(i,j,k)=y3;
            sto_z3(i,j,k)=z3;
            
        end;
    end;
   
end;
%%%%%
[sz_n sz_l sz_p] =size(sto_x2); 
mat_st=[];
coun=1;coun2=1;br=1; counta=1;
for i=1:length(node2)
    if node2(i).ep==0
    for j=1:length(node2(i).links)
        clear tr
        tr=ismember(mat_st, node2(i).links(j));
        if sum(tr)==0
        mat_st = [mat_st node2(i).links(j)];
        bran_len=0;
        [xf3,yf3,zf3]=ind2sub([w,l,h],link2(node2(i).links(j)).point(1));
        [xl2,yl2,zl2]=ind2sub([w,l,h],link2(node2(i).links(j)).point(length(link2(node2(i).links(j)).point)));
        bran_len_str(br)=bran_len+sqrt((xf3-xl2)^2+(yf3-yl2)^2+(zf3-zl2)^2);
        br=br+1;
            for k=1:length(link2(node2(i).links(j)).point)-1            
            [xe3,ye3,ze3]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k));
            [xe2,ye2,ze2]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k+1));
            step_l(coun)=sqrt((xe3-xe2)^2+(ye3-ye2)^2+(ze3-ze2)^2);
            bran_len=bran_len+sqrt((xe3-xe2)^2+(ye3-ye2)^2+(ze3-ze2)^2);
            coun=coun+1;
            end
            branch_l(coun2)=bran_len;
            coun2=coun2+1;
            for k=1:length(link2(node2(i).links(j)).point)-2            
            [xa3,ya3,za3]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k));
            [xa2,ya2,za2]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k+1));
            [xa1,ya1,za1]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k+2)); %farthest point
            vec1    = [xa3,ya3,za3];
            vec2     = [xa2,ya2,za2];
           vec3    = [xa1,ya1,za1];
            dv1     = vec1 - vec2;
            dv2     = vec3 - vec2;
angle(j) = radtodeg(acos(max(-1,min(1,dot(dv1,dv2)/(norm(dv1)*norm(dv2))))));
            step_ang(counta)=180-radtodeg(acos(max(-1,min(1,dot(dv1,dv2)/(norm(dv1)*norm(dv2))))));
            counta=counta+1;
            end
        end
    end
    end
end

cu=1;cu2=1;
for i=1:length(node2)
    if node2(i).ep==0
        num_nodes(cu)=length(node2(i).links); 
        cu=cu+1;
    for j=1:length(node2(i).links)-1
         if length(link2(node2(i).links(j+1)).point)>1 & length(link2(node2(i).links(j)).point)>1
        [xb3,yb3,zb3]=ind2sub([w,l,h],link2(node2(i).links(j)).point(2));
        [xb2,yb2,zb2]=ind2sub([w,l,h],link2(node2(i).links(j)).point(1));
        [xb1,yb1,zb1]=ind2sub([w,l,h],link2(node2(i).links(j+1)).point(2));
            vec1    = [xb3,yb3,zb3];
            vec2     = [xb2,yb2,zb2];
            vec3    = [xb1,yb1,zb1];
            dv1     = vec1 - vec2;
            dv2     = vec3 - vec2;
            br_ang2(cu2)=radtodeg(acos(max(-1,min(1,dot(dv1,dv2)/(norm(dv1)*norm(dv2))))));
            cu2=cu2+1;
         end
    end
    end
end

cou=1;
for i=1:length(branch_l)
    if branch_l(i)>0 
branch_lf(cou)=branch_l(i);
bran_len_str_f(cou)=bran_len_str(i);
cou=cou+1;
    end
end
dt2=1;

cou=1;
for i=1:length(br_ang2)   
    if br_ang2(i)~=90 
br_ang2_f(cou)=br_ang2(i);
cou=cou+1;
    end
end


% plotting figure (4 panels)

dt=0.1;
x=0.01:dt:200;
yd=branch_lf;
[sims xxx]=size(yd);
data=yd+0.0001; 
global data  
figure(1)

%panel A - setting up plot
subplot(1,4,1); hold on
xlabel('Branch length d (\mum)'); ylabel('Prb (d)')
dt2=1;
h=histogram(data,'Normalization','probability','BinWidth',dt2); %
h.FaceColor=[0.5 0.5 0.5]
h.EdgeColor='none'
options = optimset('Display','iter');

%%%%%%%%%%%%%%%%%%%%%%%%%%Stable Caushy
sig_022=1.2;
[x_3,fval_3]=fminsearch('objectivefn_Stable_Caushy',sig_022);
pd_c= makedist('Stable','alpha',1,'beta',0,'gam',x_3(1)^2,'delta',0);
y_est_3 = pdf(pd_c,x);
z_est_3=y_est_3;
plot(x,z_est_3, '--','color', 'b',...
    'LineWidth', 4)
nLL_Caushy=fval_3
para_Caushy=x_3(1)^2
[aic_Cau,bicCau] = aicbic(-nLL_Caushy,1,length(yd))

%%%%%%%%%%%%%%%%%%%%%%%%%%Stable Levy
sig_0=0.02;
[x_4L,fval_4L]=fminsearch('objectivefn_Stable_Levy',sig_0);
pd_lL= makedist('Stable','alpha',0.5,'beta',1,'gam',x_4L(1)^2,'delta',0);
y_est_4L = pdf(pd_lL,x);
z_est_4=y_est_4L;
plot(x,z_est_4, ':','color',  'r',...
    'LineWidth', 4)
nLL_Levy=fval_4L
para_Levy=x_4L^2
[aic_Lv,bicLv] = aicbic(-nLL_Levy,1,length(yd))


%%%%%%%%%%%%%%%%%%%%%%%%%%Gen_Perato_FULL 5 2 parameter, location=0;
pd_gp=fitdist(data','GeneralizedPareto');
fval_gp=negloglik(pd_gp);
y_est_8gp = gppdf(x,pd_gp.k,pd_gp.sigma,pd_gp.theta);
z_est_8gp=y_est_8gp;
plot(x,z_est_8gp, '-','color', 'r',......
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
plot(x,z_est_5hn, 'square','color',  'g',...
'LineWidth', 2)
nLL_HfN=fval_HfN
para_HfN=pd_nhn.sigma
[aic_HfN,bicHfN] = aicbic(-nLL_HfN,1,length(yd))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot
lgd=legend('Data','Cauchy','Levy',...
   'Generalized Pareto',...
   'Half Normal')
lgd.Location='southwest'
xlabel('Branch length d (\mum)'); ylabel('Prb(d)')

xlim([1.0 100]);ylim([0.00001 0.2]); %non scaled

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
aic=[aic_Lv aic_Cau aic_GP aic_Exp aic_HfN]
uint64(aic)
% weights
w=exp((min(aic)-aic)/2)/sum(exp((min(aic)-aic)/2))


% panel B: tail fit
[mu, xmin, L]=tailfit(yd)
subplot(1,4,2)
tailplot(yd, xmin, mu);

mu_tail=mu
r_min=xmin


% panel C - branching angle distribution
subplot(1,4,3)
dt5=10;
angles1=real(br_ang2_f)
% mean angle for acute angles
mean(angles1(angles1<90))

h1=histogram(real(br_ang2_f),'Normalization','probability','BinWidth',dt5); 
xlabel('Angle between brances at nodes (\theta)');ylabel('Prb(\theta)')
xticks([0 45 90 135 180])
xline(90,'--')
% fraction > 90
sum(real(br_ang2_f)>90)/length(real(br_ang2_f))


%panel D - straightness index
subplot(1,4,4)
dt4=0.05;
histogram(bran_len_str_f./branch_lf, 'Normalization','probability','BinWidth',dt4);
xlabel('Branch straightness index');ylabel('Prb(straightness)')
% mean SI
mean(bran_len_str_f./branch_lf)





