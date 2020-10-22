% Plot CBW and LW on structure & compute MSD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CBM
clear all; hold on;clf

gp=200; %cells
sims=50; %secs


%load('skel_LN.mat')
load('IMG_OPENED_3D_LN.mat');
skel_LN = Skeleton3D(IMG_OPENED_3D_LN); 

skel=skel_LN;

mean_R=2.1;
alphR=4;; 
k_gpR=1/alphR;
xminR=mean_R*((alphR-1)/alphR);
sig_gpR=xminR/alphR; 
mean_PaR=alphR*xminR/(alphR-1);
sig_GR=mean_PaR*sqrt(pi()/2);
mean_GuR=sig_GR*sqrt(2/pi());
max_stepR=300;

[max_x,max_y,max_z]=size(skel);

vec=zeros(max_x*max_y,2);
for i=1:max_z
[row,col]=find(skel(:,:,i));
sz=size(row);
sz_st(i)=sz(1);
vec(1:sz(1),:,i)=[row,col];
end 
  
c_count=1;
group=1;        
for g=1:gp
    g
clear st_walk_x st_walk_y st_walk_z st_walk_x_f st_walk_y_f st_walk_z_f 
init_point=[0 0];
while init_point==[0 0]
dr=round(1+ (sz_st(1)-1).*rand());
dr_z=round(1+ (44-1).*rand());
init_point=vec(dr,:,dr_z);
x=init_point(1);xi=x;
y=init_point(2);yi=y;
z=dr_z;zi=z;
x2=x;y2=y;z2=z;
end
st_walk_x(1)=x; st_walk_y(1)=y; st_walk_z(1)=z;
st_walk_x_f(1)=x;st_walk_y_f(1)=y; st_walk_z_f(1)=z; 

st_walk_xg((group-1)*sims+1)=x; st_walk_yg((group-1)*sims+1)=y; st_walk_zg((group-1)*sims+1)=z;
st_walk_x_fg((group-1)*sims+1)=x;st_walk_y_fg((group-1)*sims+1)=y; st_walk_z_fg((group-1)*sims+1)=z; 

count=1;

for c=1:sims

brlen1R=gprnd(k_gpR,sig_gpR,sig_gpR/k_gpR,1);
step=min(brlen1R,max_stepR);

    
for m=1:step 
    clear xn xp yn yp zn zp
if x-1<1; xn=1;else; xn=x-1;end
if x+1>max_x; xp=x; else; xp=x+1; end
if y-1<1; yn=1; else; yn=y-1; end
if y+1>max_y; yp=y; else; yp=y+1; end
if z-1<1; zn=1; else; zn=z-1; end
if z+1>max_z; zp=z; else; zp=z+1; end

%st_i(1)=nan; st_j(1)=nan; st_k(1)=nan;
cou=1; clear st_i st_j st_k
    for k=zn:zp
        for j=yn:yp
            for i=xn:xp
                for p=1:max_x*max_y
                    if i==vec(p,1,k)  
                    if j==vec(p,2,k) 
                        if i==x & j==y & k==z
                        else
                    st_i(cou)=i; st_j(cou)=j; st_k(cou)=k;
                    cou=cou+1;
                        end
                    else
                    end
                    end
                end
            end
        end
    end
    st_fl(c)=cou;
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    flag=0;
    iter=0;
    while flag==0
    dr2=round(1+ (cou-1-1).*rand()); 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if cou>1
    x2=st_i(dr2); y2=st_j(dr2); z2=st_k(dr2);
    else
    x2=x2; y2=y2; z2=z2;
    end  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    x2=st_i(dr2); y2=st_j(dr2); z2=st_k(dr2);
        if st_walk_x_f(count)==x2 & st_walk_y_f(count)==y2 & st_walk_z_f(count)==z2
                 iter=iter+1;
        else if count>2 & ...
                x2== st_walk_x_f(count-1) & y2== st_walk_y_f(count-1) &  z2==st_walk_z_f(count-1)
                     iter=iter+1;
            else if count>3 & ...
                    x2== st_walk_x_f(count-2) & y2== st_walk_y_f(count-2) &  z2==st_walk_z_f(count-2)
                     iter=iter+1;                  
                    else 
                        flag=1;
                  end
             end
        end
              if iter>2
                  flag=1;
              end
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    x=x2;y=y2;z=z2;
    st_walk_x_f(count+1)=x2;st_walk_y_f(count+1)=y2; st_walk_z_f(count+1)=z2;
    st_walk_x_fg((group-1)*sims+count+1)=x2;st_walk_y_fg((group-1)*sims+count+1)=y2; st_walk_z_fg((group-1)*sims+count+1)=z2; 
    count=count+1;
end
    
   st_walk_x(c+1)=x2;st_walk_y(c+1)=y2; st_walk_z(c+1)=z2;
   st_walk_xg((group-1)*sims+c+1)=x2;st_walk_yg((group-1)*sims+c+1)=y2; st_walk_zg((group-1)*sims+c+1)=z2;

end

for i=1:sims
disp1((group-1)*sims+i)=sqrt((st_walk_xg((group-1)*sims+i+1)-st_walk_xg((group-1)*sims+i))^2+....
    (st_walk_yg((group-1)*sims+i+1)-st_walk_yg((group-1)*sims+i))^2+...
    (st_walk_zg((group-1)*sims+i+1)-st_walk_zg((group-1)*sims+i))^2);
cell1((group-1)*sims+i)=group;
end

con=2;
for i=1:sims-con
disp3((group-1)*sims+i)=sqrt((st_walk_xg((group-1)*sims+i+1+con)-st_walk_xg((group-1)*sims+i))^2+....
    (st_walk_yg((group-1)*sims+i+1+con)-st_walk_yg((group-1)*sims+i))^2+...
    (st_walk_zg((group-1)*sims+i+1+con)-st_walk_zg((group-1)*sims+i))^2);
cell2((group-1)*sims+i)=group;
end

con=7;
for i=1:sims-con
disp8((group-1)*sims+i)=sqrt((st_walk_xg((group-1)*sims+i+1+con)-st_walk_xg((group-1)*sims+i))^2+....
    (st_walk_yg((group-1)*sims+i+1+con)-st_walk_yg((group-1)*sims+i))^2+...
    (st_walk_zg((group-1)*sims+i+1+con)-st_walk_zg((group-1)*sims+i))^2);
cell6((group-1)*sims+i)=group;
end

group=group+1;

end
w=size(skel,1);
l=size(skel,2); 
h=size(skel,3);
[x,y,z]=ind2sub([w,l,h],find(skel(:)));
 
 step_1=sims;% number of steps
sim=gp; %number of cells
figure(1); %clf
cou=1;
for i=1:gp
    plot3(st_walk_yg((cou-1)*step_1+1:cou*step_1),st_walk_xg((cou-1)*step_1+1:cou*step_1),st_walk_zg((cou-1)*step_1+1:cou*step_1),'MarkerSize',14,'MarkerFaceColor','g', 'LineWidth',2)
   
    xst=st_walk_yg((cou-1)*step_1+1:cou*step_1);
    yst=st_walk_xg((cou-1)*step_1+1:cou*step_1);
    zst=st_walk_zg((cou-1)*step_1+1:cou*step_1);
    
    xstm=xst-xst(1);
    ystm=yst-yst(1);
    zstm=zst-zst(1);
    
    pmsde=sqrt(xstm.^2+ystm.^2+zstm.^2);    
    pmsd(i,:)=pmsde(2:sims);
    hold on
    cou=cou+1;
end
grid on
view(150,75)

axis([1 512 1 512 1 44]);

pmsd=(300./512).*pmsd;
% msd2=log(mean(pmsd.^2)+0.0001);
msd2=(mean(pmsd.^2));

figure(3)
%    time2=log([1:sims-1]);
    time2=([1:sims-1]);

   
hold on; plot(time2',msd2','square') 
xlabel('log(Time delay)')
ylabel('log(MSD)')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Levy
clear all; hold on;%clf


gp=200; %cells
sims=50; %secs 1/2 hr


%load('skel.mat')
load('IMG_OPENED_3D_LN.mat');
skel_LN = Skeleton3D(IMG_OPENED_3D_LN); 

skel=skel_LN;

mean_R=2.1;
alphR=1.5;%alpha=4.5;  1
k_gpR=1/alphR;
xminR=mean_R*((alphR-1)/alphR);
sig_gpR=xminR/alphR; 
mean_PaR=alphR*xminR/(alphR-1);
sig_GR=mean_PaR*sqrt(pi()/2);
mean_GuR=sig_GR*sqrt(2/pi());
max_stepR=300;

[max_x,max_y,max_z]=size(skel);

vec=zeros(max_x*max_y,2);
for i=1:max_z
[row,col]=find(skel(:,:,i));
sz=size(row);
sz_st(i)=sz(1);
vec(1:sz(1),:,i)=[row,col];
end 
  
c_count=1;
group=1;        
for g=1:gp
    g
clear st_walk_x st_walk_y st_walk_z st_walk_x_f st_walk_y_f st_walk_z_f 
init_point=[0 0];
while init_point==[0 0]
dr=round(1+ (sz_st(1)-1).*rand());
dr_z=round(1+ (44-1).*rand());
init_point=vec(dr,:,dr_z);
x=init_point(1);xi=x;
y=init_point(2);yi=y;
z=dr_z;zi=z;
x2=x;y2=y;z2=z;
end
st_walk_x(1)=x; st_walk_y(1)=y; st_walk_z(1)=z;
st_walk_x_f(1)=x;st_walk_y_f(1)=y; st_walk_z_f(1)=z; 

st_walk_xg((group-1)*sims+1)=x; st_walk_yg((group-1)*sims+1)=y; st_walk_zg((group-1)*sims+1)=z;
st_walk_x_fg((group-1)*sims+1)=x;st_walk_y_fg((group-1)*sims+1)=y; st_walk_z_fg((group-1)*sims+1)=z; 

count=1;

for c=1:sims

brlen1R=gprnd(k_gpR,sig_gpR,sig_gpR/k_gpR,1);
step=min(brlen1R,max_stepR);

    
for m=1:step 
    clear xn xp yn yp zn zp
if x-1<1; xn=1;else; xn=x-1;end
if x+1>max_x; xp=x; else; xp=x+1; end
if y-1<1; yn=1; else; yn=y-1; end
if y+1>max_y; yp=y; else; yp=y+1; end
if z-1<1; zn=1; else; zn=z-1; end
if z+1>max_z; zp=z; else; zp=z+1; end

%st_i(1)=nan; st_j(1)=nan; st_k(1)=nan;
cou=1; clear st_i st_j st_k
    for k=zn:zp
        for j=yn:yp
            for i=xn:xp
                for p=1:max_x*max_y
                    if i==vec(p,1,k)  
                    if j==vec(p,2,k) 
                        if i==x & j==y & k==z
                        else
                    st_i(cou)=i; st_j(cou)=j; st_k(cou)=k;
                    cou=cou+1;
                        end
                    else
                    end
                    end
                end
            end
        end
    end
    st_fl(c)=cou;
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    flag=0;
    iter=0;
    while flag==0
    dr2=round(1+ (cou-1-1).*rand()); 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if cou>1
    x2=st_i(dr2); y2=st_j(dr2); z2=st_k(dr2);
    else
    x2=x2; y2=y2; z2=z2;
    end  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    x2=st_i(dr2); y2=st_j(dr2); z2=st_k(dr2);
        if st_walk_x_f(count)==x2 & st_walk_y_f(count)==y2 & st_walk_z_f(count)==z2
                 iter=iter+1;
        else if count>2 & ...
                x2== st_walk_x_f(count-1) & y2== st_walk_y_f(count-1) &  z2==st_walk_z_f(count-1)
                     iter=iter+1;
            else if count>3 & ...
                    x2== st_walk_x_f(count-2) & y2== st_walk_y_f(count-2) &  z2==st_walk_z_f(count-2)
                     iter=iter+1;                  
                    else 
                        flag=1;
                  end
             end
        end
              if iter>2
                  flag=1;
              end
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    x=x2;y=y2;z=z2;
    st_walk_x_f(count+1)=x2;st_walk_y_f(count+1)=y2; st_walk_z_f(count+1)=z2;
    st_walk_x_fg((group-1)*sims+count+1)=x2;st_walk_y_fg((group-1)*sims+count+1)=y2; st_walk_z_fg((group-1)*sims+count+1)=z2; 
    count=count+1;
end
    
   st_walk_x(c+1)=x2;st_walk_y(c+1)=y2; st_walk_z(c+1)=z2;
   st_walk_xg((group-1)*sims+c+1)=x2;st_walk_yg((group-1)*sims+c+1)=y2; st_walk_zg((group-1)*sims+c+1)=z2;

end

for i=1:sims
disp1((group-1)*sims+i)=sqrt((st_walk_xg((group-1)*sims+i+1)-st_walk_xg((group-1)*sims+i))^2+....
    (st_walk_yg((group-1)*sims+i+1)-st_walk_yg((group-1)*sims+i))^2+...
    (st_walk_zg((group-1)*sims+i+1)-st_walk_zg((group-1)*sims+i))^2);
cell1((group-1)*sims+i)=group;
end

con=2;
for i=1:sims-con
disp3((group-1)*sims+i)=sqrt((st_walk_xg((group-1)*sims+i+1+con)-st_walk_xg((group-1)*sims+i))^2+....
    (st_walk_yg((group-1)*sims+i+1+con)-st_walk_yg((group-1)*sims+i))^2+...
    (st_walk_zg((group-1)*sims+i+1+con)-st_walk_zg((group-1)*sims+i))^2);
cell2((group-1)*sims+i)=group;
end

con=7;
for i=1:sims-con
disp8((group-1)*sims+i)=sqrt((st_walk_xg((group-1)*sims+i+1+con)-st_walk_xg((group-1)*sims+i))^2+....
    (st_walk_yg((group-1)*sims+i+1+con)-st_walk_yg((group-1)*sims+i))^2+...
    (st_walk_zg((group-1)*sims+i+1+con)-st_walk_zg((group-1)*sims+i))^2);
cell6((group-1)*sims+i)=group;
end

group=group+1;

end
w=size(skel,1);
l=size(skel,2); 
h=size(skel,3);
[x,y,z]=ind2sub([w,l,h],find(skel(:)));
 
 step_1=sims;% number of steps
sim=gp; %number of cells
figure(2); %clf
cou=1;
for i=1:gp
    plot3(st_walk_yg((cou-1)*step_1+1:cou*step_1),st_walk_xg((cou-1)*step_1+1:cou*step_1),st_walk_zg((cou-1)*step_1+1:cou*step_1),'MarkerSize',14,'MarkerFaceColor','g', 'LineWidth',2)

    xst=st_walk_yg((cou-1)*step_1+1:cou*step_1);
    yst=st_walk_xg((cou-1)*step_1+1:cou*step_1);
    zst=st_walk_zg((cou-1)*step_1+1:cou*step_1);
    
    xstm=xst-xst(1);
    ystm=yst-yst(1);
    zstm=zst-zst(1);
    
    pmsde=sqrt(xstm.^2+ystm.^2+zstm.^2);    
    pmsd(i,:)=pmsde(2:sims);
    hold on
    cou=cou+1;
end
grid on
view(150,75)

axis([1 512 1 512 1 44]);


pmsd=(300./512).*pmsd;
% msd2=log(mean(pmsd.^2)+0.0001);
 msd2=(mean(pmsd.^2));

figure(3); hold on
%    time2=log([1:sims-1]);
    time2=([1:sims-1]);

    hold on; plot(time2',msd2','o') 
xlabel('log(Time delay)')
ylabel('log(MSD)')

legend('CRW','LW')


     
     

     
     

     
     