%  Creates digital sinusoidal structure; 5120x512x44

clear all; clf
close all;
load('IMG_OPENED_3D_Evans_old.mat');
IMG=IMG_OPENED_3D_Evans_old;
skel = Skeleton3D(IMG); %testvol is the original 0-1 data
col=[.7 .7 .8];
hiso = patch(isosurface(IMG,0),'FaceColor',col,'EdgeColor','none');
axis equal;axis off;
lighting phong;
hold on
w=size(skel,1);
l=size(skel,2);
h=size(skel,3);
[x,y,z]=ind2sub([w,l,h],find(skel(:)));
plot3(y,x,z,'square','Markersize',4,'MarkerFaceColor','r','Color','r');            
set(gcf,'Color','white');
view(140,80)
isonormals(IMG,hiso);
alpha(0.5);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
hold on;







