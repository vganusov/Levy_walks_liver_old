% Creates digital sinusoidal structure; 100x100x44 

clear all; clf
close all;
load('IMG_OPENED_3D_LN');
IMG=IMG_OPENED_3D_LN(1:100,1:100,:);
skel = Skeleton3D(IMG);
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
view(140,50)
isonormals(IMG,hiso);
alpha(0.5);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
hold on;







