% Creates the digital graph; 100x100x44 

clear all;
close all;
load('IMG_OPENED_3D_Evans_old.mat');
skel = Skeleton3D(IMG_OPENED_3D_Evans_old); 
skel_2=(skel(1:100,1:100,:));
clear skel
skel=skel_2;
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
figure();
hold on;
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
            line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',2);
        end;
    end;

    plot3(y1,x1,z1,'o','Markersize',9,...
        'MarkerFaceColor',ncol,...
        'Color','k');
end;
axis image;axis off;
set(gcf,'Color','white');
drawnow;
view(140,50)


