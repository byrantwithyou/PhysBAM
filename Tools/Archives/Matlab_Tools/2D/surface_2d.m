% surface_2d('directory',filenumber,'variable')
function surface_2d(directory,filenumber,variable)
% get dimensions and grid
filename=sprintf('%s/header.%d',directory,filenumber);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');n=fread(fid,1,'int');
% get variable
filename=sprintf('%s/%s.%d',directory,variable,filenumber);fid=fopen(filename,'rb','l');
u=fread(fid,[n,m],'double')';
% draw graph
clf;
mesh(u');
title('surface');
axis equal;
colormap winter;
% close files
fclose('all'); 