% object_contour_2d('directory',filenumber)
function object_contour_2d(directory,filenumber)
% get dimensions
filename=sprintf('%s/header.%d',directory,filenumber);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');n=fread(fid,1,'int');
% get level set
filename=sprintf('%s/object.%d',directory,filenumber);fid=fopen(filename,'rb','l');
phi=fread(fid,[n,m],'double');
% draw graph
clf;orient tall;hold on;
contour(phi,0,'r')
axis equal;axis([1 m 1 n])
title('object')
% close files
fclose('all'); 