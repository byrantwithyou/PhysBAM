% levelset_object_contour_2d('directory',filenumber)
function levelset_object_contour_2d(directory,filenumber)
% get dimensions
filename=sprintf('%s/header.%d',directory,filenumber);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');n=fread(fid,1,'int');
% get level set
filename=sprintf('%s/levelset.%d',directory,filenumber);fid=fopen(filename,'rb','l');
phi=fread(fid,[n,m],'double');
% get object
filename=sprintf('%s/object.%d',directory,filenumber);fid=fopen(filename,'rb','l');
psi=fread(fid,[n,m],'double');
% draw graph
clf;orient tall;hold on;
contour(phi,0,'b')
contour(psi,0,'r')
axis equal;axis([1 m 1 n])
title('levelset and object')
% close files
fclose('all'); 