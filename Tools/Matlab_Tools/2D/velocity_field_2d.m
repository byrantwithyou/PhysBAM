% velocity_field_2d('directory',arrow_size,filenumber)
function velocity_field_2d(directory,arrow_size,filenumber)
% get dimensions and grid
filename=sprintf('%s/header.%d',directory,filenumber);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');n=fread(fid,1,'int');
x=fread(fid,m*n,'double');y=fread(fid,m*n,'double');
a=min(min(x));b=max(max(x));c=min(min(y));d=max(max(y));
% get velocity
filename=sprintf('%s/velocity1.%d',directory,filenumber);fid=fopen(filename,'rb','l');
u=fread(fid,m*n,'double');
filename=sprintf('%s/velocity2.%d',directory,filenumber);fid=fopen(filename,'rb','l');
v=fread(fid,m*n,'double');
% draw graph
clf;hold on;
quiver_ron(x,y,u,v,arrow_size,'b');
axis([a b c d]);axis image;
title('velocity field');
% close files
fclose('all'); 