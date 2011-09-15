% levelset__normals_2d('directory',arrow_size,filenumber)
function levelset_normals_2d(directory,arrow_size,filenumber)
% get dimensions and grid
filename=sprintf('%s/header.%d',directory,filenumber);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');n=fread(fid,1,'int');
x=fread(fid,m*n,'double');y=fread(fid,m*n,'double');
a=min(min(x));b=max(max(x));c=min(min(y));d=max(max(y));
% get level set
filename=sprintf('%s/levelset.%d',directory,filenumber);fid=fopen(filename,'rb','l');
phi=fread(fid,m*n,'double');
% get velocity
filename=sprintf('%s/n1.%d',directory,filenumber);fid=fopen(filename,'rb','l');
n1=fread(fid,m*n,'double');
filename=sprintf('%s/n2.%d',directory,filenumber);fid=fopen(filename,'rb','l');
n2=fread(fid,m*n,'double');
% set up different regions
x_phi=x;y_phi=y;n1_phi=n1;n2_phi=n2;
for k=1:m*n;  
   if(phi(k) > 0) x_phi(k)=NaN;y_phi(k)=NaN;
   else x(k)=NaN;y(k)=NaN;
   end
end
% draw graph
clf;hold on;
quiver_ron(x,y,n1,n2,arrow_size,'r');
quiver_ron(x_phi,y_phi,n1_phi,n2_phi,arrow_size,'b');
axis([a b c d]);axis image;
title('velocity field')
% close files
fclose('all'); 