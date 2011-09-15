% scalar_levelset_2d('directory',filenumber,'variable')
function scalar_levelset_2d(directory,filenumber,variable)
% get dimensions and grid
filename=sprintf('%s/header.%d',directory,filenumber);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');n=fread(fid,1,'int');
x=fread(fid,[n,m],'double')';y=fread(fid,[n,m],'double')';
a=min(min(x));b=max(max(x));c=min(min(y));d=max(max(y));
% get level set
filename=sprintf('%s/levelset.%d',directory,filenumber);fid=fopen(filename,'rb','l');
phi=fread(fid,[n,m],'double')';
% get variable
filename=sprintf('%s/%s.%d',directory,variable,filenumber);fid=fopen(filename,'rb','l');
u=fread(fid,[n,m],'double')';
minimum=min(min(u));maximum=max(max(u));difference=abs(maximum-minimum);
minimum=minimum-.1*max(difference,.0001*abs(minimum));maximum=maximum+.1*max(difference,.0001*abs(maximum));
if(minimum == 0 & maximum == 0) minimum=-.01;maximum=.01;end;
% set up different regions
x_phi=x;y_phi=y;u_phi=u;
for i=1:m; for j=1:n; 
   if(phi(i,j) > 0) x_phi(i,j)=NaN;y_phi(i,j)=NaN;
   else x(i,j)=NaN;y(i,j)=NaN;
   end
end;end
% draw graph
clf;axis([a b c d minimum maximum]);hold on;
plot3(x,y,u,'ro');
plot3(x_phi,y_phi,u_phi,'bo');
title(variable);
% close files
fclose('all'); 
