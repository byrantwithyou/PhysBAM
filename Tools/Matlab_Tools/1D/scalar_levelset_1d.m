% scalar_levelset_1d('directory',filenumber,'variable')
function scalar_levelset_1d(directory,filenumber,variable)
% get dimensions and grid
filename=sprintf('%s/header.%d',directory,filenumber);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');
x=fread(fid,[m],'double')';
a=min(min(x));b=max(max(x));
% get level set
filename=sprintf('%s/%s.%d',directory,variable,filenumber);fid=fopen(filename,'rb','l');
phi=fread(fid,[m],'double')';
% get variable
filename=sprintf('%s/%s.%d',directory,variable,filenumber);fid=fopen(filename,'rb','l');
u=fread(fid,[m],'double')';
minimum=min(min(u));maximum=max(max(u));difference=abs(maximum-minimum);
minimum=minimum-.1*max(difference,.0001*abs(minimum));maximum=maximum+.1*max(difference,.0001*abs(maximum));
if(minimum == 0 & maximum == 0) minimum=-.01;maximum=.01;end;
% set up different regions
x_phi=x;u_phi=u;
for i=1:m; 
   if(phi(i) > 0) x_phi(i)=NaN;
   else x(i)=NaN;
   end
end
% draw graph
clf;axis([a b minimum maximum]);hold on;
plot(x,u,'ro');
plot(x_phi,u_phi,'bo');
title(variable);
% close files
fclose('all'); 