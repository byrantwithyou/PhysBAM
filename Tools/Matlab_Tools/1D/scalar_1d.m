% scalar_1d('directory',filenumber,'variable')
function scalar_1d(directory,filenumber,variable,dont_clear)
% get dimensions and grid
filename=sprintf('%s/header.%d',directory,filenumber);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');
x=fread(fid,[m],'double')';
a=min(min(x));b=max(max(x));
% get variable
filename=sprintf('%s/%s.%d',directory,variable,filenumber);fid=fopen(filename,'rb','l');
u=fread(fid,[m],'double')';
minimum=min(min(u));maximum=max(max(u));difference=abs(maximum-minimum);
minimum=minimum-.1*max(difference,.0001*abs(minimum));maximum=maximum+.1*max(difference,.0001*abs(maximum));
if(minimum == 0 & maximum == 0) minimum=-.01;maximum=.01;end;
% draw graph
if nargin~=4
    clf;
end
axis([a b minimum maximum]);hold on;
plot(x,u,'bo');
title(variable);
% close files
fclose('all'); 