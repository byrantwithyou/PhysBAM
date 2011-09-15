% scalar4_1d('directory',filenumber,'variable_1','variable_2','variable_3','variable_4')
function scalar4_1d(directory,filenumber,variable_1,variable_2,variable_3,variable_4)
% get dimensions and grid
filename=sprintf('%s/header.%d',directory,filenumber);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');
x=fread(fid,[m],'double')';
a=min(min(x));b=max(max(x));
% get level set
filename=sprintf('%s/levelset.%d',directory,filenumber);fid=fopen(filename,'rb','l');
phi=fread(fid,[m],'double')';
% main loop
clf;
for count=1:4
   if(count==1) variable=variable_1;
   elseif(count==2) variable=variable_2;
   elseif(count==3) variable=variable_3;
   elseif(count==4) variable=variable_4;
   end
   % get variable
   filename=sprintf('%s/%s.%d',directory,variable,filenumber);fid=fopen(filename,'rb','l');
  	u=fread(fid,[m],'double')';
	minimum=min(min(u));maximum=max(max(u));difference=abs(maximum-minimum);
   minimum=minimum-.1*max(difference,.0001*abs(minimum));maximum=maximum+.1*max(difference,.0001*abs(maximum));
   if(minimum == 0 & maximum == 0) minimum=-.01;maximum=.01;end;
   % set up different regions
   x_phi=x;u_phi=u;
   x_phi2=x;u_phi2=u;
	for i=1:m; 
   	if(phi(i) > 0) x_phi(i)=NaN;
   	else x_phi2(i)=NaN;
   	end
	end
	% draw graph
	subplot(2,2,count);axis([a b minimum maximum]);hold on;
	plot(x_phi2,u_phi2,'ro');
	plot(x_phi,u_phi,'bo');
   title(variable);
end
% close files
fclose('all'); 


