% scalar_levelset_slice_2d('directory','xy_constant',slice_number,filenumber,'variable_1','variable_2','variable_3','variable_4')
function scalar_levelset_slice_2d(directory,xy_constant,slice_number,filenumber,variable_1,variable_2,variable_3,variable_4)
% get dimensions and grid
filename=sprintf('%s/header.%d',directory,filenumber);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');n=fread(fid,1,'int');
x=fread(fid,[n,m],'double')';y=fread(fid,[n,m],'double')';
a=min(min(x));b=max(max(x));c=min(min(y));d=max(max(y));
% get level set
filename=sprintf('%s/levelset.%d',directory,filenumber);fid=fopen(filename,'rb','l');
phi=fread(fid,[n,m],'double')';
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
	u=fread(fid,[n,m],'double')';
	minimum=min(min(u));maximum=max(max(u));difference=abs(maximum-minimum);
   minimum=minimum-.1*max(difference,.0001*abs(minimum));maximum=maximum+.1*max(difference,.0001*abs(maximum));
   if(minimum == 0 & maximum == 0) minimum=-.01;maximum=.01;end;
	% set up different regions
	x_phi=x;y_phi=y;u_phi=u;
	x_phi2=x;y_phi2=y;u_phi2=u;
	for i=1:m; for j=1:n; 
   	if(phi(i,j) > 0) x_phi(i,j)=NaN;y_phi(i,j)=NaN;
   	else x_phi2(i,j)=NaN;y_phi2(i,j)=NaN;
   	end
   end;end
	% draw graph
	subplot(2,2,count);	
	if(xy_constant == 'x') 
   	for j=1:n;
      	u_1d(j)=u_phi2(slice_number,j);y_1d(j)=y_phi2(slice_number,j);
      	u_phi_1d(j)=u_phi(slice_number,j);y_phi_1d(j)=y_phi(slice_number,j);
   	end;
   	axis([c d minimum maximum]);hold on;plot(y_1d,u_1d,'ro');plot(y_phi_1d,u_phi_1d,'bo');
	elseif(xy_constant == 'y') 
   	for i=1:m;
      	u_1d(i)=u(i,slice_number);x_1d(i)=x(i,slice_number);
      	u_phi_1d(i)=u_phi(i,slice_number);x_phi_1d(i)=x_phi(i,slice_number);
   	end;
 		axis([a b minimum maximum]);hold on;plot(x_1d,u_1d,'ro');plot(x_phi_1d,u_phi_1d,'bo');
	end
   title(variable);
end
% close files
fclose('all'); 
   

