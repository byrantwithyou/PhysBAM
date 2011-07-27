% levelset_object_3d('directory',filenumber)
function levelset_object_3d(directory,filenumber)
% get dimensions
filename=sprintf('%s/header.%d',directory,filenumber);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');n=fread(fid,1,'int');mn=fread(fid,1,'int');
% get level set
filename=sprintf('%s/levelset.%d',directory,filenumber);fid=fopen(filename,'rb','l');
phi=zeros(m,n,mn);for i=1:m;phi(i,:,:)=fread(fid,[mn,n],'double')';end
% get object
filename=sprintf('%s/object.%d',directory,filenumber);fid=fopen(filename,'rb','l');
psi=zeros(m,n,mn);for i=1:m;psi(i,:,:)=fread(fid,[mn,n],'double')';end
% draw graph
phi_new=zeros(mn,m,n);for i=1:m;for j=1:n;for ij=1:mn;phi_new(ij,i,j)=phi(i,j,ij);end;end;end;
psi_new=zeros(mn,m,n);for i=1:m;for j=1:n;for ij=1:mn;psi_new(ij,i,j)=psi(i,j,ij);end;end;end;
clf;orient tall;
axis([1 m 1 mn 1 n]);hold on;axis equal;view(0,10);
isonormals(phi_new,patch(isosurface(phi_new,0),'FaceColor','blue','EdgeColor','none'))
isonormals(psi_new,patch(isosurface(psi_new,0),'FaceColor','red','EdgeColor','none'))
camlight;camlight(-80,-10);lighting phong;
title('levelset and object');   
% close files
fclose('all'); 