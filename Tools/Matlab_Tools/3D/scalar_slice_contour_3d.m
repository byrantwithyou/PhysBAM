% scalar_slice_contour_3d('directory','xyz_constant',slice_number,filenumber,'variable')
function scalar_slice_contour_3d(directory,xyz_constant,slice_number,filenumber,variable)
% get dimensions
filename=sprintf('%s/header.%d',directory,filenumber);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');n=fread(fid,1,'int');mn=fread(fid,1,'int');
% get variable
filename=sprintf('%s/%s.%d',directory,variable,filenumber);fid=fopen(filename,'rb','l');
if(xyz_constant=='x') fseek(fid,8*mn*n*(slice_number-1),'bof');u=fread(fid,[mn,n],'double');
elseif(xyz_constant=='y') u=zeros(m,mn);fseek(fid,8*mn*(slice_number-1),'bof');
   for i=1:m
      u(i,:)=fread(fid,mn,'double')';
      fseek(fid,8*mn*(n-1),'cof');
   end
elseif(xyz_constant=='z') fseek(fid,8*(slice_number-1),'bof');u=fread(fid,[n,m],'double',8*(mn-1));   
end
% draw graph
clf;orient tall;hold on;
if(xyz_constant=='x') 
   contourf(u');axis equal;axis([1 mn 1 n]);
   titlewords=sprintf('levelset zy cross section x=%d',slice_number);title(titlewords);   
elseif(xyz_constant=='y') 
   contourf(u);axis equal;axis([1 m 1 mn]);
   titlewords=sprintf('levelset xz cross section y=%d',slice_number);title(titlewords);
elseif(xyz_constant=='z') 
   contourf(u);axis equal;axis([1 m 1 n]);
   titlewords=sprintf('levelset xy cross section z=%d',slice_number);title(titlewords);
end  
% close files
fclose('all'); 