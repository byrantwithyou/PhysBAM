% levelset_object_movie_2d('directory',number_of_frames)
function levelset_object_movie_2d(directory,number_of_frames)
colormap winter
% get dimensions
filename=sprintf('%s/header.%d',directory,0);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');n=fread(fid,1,'int');
% main loop
for frame_number=0:number_of_frames
   % get level set
   filename=sprintf('%s/levelset.%d',directory,frame_number);fid=fopen(filename,'rb','l');
	phi=fread(fid,[n,m],'double');
   % get object
	filename=sprintf('%s/object.%d',directory,frame_number);fid=fopen(filename,'rb','l');
	psi=fread(fid,[n,m],'double');
	% draw graph
   clf;orient tall;hold on;
	contour(phi,0,'b');
   contour(psi,0,'r');
  	axis equal;axis([1 m 1 n]);
   % make movie
 	M=getframe;
end
% close files
fclose('all'); 

