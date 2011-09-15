% velocity_movie_2d('directory',arrrow_size,number_of_frames)
function velocity_movie_2d(directory,arrow_size,number_of_frames)
% get dimensions and grid
filename=sprintf('%s/header.%d',directory,0);fid=fopen(filename,'rb');
m=fread(fid,1,'int');n=fread(fid,1,'int');
x=fread(fid,m*n,'double');y=fread(fid,m*n,'double');
a=min(min(x));b=max(max(x));c=min(min(y));d=max(max(y));
% main loop 
for frame_number=0:number_of_frames   
   % get velocity
	filename=sprintf('%s/velocity1.%d',directory,frame_number);fid=fopen(filename,'rb');
	u=fread(fid,m*n,'double');
	filename=sprintf('%s/velocity2.%d',directory,frame_number);fid=fopen(filename,'rb');
	v=fread(fid,m*n,'double');
	% draw graph
	clf;hold on;
	quiver_ron(x,y,u,v,arrow_size,'b');
   axis([a b c d]);axis image;
   % make movie
 	M=getframe;
end
% close files
fclose('all'); 