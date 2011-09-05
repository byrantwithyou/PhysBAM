% levelset_movie_2d('directory',number_of_frames)
function levelset_movie_2d(directory,number_of_frames)

if (nargin < 1) directory = '.'; end
if (nargin < 2)
	% automatically calculate number of frames
	number_of_frames = 0;
	while (fopen(sprintf('%s/levelset.%d',directory,number_of_frames + 1))) >= 0,
		number_of_frames = number_of_frames + 1;
	end
end

% get dimensions
filename=sprintf('%s/header.%d',directory,0);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');n=fread(fid,1,'int');
% main loop
for frame_number=0:number_of_frames
   % get level set
   filename=sprintf('%s/levelset.%d',directory,frame_number);fid=fopen(filename,'rb','l');
	phi=fread(fid,[n,m],'double');
   % draw graph
   clf;orient tall;hold on;
	contour(phi,0,'b');
   axis equal;axis([1 m 1 n]);
   % make movie
 	M=getframe;
end
% close files
fclose('all');

