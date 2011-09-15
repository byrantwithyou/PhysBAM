% scalar_1d_movie('directory','variable',number_of_frames)
function scalar_1d(directory,variable,number_of_frames)

% get dimensions and grid
filename=sprintf('%s/header.%d',directory,0);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');
x=fread(fid,[m],'double')';
a=min(min(x));b=max(max(x));

frame_number=0;
filename=sprintf('%s/%s.%d',directory,variable,frame_number);fid=fopen(filename,'rb','l');
u=fread(fid,[m],'double')';
minimum=min(min(u));maximum=max(max(u));difference=abs(maximum-minimum);
minimum=minimum-.2*max(difference,.0001*abs(minimum));maximum=maximum+.2*max(difference,.0001*abs(maximum));
if(minimum == 0 & maximum == 0) minimum=-.01;maximum=.01;end;    

minimum=-300;maximum=300;

for frame_number=0:number_of_frames
    % get variable
    filename=sprintf('%s/%s.%d',directory,variable,frame_number);fid=fopen(filename,'rb','l');
    u=fread(fid,[m],'double')';
    % draw graph
    clf;axis([a b minimum maximum]);hold on;
    plot(x,u,'bo');
    %axis equal;
    % make movie
 	M=getframe;
end

% close files
fclose('all'); 