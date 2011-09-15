% particle_levelset_contour_2d('directory',filenumber)
function particle_levelset_contour_2d(directory,filenumber)
% get dimensions
filename=sprintf('%s/header.%d',directory,filenumber);fid=fopen(filename,'rb','l');
m=fread(fid,1,'int');n=fread(fid,1,'int');
x=fread(fid,[n,m],'double')';y=fread(fid,[n,m],'double')';
a=min(min(x));b=max(max(x));c=min(min(y));d=max(max(y));
% get level set
filename=sprintf('%s/levelset.%d',directory,filenumber);fid=fopen(filename,'rb','l');
phi=fread(fid,[n,m],'double');
% get particles
filename=sprintf('%s/positive_particles.%d',directory,filenumber);fid=fopen(filename,'rb','l');
number=fread(fid,1,'int');
positive_particle_x=fread(fid,number,'double');
positive_particle_y=fread(fid,number,'double');
filename=sprintf('%s/negative_particles.%d',directory,filenumber);fid=fopen(filename,'rb','l');
number=fread(fid,1,'int');
negative_particle_x=fread(fid,number,'double');
negative_particle_y=fread(fid,number,'double');
% get escaped particles
filename=sprintf('%s/escaped_positive_particles.%d',directory,filenumber);fid=fopen(filename,'rb','l');
number=fread(fid,1,'int');
escaped_positive_particle_x=fread(fid,number,'double');
escaped_positive_particle_y=fread(fid,number,'double');
filename=sprintf('%s/escaped_negative_particles.%d',directory,filenumber);fid=fopen(filename,'rb','l');
number=fread(fid,1,'int');
escaped_negative_particle_x=fread(fid,number,'double');
escaped_negative_particle_y=fread(fid,number,'double');
% draw graph
clf;orient tall;hold on;contour(x,y,phi',[0 0],'k');
axis equal;axis([a b c d]);
plot(positive_particle_x,positive_particle_y,'r.');
plot(negative_particle_x,negative_particle_y,'b.');
length(escaped_positive_particle_x)
length(escaped_negative_particle_x)
plot(escaped_positive_particle_x,escaped_positive_particle_y,'g.');
plot(escaped_negative_particle_x,escaped_negative_particle_y,'c.');
title('particle levelset')
% close files
fclose('all');
