tmax=2;
samples=1000;

ZERO=zeros(1,samples);
ONE=ones(1,samples);

T=(1/(samples-1))*(0:(samples-1));
track.time=tmax*T;
track.trajectory.t=zeros(3,samples);

%track.trajectory.r=quaternion_from_euler([pi/2*ONE;pi/4*cos(2*pi*T);ZERO]);
track.trajectory.r=quatmult(quaternion_from_euler([3*pi/8*ONE;ZERO;pi/8*ONE]),quaternion_from_euler([ZERO;ZERO;pi/4*cos(2*pi*T)]));
write_track('ChainSwing/glenohumeral.track',track);

%track.trajectory.r=quaternion_from_euler([sin(2*pi*track.time/tmax);ZERO;ZERO]);
track.trajectory.r=quaternion_from_euler([pi/2+(pi/4)*sin(2*pi*T);ZERO;ZERO]);
write_track('ChainSwing/humeroulnar.track',track);

%track.trajectory.r=quaternion_from_euler([cos(2*pi*track.time/tmax);ZERO;ZERO]);
track.trajectory.r=quaternion_from_euler([pi/2*ONE;ZERO;ZERO]);
write_track('ChainSwing/radioulnar.track',track);

track.trajectory.r=quaternion_from_euler([ZERO;0.4*ONE;ZERO]);
write_track('ChainSwing/hand.track',track);
