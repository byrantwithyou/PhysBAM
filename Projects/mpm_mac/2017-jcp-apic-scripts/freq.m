function freq=freq(name,start)
a=importdata(strcat(name,'-time'));
b=importdata(strcat(name,'-vort'));
a=a(start:end);
b=b(start:end);
dts=a(2:end)-a(1:end-1);
dtsnz=dts(dts>0);
dt=min(dtsnz);
t=a(1):dt:a(end);
v=interp1(a,b,t);
sp=abs(fft(v));
[~,idx]=max(sp(1:end/2));
f=(idx-1)/dt/size(t,2);
fprintf('Maximum occurs at %d Hz.\n',f);
fprintf('Period: %f.\n',1/f);
