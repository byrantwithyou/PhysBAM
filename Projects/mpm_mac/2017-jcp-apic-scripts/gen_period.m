function gen_period = gen_period(name,start)
a = importdata(strcat(name,'-time'));
b = importdata(strcat(name,'-vort'));
a = a(start:end);
b = b(start:end);
z = data_zeros(a,b);
r1 = z(2:end)-z(1:end-1);
r = sort(r1(1:end-1)+r1(2:end));
dlmwrite(strcat(name,'-period'),r);
