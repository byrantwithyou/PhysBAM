#!/usr/bin/octave -q

n=size(argv())(1);

legs={};
hold on
for i = 1:n/3
  filename=argv(){3*i-2};
  leg=argv(){3*i-1};
  form=argv(){3*i};
  x=load(["am-" filename ".txt"]);
  plot(x(:,1),x(:,2),form);
  legs{i}=leg;
end
legend(legs,"location","southeast");
xlabel('Time');
ylabel('Angular momentum');
title('Angular momentum for varous schemes');
print -color -deps "all-am.pdf";
