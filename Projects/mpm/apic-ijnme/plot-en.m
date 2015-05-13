#!/usr/bin/octave -q

n=size(argv())(1);

legs={};
hold on
axis([0,5,0.0018,0.0021]);
for i = 1:n/3
  filename=argv(){3*i-2};
  leg=argv(){3*i-1};
  form=argv(){3*i};
  ["en-" filename ".txt"]
  x=load(["en-" filename ".txt"]);
  plot(x(:,1),x(:,2),form);
  legs{i}=leg;
end
legend(legs,"location","southeast");
xlabel('Time');
ylabel('Energy');
title('Total energy for varous schemes');
print -color -deps "all-en.pdf";
