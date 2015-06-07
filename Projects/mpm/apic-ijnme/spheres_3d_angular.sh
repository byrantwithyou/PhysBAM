#!/bin/bash

NAME=`basename $0 | sed 's/\.sh$//'`
 ./mpm -3d 11 -print_stats -resolution 256 -framerate 12 -last_frame 250 -affine -midpoint -w1 1 -w2 1 -scale_E 8 

for i in {apic}{,} ; do
    grep 'particle state angular' $NAME-$i/common/log.txt | awk '{print $5 " " $7}' | sed 's/[()]//g' > am-$i.txt
    grep 'after particle to grid particle total energy' $NAME-$i/common/log.txt | awk '{print $9 " " $11}' | sed 's/<.*//g' > en-p-$i.txt
    grep 'after particle to grid total energy' $NAME-$i/common/log.txt | awk '{print $8 " " $10}' | sed 's/<.*//g' > en-g-$i.txt
done

cat <<EOF > st.m
n=size(argv())(1);

legs={};
hold on
# axis([0,20,-0.02,0.035]);
for i = 1:n/3
  filename=argv(){3*i-2};
  leg=argv(){3*i-1};
  form=argv(){3*i};
  x=load(["am-" filename ".txt"]);
  plot(x(:,1),x(:,2),form);
  legs{i}=leg;
end
legend(legs,"location","northeast");
xlabel('Time');
ylabel('Angular momentum');
title('Angular momentum for varous schemes');
print -color -deps "all-am-$NAME.pdf";
EOF
octave -q st.m \
apic "APIC Midpoint Rule"  r- \

cat <<EOF > st.m
n=size(argv())(1);

legs={};
hold on
# axis([0,20,60,95]);
for i = 1:n/3
  filename=argv(){3*i-2};
  leg=argv(){3*i-1};
  form=argv(){3*i};
  ["en-" filename ".txt"]
  x=load(["en-" filename ".txt"]);
  plot(x(:,1),x(:,2),form);
  legs{i}=leg;
end
legend(legs,"location","southwest");
xlabel('Time');
ylabel('Energy');
title('Total energy for varous schemes');
print -color -deps "all-en-$NAME.pdf";
EOF
octave -q st.m \
g-apic "APIC Midpoint Rule"  r- \
