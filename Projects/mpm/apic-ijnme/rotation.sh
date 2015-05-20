#!/bin/bash

NAME=`basename $0 | sed 's/\.sh$//'`
ARGS="1 -max_dt 1e-2 -newton_tolerance 1e-4 -regular_seeding -print_stats -last_frame 120"
if : ; then
../mpm -affine -midpoint -o $NAME-apic $ARGS >/dev/null &
../mpm -midpoint -o $NAME-pic $ARGS >/dev/null &
../mpm -flip 1 -midpoint -o $NAME-flip $ARGS >/dev/null &

../mpm -affine -o $NAME-apic-be $ARGS >/dev/null &
../mpm -o $NAME-pic-be $ARGS >/dev/null &
../mpm -flip 1 -o $NAME-flip-be $ARGS >/dev/null &

../mpm -affine -symplectic_euler -o $NAME-apic-symp $ARGS -max_dt 1e-3 >/dev/null &
../mpm -symplectic_euler -o $NAME-pic-symp $ARGS -max_dt 1e-3 >/dev/null &
../mpm -flip 1 -symplectic_euler -o $NAME-flip-symp $ARGS -max_dt 1e-3 >/dev/null &
wait
fi

for i in {apic,pic,flip}{,-be,-symp} ; do
    grep 'particle state angular' $NAME-$i/common/log.txt | awk '{print $5 " " $7}' | sed 's/[()]//g' > am-$i.txt
    grep 'after particle to grid particle total energy' $NAME-$i/common/log.txt | awk '{print $9 " " $11}' | sed 's/<.*//g' > en-p-$i.txt
    grep 'after particle to grid total energy' $NAME-$i/common/log.txt | awk '{print $8 " " $10}' | sed 's/<.*//g' > en-g-$i.txt
done

cat <<EOF > st.m
n=size(argv())(1);

legs={};
hold on
axis([0,5,0.00980,0.01035]);
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
print -color -deps "all-am-$NAME.pdf";
EOF
octave -q st.m \
apic "APIC Midpoint Rule"  r- \
apic-be "APIC Backward Euler"  r-- \
apic-symp "APIC Symplectic Euler"  r-. \
pic "PIC Midpoint Rule"  g- \
pic-be "PIC Backward Euler"  g-- \
pic-symp "PIC Symplectic Euler"  g-. \
flip "FLIP Midpoint Rule"  b- \
flip-be "FLIP Backward Euler"  b-- \
flip-symp "FLIP Symplectic Euler"  b-.

cat <<EOF > st.m
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
print -color -deps "all-en-$NAME.pdf";
EOF
octave -q st.m \
g-apic "APIC Midpoint Rule"  r- \
g-apic-be "APIC Backward Euler"  r-- \
g-apic-symp "APIC Symplectic Euler"  r-. \
g-pic "PIC Midpoint Rule"  g- \
g-pic-be "PIC Backward Euler"  g-- \
g-pic-symp "PIC Symplectic Euler"  g-. \
g-flip "FLIP Midpoint Rule (grid)"  b- \
g-flip-be "FLIP Backward Euler (grid)"  b-- \
g-flip-symp "FLIP Symplectic Euler (grid)"  b-. \
p-flip "FLIP Midpoint Rule (particle)"  k- \
p-flip-be "FLIP Backward Euler (particle)"  k-- \
p-flip-symp "FLIP Symplectic Euler (particle)"  k-.
