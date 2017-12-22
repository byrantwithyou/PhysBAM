#!/bin/bash

NAME=results

restitutions=("0" "0.5" "1")
frictions=("0" "0.5" "1")
pattern="'Objective...|Number of Iterations|sol:|EXIT|initial_guess:|velocity:|friction:|restitution:|ERROR:'"

velocities=('"0 -1"' '"0 0"' '"1 0"' '"1 -1"')
initial_guesses=('"0 -1"' '"0 -1"' '"0 -1"' '"1 -1"')
prefices=("n" "zero" "t" "inclined")

rm -rf $NAME
mkdir -p $NAME

for t in `seq 0 $((${#prefices[@]}-1))` ; do
    velocity=${velocities[$t]}
    initial_guess=${initial_guesses[$t]}
    ARGS="./opt_proj -s $initial_guess -v $velocity"
    for rs in `seq 0 $((${#restitutions[@]}-1))` ; do
        for f in `seq 0 $((${#frictions[@]}-1))` ; do
            output=$NAME/${prefices[$t]}-f${frictions[$f]}-rs${restitutions[$rs]}.txt
            echo "$ARGS -friction ${frictions[$f]} -restitution ${restitutions[$rs]} | egrep $pattern > $output"
        done
    done
done | xargs -P 16 -n 1 -d '\n' bash -c

