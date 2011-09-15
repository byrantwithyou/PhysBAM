#!/usr/bin/python

import os;
import sys;
import re;

step="39"
# 25 29

step=sys.argv[2]

#cmd="grep 'Step %s' %s -A 1"%(step,sys.argv[1])
cmd="cat %s"%(sys.argv[1])

stdin,stdout=os.popen2(cmd)
energy_regex = re.compile('total energy = [0-9\.e\-]+')
time_regex = re.compile('\s[0-9\.e\-]+')
number_regex = re.compile('[0-9\.e\-]+$')
check_time_regex = re.compile("<print>time")
step_regex = re.compile('Step %s'%step)

empty_line="seed"
energy=0.0
time=0.0
while empty_line:
    early_exit=0
    # find step line
    energy_line=stdout.readline()
    while step_regex.search(energy_line)==None:
        early_exit=early_exit+1
        if early_exit>1000:
            sys.exit(0)
        energy_line=stdout.readline()
    time_line=stdout.readline()
    while check_time_regex.search(time_line)==None:
        early_exit=early_exit+1
        if early_exit>1000:
            sys.exit(0)
        time_line=stdout.readline()
    empty_line=stdout.readline()
    energy_extract=energy_regex.search(energy_line)
    if energy_extract:
        energy_number_group=energy_extract.group()
        energy_number=number_regex.search(energy_number_group).group()
        energy=float(energy_number)
    time_extract=time_regex.search(time_line)
    if time_extract:
        time_number=time_extract.group()
        time=float(time_number)
    #energy=float(number_regex.search(energy_extract).group())
    #time=float(number_regex.search(time_line.group()))
    print "%f %f"%(time,energy)
