#!/bin/bash

g++ -fPIC libmain.cpp -I./../../Public_Library -c -o libmain.os -g
g++ -fPIC DEFORMABLE_EXAMPLE.cpp -I./../../Public_Library -c -o DEFORMABLE_EXAMPLE.os -g
g++ -shared DEFORMABLE_EXAMPLE.os libmain.os -o libPhysBAM_Wrapper.so -lPhysBAM -L../../build/nocona/debug/Public_Library -g

g++ -o test_sim main.cpp -ldl -g

