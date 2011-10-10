#!/bin/bash

g++ -fPIC libmain.cpp -I./../../Public_Library -c -o libmain.os -g
g++ -fPIC DEFORMABLE_EXAMPLE.cpp -I./../../Public_Library -c -o DEFORMABLE_EXAMPLE.os -g
g++ -fPIC ACCESSORS.cpp -I./../../Public_Library -c -o ACCESSORS.os -g
g++ -shared DEFORMABLE_EXAMPLE.os libmain.os ACCESSORS.os -o libPhysBAM_Wrapper.so -lPhysBAM_PhysBAM_Solids -lPhysBAM_PhysBAM_Tools -lPhysBAM_PhysBAM_Geometry -lPhysBAM_PhysBAM_Dynamics -lPhysBAM_PhysBAM_Fluids -lz -L../../build/nocona/debug/Public_Library -g
g++ -o test_sim main.cpp -ldl -g

