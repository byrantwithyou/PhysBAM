#!/bin/bash

WRAPPER_INC=-I./../../Public_Library -g -fPIC
WRAPPER_LIBS=-lPhysBAM_PhysBAM_Solids -lPhysBAM_PhysBAM_Tools -lPhysBAM_PhysBAM_Geometry -lPhysBAM_PhysBAM_Dynamics -lPhysBAM_PhysBAM_Fluids -lz -L../../build/nocona/debug/Public_Library

default : libPhysBAM_Wrapper.so test_sim

.PHONY : default

libmain.os : libmain.cpp
	g++ -c libmain.cpp -o libmain.os $(WRAPPER_INC)

DEFORMABLE_EXAMPLE.os : DEFORMABLE_EXAMPLE.cpp
	g++ -c DEFORMABLE_EXAMPLE.cpp -o DEFORMABLE_EXAMPLE.os $(WRAPPER_INC)

ACCESSORS.os : ACCESSORS.cpp
	g++ -c ACCESSORS.cpp -o ACCESSORS.os $(WRAPPER_INC)

libPhysBAM_Wrapper.so : DEFORMABLE_EXAMPLE.os libmain.os ACCESSORS.os
	g++ -shared DEFORMABLE_EXAMPLE.os libmain.os ACCESSORS.os -o libPhysBAM_Wrapper.so $(WRAPPER_INC) $(WRAPPER_LIBS)

test_sim : main.cpp
	g++ -o test_sim main.cpp -ldl -g

