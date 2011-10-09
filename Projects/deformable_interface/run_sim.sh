#!/bin/bash

export LD_LIBRARY_PATH=".:/home/cas43/PhysBAM-2011//build/nocona/debug/Public_Library"

LD_PRELOAD=libPhysBAM_Wrapper.so LD_LIBRARY_PATH=".:/home/cas43/PhysBAM-2011//build/nocona/debug/Public_Library" gdb --args test_sim

