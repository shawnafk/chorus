#!/bin/bash
cc=gnu
make -B comp=$cc -C mathlib/quadpack 
make -B comp=$cc -C mathlib/deriv 
make -Bn comp=$cc -C mathlib/minpack
make -B comp=$cc
