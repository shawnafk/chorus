#!/bin/bash
cc=gnu
make -B comp=$cc -C mathlib/quadpack 
make -B comp=$cc -C mathlib/deriv 
make -B comp=$cc -C mathlib/minpack
