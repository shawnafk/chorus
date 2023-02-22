#!/bin/bash
cc=intel
make -C mathlib/quadpack clean 
make -C mathlib/deriv clean
make -C mathlib/minpack clean
make clean
