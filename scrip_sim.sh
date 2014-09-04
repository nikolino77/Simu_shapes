#!/bin/bash
for i1 in 1 3 5 10
do
    for i2 in 10e-9 40e-09 100e-09 300e-09
    do
      for i3 in 0 10e-12 50e-12 100e-12 300e-12
      do
        ./temp 1000 $i1 $i2 $i3 66e-12 10e-12 0  
      done
    done
done