#!/bin/bash
for file in $(ls plot/*.gp); do
  gnuplot $file
done 