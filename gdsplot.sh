#!/bin/bash
#./GDT-4.0.4/gds2gdt.Linux ./4004.gds ./4004.GDT
#sleep .1s
#mkdir 4004_recreate
./TestInterface ./nand2.gdt
sleep .5s 
gnuplot nand2_recreate.gnu
