#!/usr/bin/env python

import os
import sys
if len(sys.argv)==2:
	j=int(sys.argv[1])
else:
	j=0
i=0
k=3

dir = "./"
program = ("galsim","omp_gal")
N = (10,100,500,1000,2000)
filedir = "input_data/"
file = ("ellipse_N_00010.gal","ellipse_N_00100.gal","ellipse_N_00500.gal",
	"ellipse_N_01000.gal","ellipse_N_02000.gal")
nsteps = (1,2,10,200)
dt = 1e-05
theta_max = 0
graphics=0

os.system("clear; make "+program[j])
#exit()
cmd = dir+program[j]+" "+str(N[i])+" "+filedir+file[i]+" "+str(nsteps[k])+" "+str(dt)
cmd = cmd+" "+str(theta_max)+" "+str(graphics)
os.system(cmd)

compare = "compare_gal_files"
cmpdir = "compare_gal_files/"
refdata = ("ellipse_N_00010_after200steps.gal","ellipse_N_00100_after200steps.gal",
	"ellipse_N_00500_after200steps.gal","ellipse_N_01000_after200steps.gal",
	"ellipse_N_02000_after200steps.gal")
refdir = "ref_output_data/"

if nsteps[k]!=200:
	exit()

cmd = cmpdir+compare+" "+str(N[i])+" "+"result.gal"+" "+refdir+refdata[i]
os.system(cmd)
