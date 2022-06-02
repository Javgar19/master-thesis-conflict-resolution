#!/bin/bash

radius = {0.5, 1, 1.5}
n_ac = {50, 70, 100}
threshold = {35, 25, 20, 15}

for r in radius
do
	for n in n_ac
	do
		for thr in threshold
		do
			echo `date +"%Y-%m-%d %T"` "with r=${r} n=${n} thr=${thr}"
			python sacre_sim.py with radius=r n_ac=n -m db_test