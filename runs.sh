#!/bin/bash

radiuses = {0.5, 1, 1.5}
n_acs = {50, 70, 100}
thresholds = {35, 25, 20, 15}

for r in ${radiuses[@]}
do
	for n in ${n_acs[@]}
	do
		for thr in ${thresholds[@]}
		do
			echo `date +"%Y-%m-%d %T"` "with r=${r} n=${n} thr=${thr}"
			python sacre_sim.py with radius=$r n_ac=$n tcpa_thresh=$thr -m db_test