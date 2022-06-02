import os

radiuses = [1.5]
n_acs = [50, 70, 100]
threshold = [35.0, 25.0, 20.0, 15.0]

for r in radiuses:
    for n in n_acs:
        for thr in threshold:
            print(f"radius={r} n_ac={n} tcpa_thresh={thr}")
            os.system(f"python sacred_sim.py with radius={r} n_ac={n} tcpa_thresh={thr} -F results")