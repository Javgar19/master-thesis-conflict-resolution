import os
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()


parser.add_argument(
	"-file",
	default = "main.c",
)

parser.add_argument(
	"-o",
	default="main.exe"
)

parser.add_argument(
	"-results"
	default="output.txt"
)


if __name__ == "__main__":
	args = parser.parse_args()

	os.system(f"gcc -Wall -O3 -o {args.o} {args.file} RKF78-2.2.c/RKF78.c -lm")
	os.system(f"./{args.o}")

	xt = []
	t = range(11)
	with open(args.results) as f:
		
