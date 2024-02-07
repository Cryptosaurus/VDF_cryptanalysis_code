# Smoothness simulation

This code evaluates the probability of almost-smoothness, by
generating random factored numbers, and counting how many are
almost-smooth.

It uses Bach's algorithm implementation from:
https://github.com/ethereum/research/blob/master/rsa_moduli/bach_random_factored_numbers.py

Requirements:

	pip install git+https://github.com/elliptic-shiho/primefac-fork@master

Example:
To evaluate the probability of a 256-bit to be (2^32, 2^65)-almost-smooth
with 2^20 trials (it should take about one hour):

	python3 smooth_proba.py 256 32 65 20
