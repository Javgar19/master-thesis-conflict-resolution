from tabnanny import check
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def check_completitude(df):
    df = df.groupby(['radius','n_ac','threshold']).count()
    print(df)


def main():
    df = pd.read_csv('./results/simulations.csv')
    print(df.corr())



if __name__ == '__main__':
    main()
