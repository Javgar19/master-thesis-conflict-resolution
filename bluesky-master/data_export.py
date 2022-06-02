from dataclasses import dataclass
from unicodedata import name
import pandas as pd
import numpy as np
import os 
import json

def load_data(rootdir = "./results"):
    """
    Read the metric.json files and create a dataframe with all the simulations data 
    """
    dataframes = []
    for subdir, directories, files in os.walk("./results"):
        for directory in directories:
            if directory.isnumeric():
                filename = os.path.join(subdir, directory, "metrics.json")
                with open(filename, "r") as f:
                    dataframes.append(json_to_df(json.load(f)))
                    print(f"Run {directory} loaded")

    print("All runs loaded. Concatenating data...")
    df = dataframes.pop(0)
    for dataframe in dataframes:
        df = pd.concat([df, dataframe], ignore_index=True)
                
    return df

def json_to_df(metrics):
    """
    Convert the json file once readed to a pandas dataframe with the relevant variables
    """
    columns = list(metrics.keys())
    data = []
    
    for key in columns:
        data.append(metrics[key]["values"])
    
    data = np.array(data)
    return pd.DataFrame(data=data.transpose(), columns=columns)


def main():
    data = load_data()
    print("Dataframe created. Exporting it as csv...")

    filename = './results/simulations.csv'
    data.to_csv(filename)
    print("Data succesfully exported")
    
if __name__ == '__main__':
    main()