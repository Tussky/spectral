import pandas as pd
import numpy as np 


class channel:
    energy: pd.Series
    edges: pd.Series
    midpoints: pd.Series
    counts: pd.Series
    peaks: pd.Series
    prominent_peaks: pd.Series

    def __init__(energy: pd.Series):
        assert energy.shape[0] == 1, "energy needs to be a one dimensional pd.Series"
        
       self.energy = energy[(0 < energy) && (energy < np.percentile(energy, 97))] 





