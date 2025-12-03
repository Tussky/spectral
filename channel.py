import pandas as pd
import numpy as np 
from scipy.signal import find_peaks

class channel:
    energy: pd.Series
    edges: pd.Series
    midpoints: pd.Series
    counts: pd.Series
    prominent_peak_indices: pd.Series    


    def __init__(self, energies: pd.Series):
       assert len(energies.shape) == 1, "energy needs to be a one dimensional pd.Series"

       self.energy = energies[(0 < energies) & (energies < np.percentile(energies, 97))] 
       self.edges, self.counts = np.histogram(self.energy, 10_000) 
       self.midpoints = 0.5 * (self.edges[1:] + self.edges[:1])

    def scipy_peaks(self, inplace: bool = False, scipy_prominence = 4, how_many_peaks = 10):
        """
        finds channel's peaks using scipy.find_peaks()
        Params:
            inplace (bool): if True, will save peak_indices and prominences to the channel object

        Returns:
            None or pd.Series
        """
        if inplace:
            assert self.prominent_peak_indices == None, "self.prominent_peak_indices is already defined, cannot overwrite - true inplace = False"
        
        peak_indices, peak_data = find_peaks(self.counts, prominence=scipy_prominence)
        prominences = peak_data['prominences']

        prominent_peak_indices = peak_indices[prominences.argsort()[-8:]]

        if inplace:
            self.prominent_peak_indices
            return None 
        return prominent_peak_indices
