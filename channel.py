import pandas as pd
import numpy as np 
from scipy.signal import find_peaks
from scipy.interpolate import  make_interp_spline
from scipy.signal import savgol_filter

class channel:
    energy: pd.Series
    edges: pd.Series
    midpoints: pd.Series
    counts: pd.Series
    prominent_peak_indices: pd.Series  
    flat_counts: pd.Series


    def __init__(self, energies: pd.Series):
       assert len(energies.shape) == 1, "energy needs to be a one dimensional pd.Series"
       self.energy = energies[(0 < energies) & (energies < np.percentile(energies, 99.5))]
       self.counts, self.edges = np.histogram(self.energy, 10_000) 
       self.midpoints = 0.5 * (self.edges[1:] + self.edges[:-1])
       

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

    def deriv_peaks(self):
        d2 = savgol_filter(self.energy, 10, 3)
        prominent_peak_indices = peak_indices[prominences.argsort()[-8:]]


    def spline_baseline_subtract(self):
        '''Purpose: Subtract the baseline of the channel using low-binned spline
        ---Input: channel self
        ---Output: pd.Series of flattened counts'''
        
        #---Create low-binned dataframe for spline   
        counts_low, bin_edges_low = np.histogram(self.energy, bins = 50) #create low-binned histogram
        bin_edges_low = bin_edges_low[:-1]  #drop last bin
        counts_low = counts_low / 100 #adjust for difference in binning
    
        #---Create Spline to fit low-binned dataframe
        spl = make_interp_spline(bin_edges_low, counts_low, k=3) #create spline object
        x = np.linspace(0, max(self.energy), len(self.edges))
        spline = spl(x) #make spline line
        spline = np.maximum(spline, 0)

        
        #---Remove the baseline of data using spline
        counts = self.counts.copy()
        flat_counts= counts - spline
        flat_counts = pd.Series(flat_counts)
        self.flat_counts = flat_counts
        return flat_counts

    def savgol_baseline_subtract(self):
        '''Purpose: Subtract baseline of channel using savgol_filter
        ---Input: channel self
        ---Output: pd.Series of flattened counts'''
        
        counts = self.counts.copy()
        baseline = savgol_filter(counts, 200, 3)
        flat_counts = counts - baseline
        flat_counts = np.maximum(flat_counts, 0)
        flat_counts = pd.Series(flat_counts)
        self.flat_counts = flat_counts
        return flat_counts

    
    def plot_channel(self, flat: bool = False, with_peaks: bool = False):
        '''Purpose: Plot raw channel
        ---Input: Self
        ---Output: None'''
        
        if(flat):
            df = pd.DataFrame({
                'Count': self.flat_counts,
                'Bin' : self.midpoints
                })
            ax = df.plot(drawstyle = 'steps-pre', x = 'Bin', y =  'Count', logy = False, color = 'orange')
            if(with_peaks):
                peaks = self.peak_indicies
                ax.scatter(self.edges[peaks], self.flat_counts[peaks], color = 'r', marker = 'o')
        else:
            df = pd.DataFrame({
                'Count': self.counts,
                'Bin' : self.midpoints
                })
            ax = df.plot(drawstyle = 'steps-pre', x = 'Bin', y =  'Count', logy = True, color = 'orange')
            if(with_peaks):
                peaks = self.peak_indicies
                ax.scatter(self.edges[peaks], self.counts[peaks], color = 'r', marker = 'o')
        

        
        
        