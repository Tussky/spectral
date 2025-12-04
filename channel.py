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
       assert len(energies.shape) == 1, "energies needs to be a one dimensional pd.Series"
       assert "energy" in energies.keys(), "energies series must have a key called 'energy'"
       self.energy = energies['energy'][(0 < energies['energy']) & (energies['energy'] < np.percentile(energies['energy'], 99))] 
       self.counts, self.edges = np.histogram(self.energy, 10_000) 
       self.midpoints = 0.5 * (self.edges[1:] + self.edges[:-1])

    ####################
    ### PEAK FINDING ###
    ####################
    def scipy_peaks(self, inplace: bool = False, scipy_prominence = 4, how_many_peaks = 10):
        """
        Finds channel's peaks using scipy.find_peaks()
        Params:
            inplace (bool): if True, will save peak_indices and prominences to the channel object

        Returns:
            None or pd.Series
        """
        if inplace:
            assert self.prominent_peak_indices == None, "self.prominent_peak_indices is already defined, cannot overwrite - try inplace = False"
        
        peak_indices, peak_data = find_peaks(self.counts, prominence=scipy_prominence)
        prominences = peak_data['prominences']

        prominent_peak_indices = peak_indices[prominences.argsort()[-8:]]

        if inplace:
            self.prominent_peak_indices
            return None 
        return prominent_peak_indices


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
        
        return flat_counts

    def savgol_baseline_subtract(self):
        '''Purpose: Subtract baseline of channel using savgol_filter
        ---Input: channel self
        ---Output: pd.Series of flattened counts'''
        
        counts = self.counts.copy()
        baseline = savgol_filter(counts, 200, 3)
        flat_counts = counts - baseline
        flat_counts = pd.Series(flat_counts)

        return flat_counts

    def plot_raw_channel(self):
        '''Purpose: Plot raw channel
        ---Input: Self
        ---Output: None'''

        df = get_channel_df(channel)
        ax = df.plot(drawstyle = 'steps-pre', x = 'Bin_Left_Edge', y =  'Count', logy = True, color = 'orange', xlim = (0, 20000))
        bin_edges = df['Bin_Left_Edge']
        heights = df['Count']
        peaks, peak_heights = get_channel_peaks(channel) 
        ax.scatter(bin_edges[peaks], heights[peaks], color = 'r', marker = 'o')
        plt.title('og data with peaks')
        