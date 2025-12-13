import pandas as pd
import numpy as np 
from scipy.signal import find_peaks
from scipy.interpolate import  make_interp_spline
from scipy.signal import savgol_filter

class Channel:
    name: str
    energy: pd.Series
    edges: pd.Series
    midpoints: pd.Series
    norm_midpoints: pd.Series
    splined_midpoints: pd.Series
    alligned_counts: pd.Series
    counts: pd.Series
    prominent_peak_indices: pd.Series
    prominent_peak_edges: pd.Series  
    flat_counts: pd.Series
    WANTED_PROM_PEAKS: int
    


    def __init__(self, energies: pd.Series, wanted_prom_peaks: int = 8):
       assert len(energies.shape) == 1, "energy needs to be a one dimensional pd.Series"
       energies = energies['energy'] # to stop the a.any() error appearing, despite this being a 1D Series.
       self.energy = energies[(0 < energies) & (energies < np.percentile(energies, 99.5))]
       self.counts, self.edges = np.histogram(self.energy, 10_000) 
       self.midpoints = 0.5 * (self.edges[1:] + self.edges[:-1])
       self.WANTED_PROM_PEAKS = wanted_prom_peaks
       self.prominent_peak_indices = pd.Series()
       
    def __str__(self):
        ret_str = "class: Channel"
        ret_str = ret_str + "energy: " + str(self.energy) + "\n"
        ret_str = ret_str + "edges: " + str(self.edges) + "\n"
        ret_str = ret_str + "counts: " + str(self.counts) + "\n"
        ret_str = ret_str + "wanted_peaks: " + str(self.WANTED_PROM_PEAKS) + "\n"
        return ret_str


    def scipy_peaks(self, inplace: bool = True, scipy_prominence = 4, how_many_peaks = 10):
        """
        Finds channel's peaks using scipy.find_peaks() and returns edges as well as peaks
        Params:
            inplace (bool): if True, will save peak_indices and prominences to the channel object

        Returns:
            None or tuple(pd.Series,pd.Series)
        """
        
        peak_indices, peak_data = find_peaks(self.counts, prominence=scipy_prominence)
        prominences = peak_data['prominences']

        prominent_peak_indices = peak_indices[prominences.argsort()[-self.WANTED_PROM_PEAKS:]]
        prominent_peak_left = peak_data['left_bases'][prominences.argsort()[-self.WANTED_PROM_PEAKS:]]
        prominent_peak_right = peak_data['right_bases'][prominences.argsort()[-self.WANTED_PROM_PEAKS:]]
        # print("LEFTIES")
        # print(prominent_peak_left)
        # print("RIGHTIES")
        # print(prominent_peak_right)

        if inplace:
            self.prominent_peak_indices = prominent_peak_indices
            return None 
        return (peak_indices, prominent_peak_indices, peak_data)

    def deriv_peaks(self):
        d2 = savgol_filter(self.counts, 10, 3)
        
        prominent_peak_indices = d2.argsort()[:self.WANTED_PROM_PEAKS]
        self.prominent_peak_indices = prominent_peak_indices
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
        x = np.linspace(0, max(self.energy), len(self.edges) - 1)
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

    
    def plot_channel(self, flat: bool = False, with_peaks: bool = False, name: str = None):
        '''Purpose: Plot raw channel
        ---Input: Self
        ---Output: None'''
        if(flat):
            df = pd.DataFrame({
                'Count': self.flat_counts,
                'Energy' : self.midpoints
                })
            ax = df.plot(drawstyle = 'steps-pre', x = 'Energy', y =  'Count', logy = False, color = 'orange', label = name, legend= True)
            if(with_peaks):
                peaks = self.prominent_peak_indices
                ax.scatter(self.edges[peaks], self.flat_counts[peaks], color = 'r', marker = 'o')
        else:
            df = pd.DataFrame({
                'Count': self.counts,
                'Energy' : self.midpoints
                })
            ax = df.plot(drawstyle = 'steps-pre', x = 'Energy', y =  'Count', logy = True, color = 'orange', label = name, legend= True)
            if(with_peaks):
                peaks = self.prominent_peak_indices
                ax.scatter(self.edges[peaks], self.counts[peaks], color = 'r', marker = 'o')

        
        
    def normalize_midpoints(self):
        scaler = 20_000 / self.midpoints[-1]

        self.norm_midpoints = self.midpoints * scaler

        return self.norm_midpoints
        
        