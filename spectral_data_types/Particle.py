from .Channel import Channel
import pandas as pd
import numpy as np
import warnings 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import  make_interp_spline
import dtw
from lmfit.models import GaussianModel


class Particle: 
    channels: dict
    dtw_warped_counts: pd.Series
    peak_indices: list[pd.Series] # Will be calculated using channel allignment
    summed_heights: pd.Series # Will be calculated with sum channels
    summed_channel: Channel
    midpoints: pd.Series
    channels_s2n: list
    s2n : float
    
    
    # initalising by directly adding a dataframe with energies
    def __init__(self, spectral_dataframe: pd.DataFrame):
        """Initialise a particle object with a pd.DataFrame.
        Params:
            spectra_dataframe (pd.DataFrame) : A dataframe with loc 'channel name' and one column with the energies
        """
        self.channels = {}
        self.peak_indices = []
        self.summed_heights = pd.Series()
        self.dtw_warped_counts = pd.Series()

        for name, energies in spectral_dataframe.iterrows():
            self.channels[name] = Channel(energies)
        self.clean_channels()

    # for smooth and fast iteration over channels in the particle database
    def __iter__(self):
        return iter(self.channels)
    
    # returns how many channels particle class has
    def __len__(self):
        return len(self.channels)

    def __getitem__(self, channel_key: str):
        assert channel_key in self.channels.keys(), "Channel name not found in object's channels"
        return self.channels[channel_key]
    
    def apply(self, apply_func, inplace=False):
        """
        Apply a channel function to every channel within our particle.
        Params: 
            apply_func (func) : channel function to apply to every channel
            inplace (bool) : mutate the current list of channels or return the newly modified channels.
        Return:
            bool or None
        """
            
        try:
            new_channels = {}
            for key, channel in self.channels.items():
                new_channels[key] = apply_func(channel)
        except AttributeError:
            warnings.warn("Not a native channel method, attempting to apply anyways.")
            new_channels = {}
            for key, channel in self.channels.item():
                new_channels[key] = apply_func(channel)
        
        if inplace:
            # self.channels = new_channels
            return None
        else:
            return new_channels
                

                
    def sum_channels(self):
        '''Purpose to add all channels in particle
        ---Input: Self
        ---Output: np.array summed spectral values'''
        
        bins = np.arange(0, 18000, 2)
        summed_counts = np.zeros(len(bins) - 1)
        
        for chan in self.channels.values():
            hist, edges = np.histogram(chan.splined_midpoints, bins=bins, weights=chan.counts)
            
            summed_counts += hist
            
        self.summed_heights = summed_counts.astype(int)

        self.midpoints = bins[:-1] + 1  # midpoints of bins
        self.summed_channel = Channel(pd.Series(np.repeat(self.midpoints, self.summed_heights)))

        
        
    def dtw_warp(self, warp_on: pd.Series or str or None):
        """
        Use dynamic time warping to max channels based upon the master channel

        :param self:
        :param warp_on: pd.Series or str: Warp everything according to this channel

        :return: None
        """
        # Determining what to base the warp on.
        if type(warp_on) == str:
            master_channel = self.channels[warp_on]
        elif type(warp_on) == pd.Series:
            master_channel = warp_on
        else:
            master_channel = list(self.channels.values())[0] # NOTE: A big assumption is made here: That the first channel passed in here is the "best"
        

        warped_channel_counts = {}
        for name, chan in self.channels.items():
            allign = dtw.dtw(x=master_channel.counts, y=chan.counts)
            print(f"Done alligning chan: {name}.")

            warped = chan.counts[allign.index2]
            warped_channel_counts[name] = warped
            print(f"Saved chan: {name} to new")


        self.dtw_warped_counts = warped_channel_counts
        return
                
        
        
    def fit_peak_s2n(self, x_fit, y_fit):
        model = GaussianModel(prefix='g_')

        params = model.make_params(
            g_amplitude=y_fit.max(),
            g_center=x_fit[np.argmax(y_fit)],
            g_sigma=(x_fit[1] - x_fit[0]) * 10
        )

        result = model.fit(y_fit, params, x=x_fit)
        s2n = result.params['g_amplitude'].value / (2 * result.params['g_sigma'].value)
        return s2n

    

    def find_channels_s2n(self):
        '''
        Docstring for find_channels_s2n
        
        :param self:Calculate the height to width ratio of the peaks witin the channels.
        '''
        channel_s2n = []
        for channel in self.channels.values():
            x = channel.midpoints
            signal = channel.counts
            peak_s2n_list = []
            for peak_index in channel.prominent_peak_indices:
                fit_region_start = peak_index - 50
                fit_region_end = peak_index + 50

                x_fit = x[fit_region_start:fit_region_end]
                y_fit = signal[fit_region_start:fit_region_end]

                
                peak_s2n_list.append(self.fit_peak_s2n(x_fit, y_fit))
            avg_peak_s2n = np.mean(peak_s2n_list)
            channel_s2n.append(avg_peak_s2n)
        self.channel_s2n = channel_s2n

        return channel_s2n




    def find_particle_s2n(self):
        '''
        Docstring for find_particle_s2n
        
        :param self: find the height 2 width ratio of the peaks in the particle.
        '''

        self.summed_channel.scipy_peaks()
        signal = self.summed_heights
        
        samp_chan = self.channels['chan1']
        x = self.midpoints
        peak_s2n_list = []
        
        for peak_index in self.summed_channel.prominent_peak_indices:
            try:
                fit_region_start = peak_index - 50
                fit_region_end = peak_index + 50

                x_fit = x[fit_region_start:fit_region_end]
                y_fit = signal[fit_region_start:fit_region_end]
                
                peak_s2n_list.append(self.fit_peak_s2n(x_fit, y_fit))
            except ValueError:
                continue

        particle_s2n = np.mean(peak_s2n_list)
        self.s2n = particle_s2n
        return particle_s2n
        

    def find_s2n_zscore(self):
        '''
        Docstring for find_s2n_zscore
        
        :param self: calculate how close we are in aligning to any given channel
        '''

        try:
            mean = np.mean(self.channel_s2n)
            std = np.std(self.channel_s2n)
            particle_z = np.abs(self.s2n - mean) / std
            return particle_z
        except AttributeError:
            self.find_channels_s2n()
            self.find_particle_s2n()
            mean = np.mean(self.channel_s2n)
            std = np.std(self.channel_s2n)
            particle_z = np.abs(self.s2n - mean) / std
            return particle_z
        

    def algorithmic_aligning(self):
        '''
        Docstring for algorithmic_aligning
        
        :param self: Spline all peaks to channel1
        '''

        #find the spline points for chan1.
        master_channel = self.channels['chan1']
        master_peaks = master_channel.prominent_peak_indices
        other_master_peaks = np.delete(master_peaks, [0,1,-1])

        master_peak_heights = master_channel.flat_counts[other_master_peaks]
        master_tallest_peak = np.argmax(master_peak_heights)

        master_alignment_peaks = master_peaks[[0, 1, master_tallest_peak, -1]] #First, second, tallest, and last peak
        master_alignment_midpoints = master_channel.midpoints[master_alignment_peaks]
        master_channel.splined_midpoints = master_channel.midpoints


        #loop through all channels
        for channel_name, chan in self.channels.items():
            if channel_name == 'chan1':
                continue

            #find the spline points for chan
            chan_peaks = chan.prominent_peak_indices
            other_chan_peaks = np.delete(chan_peaks, [0, 1, -1])

            
            
            chan_peak_heights = chan.flat_counts[other_chan_peaks]
            chan_tallest_peak = np.argmax(chan_peak_heights)
            chan_alignment_peaks = chan_peaks[[0, 1, chan_tallest_peak, -1]]
            chan_alignment_midpoints = chan.midpoints[chan_alignment_peaks]

            #spline chan to channel1
            spl = make_interp_spline(chan_alignment_midpoints, master_alignment_midpoints, k = 3)
            chan_midpoints_t = spl(chan.midpoints)

            chan.splined_midpoints = chan_midpoints_t



    
    def clean_channels(self):
        '''
        Docstring for clean_channels
        
        :param self: Get rid of all bad channels
        '''
        
        channel_lengths = []
        channel_max_bin = []
        for channel_name in self.channels.keys():
            chan_x = self.channels[channel_name]
            channel_lengths.append(len(chan_x.energy))
            channel_max_bin.append(max(chan_x.midpoints))
        mean_length = np.mean(channel_lengths)
        std_length = np.std(channel_lengths)
        mean_bin = np.mean(channel_max_bin)
        std_bin = np.std(channel_max_bin)

        cutoff_length = mean_length - (2 * std_length)
        cutoff_bin = mean_bin - (2 * std_bin)


        to_remove = []
        for channel_name in self.channels.keys():
            chan_x = self.channels[channel_name]
            if len(chan_x.energy) < cutoff_length:
                to_remove.append(channel_name)
            elif max(chan_x.midpoints) < cutoff_bin:
                to_remove.append(channel_name)

        for channel_name in to_remove:
            self.channels.pop(channel_name)           
    


#-----------------Plotting Functions-----------------#

    def plot_summed_channels(self):
        '''
        Docstring for plot_summed_channels
        
        :param self: plot the summed spectral channels
        '''
        plt.plot(self.midpoints, self.summed_heights, label='Summed Channels')
        plt.xlabel('Energy')
        plt.ylabel('Counts')
        plt.title('Summed Spectral Channels')
        plt.yscale('log')
        plt.legend()
        plt.show()


    def waterfall_plot_particle(self):
        plt.figure(figsize=(10, 8))
        offset = 0
        for channel_name, chan in self.channels.items():
            plt.plot(chan.splined_midpoints, chan.counts + offset, label = False)
            if channel_name == 'chan109':
                break
            offset -= 800
        # plt.legend()
        plt.title('Aligned Spectral Channels')
        plt.xlabel('Energy')
        plt.yticks([])
        
        plt.show()



    def plot_channel(self, flat: bool = False, with_peaks: bool = False,  name: str = None):
        '''
        Docstring for plot_channel
        
        :param self: plot a specific channel by name
        '''
        chan = self.channels[name]
        chan.plot_channel(flat=flat, with_peaks=with_peaks, name=name)