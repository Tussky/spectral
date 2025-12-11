from .Channel import Channel
import pandas as pd
import numpy as np
import warnings 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class Particle: 
    channels: dict
    peak_indices: list[pd.Series] # Will be calculated using channel allignment
    summed_heights: pd.Series # Will be calculated with sum channels
    summed_channel: Channel
    channel_h2w: list
    h2w : float
    
    # intialising by adding channels
    def __init__(self, channels_to_add: Channel):
        assert type(channels_to_add)==list and type(channels_to_add[0]) == Channel, "channels_to_add must be a list of channels"
        self.channels = channels_to_add


    # initalising by directly adding a dataframe with energies
    def __init__(self, spectral_dataframe: pd.DataFrame):
        """Initialise a particle object with a pd.DataFrame.
        Params:
            spectra_dataframe (pd.DataFrame) : A dataframe with loc 'channel name' and one column with the energies
        """
        self.channels = {}
        self.peak_indices = []
        self.summed_heights = pd.Series()

        for name, energies in spectral_dataframe.iterrows():
            self.channels[name] = Channel(energies)

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
            self.channels = new_channels
            return None
        else:
            return new_channels
                

                
    def sum_channels(self):
        '''Purpose to add all channels in particle
        ---Input: Self
        ---Output: np.array summed spectral values'''
        
        summed = []
        nbins = len(self.channels[0].midpoints)
        for i in range(nbins):
            bin_total = sum(ch.counts[i] for ch in self.channels)
            summed.append(bin_total)
        self.summed_heights = pd.Series(summed)
        

    def allign_peaks(self):
        
        return
        
        
    def plot_particle(self):
        chan = self.channels[0]
        bins = chan.midpoints
        plt.plot(bins, self.summed_heights, drawstyle='steps-mid')

    def gaussian(x, amplitude, mean, stdev):
        return amplitude * np.exp(-((x- mean) / (2* stdev)) **2)

    def find_channels_h2w(self):
        '''
        Docstring for find_channels_h2w
        
        :param self:Calculate the height to width ratio of the peaks witin the channels.
        '''
        channel_h2w = []
        for i, channel in enumerate(self.channels):
            x = channel.midpoints
            signal = channel.counts
            peak_h2w_list = []
            for peak_index in self.peak_indices[i]:
                fit_region_start = peak_index - 50
                fit_region_end = peak_index + 50

                x_fit = x[fit_region_start:fit_region_end]
                y_fit = signal[fit_region_start:fit_region_end]

                init_guess = [signal[peak_index], x[peak_index], 5]

                params, covariance = curve_fit(gaussian, x_fit, y_fit, p0 = init_guess)

                amplitude_fit, mean_fit, stddev_fit = params

                peak_h2w = amplitude_fit / (2 * stddev_fit)
                peak_h2w_list.append(peak_h2w)
            avg_peak_h2w = np.mean(peak_h2w_list)
            channel_h2w.append(avg_peak_h2w)
        self.channel_h2w = channel_h2w

        return channel_h2w



    def find_particle_h2w(self):
        '''
        Docstring for find_particle_h2w
        
        :param self: find the height 2 width ratio of the peaks in the particle.
        '''

        particle_peaks = self.summed_channel.scipy_peaks
        particle_peak_indices = particle_peaks[0]
        signal = self.summed_heights
        samp_chan = self.channels['chan01']
        x = samp_chan.midpoints
        peak_h2w_list = []
        for peak_index in particle_peak_indices:
            fit_region_start = peak_index - 50
            fit_region_end = peak_index + 50

            x_fit = x[fit_region_start:fit_region_end]
            y_fit = signal[fit_region_start:fit_region_end]

            init_guess = [signal[peak_index], x[peak_index], 5]

            params, covariance = curve_fit(gaussian, x_fit, y_fit, p0 = init_guess)

            amplitude_fit, mean_fit, stddev_fit = params

            peak_h2w = amplitude_fit / (2 * stddev_fit)
            peak_h2w_list.append(peak_h2w)

        particle_h2w = np.mean(peak_h2w_list)
        self.h2w = particle_h2w
        return particle_h2w
        
    def find_h2w_zscore(self):
        '''
        Docstring for find_h2w_zscore
        
        :param self: calculate how close we are in aligning to any given channel
        '''

        try:
            mean = np.mean(self.channel_h2w)
            std = np.std(self.channel_h2w)
            particle_z = np.abs(self.h2w - mean) / std
            return particle_z
        except AttributeError:
            self.find_channels_h2w
            self.find_particle_h2w
            mean = np.mean(self.channel_h2w)
            std = np.std(self.channel_h2w)
            particle_z = np.abs(self.h2w - mean) / std
            return particle_z
        



                



    
                