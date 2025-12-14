from .Channel import Channel
import pandas as pd
import numpy as np
import warnings 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import  make_interp_spline
from dtw import dtw 


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
        self.clean_channels()


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
        
        summed = []
        nbins = len(self.channels[0].midpoints)
        for i in range(nbins):
            bin_total = sum(ch.counts[i] for ch in self.channels)
            summed.append(bin_total)
        self.summed_heights = pd.Series(summed)
        

    def allign_peaks(self):
        
        return
    
    def dtw_allignment(self):
        """
        Docstring for dtw_allignment
        
        :param self: Alligns all the channels from the particle class to channel1 - does so in place.
        """
        the_master_channel = self.channels['chan1']
        print(self.channels)
        for name, iter_channel in self.channels.items():
            if name == 'chan1':
                continue
            else:    
                allignment = dtw(the_master_channel.counts, iter_channel.counts)
                iter_channel.alligned_counts = iter_channel.counts[allignment.index2] 
        print("Done with: ", name)
                
        
    def waterfall_plot_particle(self):
        # scaled_channels
        plt.figure(figsize=(10, 100))
        offset = 0

        for midpoints, counts, channel_name in zip(scaled_channels['bin_left_edge'],scaled_channels['counts'], scaled_channels['channel']):
            counts = np.array(counts)
            plt.plot(bins, counts + offset, label = channel_name)
            offset -= 2000
        plt.legend()
        plt.show()


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
        samp_chan = self.channels['chan1']
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
        

    def algorithmic_aligning(self):
        '''
        Docstring for algorithmic_aligning
        
        :param self: Spline all peaks to channel1
        '''

        #find the spline points for chan1.
        master_channel = self.channels['chan1']
        master_peaks = master_channel.prominent_peak_indices
        master_peak_heights = master_channel.flat_counts[master_peaks]
        master_tallest_peak = np.argmax(master_peak_heights)
        master_alignment_peaks = master_peaks[[0, 1, master_tallest_peak, -1]] #First, second, tallest, and last peak
        master_alignment_midpoints = master_channel.midpoints[master_alignment_peaks]


        #loop through all channels
        for channel_name, chan in self.channels.items():
            if channel_name == 'chan1':
                continue

            #find the spline points for chan
            chan.scipy_peaks()
            chan_peaks = chan.prominent_peak_indices
            # print(chan_peaks)
            other_chan_peaks = np.delete(chan_peaks, [0, 1, -1])
            # print(chan_peaks)
            
            chan.savgol_baseline_subtract()
            chan_peak_heights = chan.flat_counts[other_chan_peaks]
            chan_tallest_peak = np.argmax(chan_peak_heights)
            chan_alignment_peaks = chan_peaks[[0, 1, chan_tallest_peak, -1]]
            chan_alignment_midpoints = chan.midpoints[chan_alignment_peaks]


            #spline chan to channel1
            print(f'{channel_name}', chan_alignment_midpoints, '\n chan1:', master_alignment_midpoints)
            spl = make_interp_spline(chan_alignment_midpoints, master_alignment_midpoints, k = 3)
            chan_midpoints_t = spl(chan.midpoints)

            chan.spline_midpoints = chan_midpoints_t



    def old_algorithmic_aligning(self):
        #normalize all midpoints
        for channel in self.channels.values():
            channel.normalize_midpoints()

        #Get chan1 peak heights, energy
        aligning_channel = self.channels['chan1']
        chan1_peak_inds = aligning_channel.prominent_peak_indices
        chan1_peak_heights = aligning_channel.counts[chan1_peak_inds]
        chan1_peak_energies = aligning_channel.norm_midpoints[chan1_peak_inds]
        chan1_peaks = zip(chan1_peak_heights, chan1_peak_energies)

        #find std of height and energy
        std_chan1_energies = np.std(chan1_peak_energies)
        std_chan1_heights = np.std(chan1_peak_heights)


        all_aligning_peaks = {}
        for channel_name, channel in self.channels.items():
            chan_min_scores = []
            for peak_ind in channel.prominent_peak_indices:
                peak_height = channel.counts[peak_ind]
                peak_energy = channel.norm_midpoints[peak_ind]
                min_chan1_i, min_score = min([(i, (chan1_peak_heights[i] - peak_height)**2 + (chan1_peak_energies[i] - peak_energy)**2) for i in range(len(chan1_peak_inds))])
                chan1_energy = chan1_peak_energies[min_chan1_i]
                #TODO figure out math so can do 1 - h/h, difference should be normalized ont sheer number.
                chan_min_scores.append((peak_energy, chan1_energy, min_score))

            
            chan_aligning_peaks = [(peak_energy, chan1_energy) for peak_energy, chan1_energy, min_score in sorted(chan_min_scores, key=lambda x: x[2])[:5]]
            all_aligning_peaks[channel_name] = chan_aligning_peaks

        for channel in self.channels.keys():
            aligment_peaks = all_aligning_peaks[channel]
<<<<<<< HEAD
=======
            current_channel_midpoints, chan1_midpoins =  [(tuple[0], tuple[1]) for tuple in aligment_peaks]
            spl = make_interp_spline(current_channel_midpoints, chan1_midpoins, k = 3)
>>>>>>> 3215c77 (Fixed tuple error - added plotting notebook)

            
            current_channel_midpoints = [t[0] for t in aligment_peaks]
            chan1_midpoints = [t[1] for t in aligment_peaks]
            print(type(current_channel_midpoints), current_channel_midpoints)
            print(type(chan1_midpoints), chan1_midpoints)
            chan1_midpoints = sorted(chan1_midpoints)
            current_channel_midpoints = sorted(current_channel_midpoints)

            try:
                spl = make_interp_spline(current_channel_midpoints, chan1_midpoints, k = 3)
            except IndexError:
                print(channel)
            channel_x = self.channels[channel]
            channel_x_midpoints = channel_x.midpoints
            channel_x_t = spl(chan1_midpoints)

            channel_x.splined_midpoints = channel_x_t




    def clean_channels(self):
        '''
        Docstring for clean_channels
        
        :param self: Get rid of all bad channels
        '''
        print('called clean channels')
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
        





        


            


            
                






                    




                



    
                