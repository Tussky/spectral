import pandas as pd
import channel
class particle: 
    channels: list[channels]
    peak_indices: list[pd.Series]
    summed_channels: pd.Series
    for channel in channels:
        self.peak_indices.append(channel.peak_indices)

    def sum_channels(self):
        '''Purpose to add all channels in particle
        ---Input: Self
        ---Output: np.array summed spectral values'''
        
        summed = []
        nbins = len(channels[0].midpoints)
        for i in range(n_bins):
            bin_total = sum(ch.energy[i] for ch in self.channels)
            summed.append(bin_total)
        summed_channels = pd.Series(summed)
        
        
    def plot_particle(self):
        chan = self.channels[0]
        bins = chan.midpoints
        plt.plot(bins, self.summed_channels, drawstyle='steps-mid')
                