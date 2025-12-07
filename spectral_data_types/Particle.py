from .Channel import Channel
import pandas as pd
import warnings 
class Particle: 
    channels: dict
    peak_indices: list[pd.Series] # Will be calculated using channel allignment
    summed_channels: pd.Series # Will be calculated with sum channels

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
        self.summed_channels = pd.Series()

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
        for i in range(n_bins):
            bin_total = sum(ch.energy[i] for ch in self.channels)
            summed.append(bin_total)
        self.summed_channels = pd.Series(summed)
        
        
    def plot_particle(self):
        chan = self.channels[0]
        bins = chan.midpoints
        plt.plot(bins, self.summed_channels, drawstyle='steps-mid')
                