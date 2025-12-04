import pandas as pd
import channel
class particle: 
    channels: list[channels]
    peak_indices: list[pd.Series]

    for channel in channels:
        self.peak_indices.append(channel.peak_indices)