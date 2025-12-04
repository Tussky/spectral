import pandas as pd
import numpy as np 
import h5py

def read_in_hdf5_spectra(path = "./../Data/Gamma/210601_NBS295-106/20210601_152616_mass-001.hdf5"):
    """
    Args:
        path (str): The path to your file. Defaults.
    Returns: 
        pd.DataFrame: A DataFrame of series for all included channels.
    
    """
    filename = "./../Data/Gamma/210601_NBS295-106/20210601_152616_mass-001.hdf5"
    with h5py.File(filename, "r") as hdf_file:
        channels = pd.DataFrame(
            columns=["energy"],
            index=hdf_file.keys(),
        )

        for channel_name in hdf_file:
            going_in = np.array(hdf_file[channel_name]["filt_value"])
            going_in = going_in[(0 < going_in) & (going_in < np.percentile(going_in, 97))]
            channels.loc[channel_name] = [going_in]
    return channels