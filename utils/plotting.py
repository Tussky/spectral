import pandas as pd
import numpy as np 
import plotly.express as px 
import plotly.graph_objects as go 

def plot_spectra_on_histogram(histogram: pd.Series, midpoints: pd.Series, peak_indices: pd.Series, should_show = True):
    counts_hist = px.line(
        x = midpoints,
        y = histogram,
        labels = {"x":"energy_level", "y":"frequency"}
    )

    highest_peaks = px.scatter(
        x = midpoints[peak_indices],
        y = histogram[peak_indices],
        color_discrete_sequence = ['purple']
    )
    highest_peaks.update_traces(marker=dict(
        size = 8, 
    ))
 
    counts_hist.add_traces(highest_peaks.data)
    return counts_hist