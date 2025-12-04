class baseline_subtraction:
    def get_channel_df(channel):
    ''' Purpose: create the dataframe for a channel
        Input: string channel
        output: dataframe df'''

    #---Get and clean channel
    data = channels_df.loc[channel]['data'] #Get channel data
    data = data[(data >=0) & (data < 24000)] #Clean edges

    # Make single data at boundaries so that bins are always same
    data = np.append(data, [0, 24000])

    
    #---create channel data frame with bins and bin_counts
    counts, bin_edges = np.histogram(data, bins = 6000) #create histogram
    bin_edges = bin_edges[:-1] #drop last bin edge so len(heights) == len(bin_edges)

    #---normalize data
    counts = counts / len(data)

    #create df
    df = pd.DataFrame({
        'Bin_Left_Edge': bin_edges,
        'Count': counts}) 

    return df

    
    def make_spline(channel):
        '''---Purpose: Create a spline line for baseline removal
        ------Input: Takes in the channel data as np.array, and bins as np.array
        ------Output: Returns the spline line as np.array---'''
        #---Get and clean channel
        data = channels_df.loc[channel]['data'] #Get channel data
        data = data[(data >=0) & (data < 24000)] #Clean edges
        
    
        #---Create low-binned dataframe for spline   
        counts_low, bin_edges_low = np.histogram(data, bins = 50) #create low-binned histogram
        bin_edges_low = bin_edges_low[:-1]  #drop last bin
        counts_low = counts_low / 100 #adjust for difference in binning
    
        #---Create Spline to fit low-binned dataframe
        spl = make_interp_spline(bin_edges_low, counts_low, k=3) #create spline object
        x = np.linspace(0, 24000, 6000)
        spline = spl(x) #make spline line
        spline = np.maximum(spline, 0)
    
        return spline


    
    
    def savgol():
        return

    def wavelets():
        return 

    def low_binning(channels_df):
        
            
        return 
