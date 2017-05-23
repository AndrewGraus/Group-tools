####
#Just some simple funtions to make modules out of

def high_low_limit(data, confidence):
    #############################################
    #
    # Input: data - a list of data for which the 
    #        confidence limits will be computed
    #        confidence - the confidence interval
    #        that will be computed so if you
    #        want the 95% confidence interval
    #        that's 0.95
    # 
    # Output: the average, and the upper and lower
    #         limits of the confidence interval
    #
    ##############################################
    n_low = 0.5-confidence/2.0
    n_high = 0.5+confidence/2.0
    n_avg = 0.5
    sorted_data = sorted(data)
    lim_low = sorted_data[int(len(data)*n_low)]
    lim_high = sorted_data[int(len(data)*n_high)]
    lim_avg = sorted_data[int(len(data)*n_avg)]
    #print len(data)*n_high, int(len(data)*n_high)
    return lim_low, lim_high, lim_avg

def anticumulative_hist(data,bins=None):
    """
    Makes an anticumulative histogram (e.g. N(>Vmax)) when fed data (e.g. Vmax
    values).
    
    :param data:
    the input array to be made into a histogram, already presliced etc.
    
    :param bins:
    if given, used as the bins for the histogram; if not given, the unique
    values in the data are used as the histogram (so the returned histogram
    is exact in that case)
    
    :returns:
    the histogram, then the bins, so to plot, do, e.g.
    hist,bins = anticum_hist(vmax[slice])
    loglog(bins,hist)
    """

    from numpy import histogram,unique,append,cumsum,zeros
    
    if bins == None:  
        temp,bins = histogram(data,unique(data))
        hist = append(cumsum(temp[::-1])[::-1],1)
    else:
        if len(data) == 0:
            return zeros(len(bins)),bins
        if max(data) > max(bins):
            print "Must have the largest bin be bigger than the largest data point."
            print max(bins)
            print max(data)
            import sys
            sys.exit(1337)
        temp,bins = histogram(data,bins)
        hist = append(cumsum(temp[::-1])[::-1],0)
        
    return hist,bins

def get_distance(position_array, center):
    #simple script that calculates a list of distances given an array of positions and a 'center'
    import numpy as np
    x, y, z = position_array[:,0], position_array[:,1], position_array[:,2]
    dx, dy, dz = x-center[0], y-center[1], z-center[2]
    return np.sqrt(dx*dx+dy*dy+dz*dz)

def get_distance_vector(position_array, center):
    #simple script that calculates a list of distances given an array of positions and a 'center'
    import numpy as np
    x, y, z = position_array[:,0], position_array[:,1], position_array[:,2]
    dx, dy, dz = x-center[0], y-center[1], z-center[2]
    dist_array = np.zeros((len(dz),3))
    dist_array[:,0], dist_array[:,1], dist_array[:,2] = dx, dy, dz
    return dist_array
    
