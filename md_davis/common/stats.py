def moving_avg(data, window=1, shift=1):
    """ Canculate the mean and standard deviation of equal length
        non-overlapping contiguous subsets of the timeseries data """
    rows = int(numpy.ceil(len(data) / window))
    mean = numpy.empty((rows,))
    index = 0
    start, stop = 0, window
    while stop < len(data):
        mean[index] = numpy.mean(data[start:stop])
        start = start + shift
        stop = stop + shift
        index += 1
    mean[index] = numpy.mean(data[start:])
    return mean


def rolling_std(data, window=1, shift=1):
    """ Canculate the mean and standard deviation of equal length
        non-overlapping contiguous subsets of the timeseries data """
    rows = int(numpy.ceil(len(data) / window))
    std = numpy.empty((rows,))
    index = 0
    start, stop = 0, window
    while stop < len(data):
        std[index] = numpy.std(data[start:stop])
        start = start + shift
        stop = stop + shift
        index += 1
    std[index] = numpy.std(data[start:])
    return std


def moving_avg_std(data, window=1, shift=1):
    """ Canculate the mean and standard deviation of equal length
        non-overlapping contiguous subsets of the timeseries data """
    rows = int(numpy.ceil(len(data) / window))
    mean = numpy.empty((rows,))
    std = numpy.empty((rows,))
    index = 0
    start, stop = 0, window
    while stop < len(data):
        mean[index] = numpy.mean(data[start:stop])
        std[index] = numpy.std(data[start:stop])
        start = start + shift
        stop = stop + shift
        index += 1
    mean[index] = numpy.mean(data[start:])
    std[index] = numpy.std(data[start:])
    return mean, std