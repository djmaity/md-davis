import pandas
import numpy

def rmsf_to_df(xvg_file):
    """ Convert . xvg file containing RMSF to Pandas Data Frame """
    data = numpy.loadtxt(xvg_file, comments=('#','@'))
    df = pandas.DataFrame(data, columns=('Residue', 'RMSF'))
    df['Residue'] = df['Residue'].astype(int)
    return df