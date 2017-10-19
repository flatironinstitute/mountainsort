import numpy as np

import sys,os
parent_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path)

import mlpy

processor_name='pyms.concatenate_firings'
processor_version='0.1'
def concatenate_firings(*,firings_list, firings_out, time_offsets, increment_labels='false'):
    """
    Combine a list of firings files to form a single firings file

    Parameters
    ----------
    firings_list : INPUT
        A list of paths of firings mda files to be concatenated
    firings_out : OUTPUT
        ...

    time_offsets : string
        An array of time offsets for each firings file. Expect one offset for each firings file.
        ...
    increment_labels : string
        ...
    """
    if time_offsets:
        time_offsets=np.fromstring(time_offsets,dtype=np.float_,sep=',')
    else:
        time_offsets=np.zeros(len(firings_list))
    if len(firings_list) == len(time_offsets):
        concatenated_firings=np.array([[],[],[],[]])
        for idx, firings in enumerate(firings_list):
            to_append=mlpy.readmda(firings)
            if increment_labels=='true':
                to_append[1,:]+=time_offsets[idx]
            concatenated_firings = np.append(concatenated_firings,to_append, axis=1)
        mlpy.writemda64(concatenated_firings,firings_out)
        return True
    else:
        print('Mismatch between number of firings files and number of offsets')
        return False

def test_concatenate_firings():
    M,N1,N2=4,2000,30000
    test_offset_str='300000,123456789'
    test_offset=[300000,123456789]
    fir1=np.around(np.random.rand(M,N1),decimals=3)
    mlpy.writemda64(fir1,'tmp.fir1.mda')
    fir2=np.around(np.random.rand(M,N2),decimals=3)
    mlpy.writemda64(fir2,'tmp.fir2.mda')
    fir1_incr=fir1
    fir2_incr=fir2
    fir12 = np.append(fir1, fir2, axis=1)
    fir12=np.around(fir12,decimals=3)
    fir1_incr[1,:]+=test_offset[0]
    fir2_incr[1,:]+=test_offset[1]
    fir12_incr = np.append(fir1_incr,fir2_incr,axis=1)
    concatenate_firings(firings_list=['tmp.fir1.mda','tmp.fir2.mda'],firings_out='tmp.test_fir12.mda',time_offsets=test_offset_str,increment_labels='false')
    concatenate_firings(firings_list=['tmp.fir1.mda', 'tmp.fir2.mda'],firings_out='tmp.test_fir12_incr.mda', time_offsets=test_offset_str, increment_labels='true')
    test_fir12=mlpy.readmda('tmp.test_fir12.mda')
    test_fir12=np.around(test_fir12,decimals=3)
    test_fir12_incr=mlpy.readmda('tmp.test_fir12_incr.mda')
    test_fir12_incr=np.around(test_fir12_incr,decimals=3)
    np.testing.assert_array_almost_equal(fir12,test_fir12,decimal=3)
    np.testing.assert_array_almost_equal(fir12_incr,test_fir12_incr,decimal=3)
    return True

concatenate_firings.test=test_concatenate_firings
concatenate_firings.name = processor_name
concatenate_firings.version = processor_version
concatenate_firings.author = 'J Chung and J Magland'

if __name__ == '__main__':
    print ('Running test')
    test_concatenate_firings()