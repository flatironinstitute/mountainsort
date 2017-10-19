import numpy as np

processor_name='pyms.extract_geom'
processor_version='0.1'
def extract_geom(*,geom,channels_array='',geom_out,channels=''):
    """
    Extract a subset of channels from a geom csv file

    Parameters
    ----------
    geom : INPUT
        Path of geom csv file having M rows and R columns where M is number of channels and R is the number of physical dimensions (1, 2, or 3)
    channels_array : INPUT 
        (optional) Path of array of channel numbers (positive integers). Either use this or the channels parameter, not both.
        
    geom_out : OUTPUT
        Path of output geom csv file containing a subset of the channels
        
    channels : string
        (Optional) Comma-separated list of channels to extract. Either use this or the channels_array input, not both.
    """    
    if channels:
        Channels=np.fromstring(channels,dtype=int,sep=',')
    elif channels_array:
        Channels=channels_array
    else:
        Channels.np.empty(0)        
    
    X=np.loadtxt(open(geom, "rb"), delimiter=",")
    X=X[(Channels-1).tolist(),:]
    np.savetxt(geom_out,X,delimiter=",",fmt="%g")
    return True
extract_geom.name=processor_name
extract_geom.version=processor_version
def test_extract_geom():
    G=np.array([[1,1],[2,1],[1,2],[2,2]])
    np.savetxt('tmp.geom.csv',G,delimiter=',',fmt='%g')
    extract_geom(geom='tmp.geom.csv',geom_out='tmp.geom2.csv',channels='1,2,4')
    G2=np.loadtxt(open('tmp.geom2.csv','rb'),delimiter=',')
    assert(G2.shape[0]==3)
    assert(G2.shape[1]==2)
    assert(np.array_equal(G[[0,1,3],:],G2))
    return True
extract_geom.test=test_extract_geom

if __name__ == '__main__':
    print ('Running test')
    test_extract_geom()