import numpy as np
import struct

class MdaHeader:
    def __init__(self, dt0, dims0):
        uses64bitdims=(max(dims0)>2e9)            
        self.uses64bitdims=uses64bitdims
        self.dt_code=_dt_code_from_dt(dt0)
        self.dt=dt0
        self.num_bytes_per_entry=get_num_bytes_per_entry_from_dt(dt0)
        self.num_dims=len(dims0)
        self.dimprod=np.prod(dims0)
        self.dims=dims0
        if uses64bitdims:
            self.header_size=3*4+self.num_dims*8
        else:
            self.header_size=(3+self.num_dims)*4

class DiskReadMda:
    def __init__(self,path,header=None):
        self._path=path
        if header:
            self._header=header
            self._header.header_size=0
        else:
            self._header=_read_header(self._path)            
    def dims(self):
        return self._header.dims
    def N1(self):
         return int(self._header.dims[0])
    def N2(self):
        return int(self._header.dims[1])
    def N3(self):
        return int(self._header.dims[2])
    def dt(self):
        return self._header.dt
    def numBytesPerEntry(self):
        return self._header.num_bytes_per_entry
    def readChunk(self,i1=-1,i2=-1,i3=-1,N1=1,N2=1,N3=1):
        #print("Reading chunk {} {} {} {} {} {}".format(i1,i2,i3,N1,N2,N3))
        if (i2<0):
            return self._read_chunk_1d(i1,N1)
        elif (i3<0):
            if N1 != self._header.dims[0]:
                print ("Unable to support N1 {} != {}".format(N1,self._header.dims[0]))
                return None
            X=self._read_chunk_1d(i1+N1*i2,N1*N2)
            return np.reshape(X,(N1,N2),order='F')
        else:
            if N1 != self._header.dims[0]:
                print ("Unable to support N1 {} != {}".format(N1,self._header.dims[0]))
                return None
            if N2 != self._header.dims[1]:
                print ("Unable to support N2 {} != {}".format(N2,self._header.dims[1]))
                return None
            X=self._read_chunk_1d(i1+N1*i2+N1*N2*i3,N1*N2*N3)
            return np.reshape(X,(N1,N2,N3),order='F')
    def _read_chunk_1d(self,i,N):
        f=open(self._path,"rb")
        try:
            f.seek(self._header.header_size+self._header.num_bytes_per_entry*i)
            ret=np.fromfile(f,dtype=self._header.dt,count=N)
            f.close()
            return ret
        except Exception as e: # catch *all* exceptions
            print (e)
            f.close()
            return None

class DiskWriteMda:
    def __init__(self,path,dims,dt='float64'):
        self._path=path
        self._header=MdaHeader(dt,dims)
        _write_header(path, self._header)
        self._write_chunk_1d(np.zeros(np.product(dims)),np.product(dims))
    def N1(self):
        return self._header.dims[0]
    def N2(self):
        return self._header.dims[1]
    def N3(self):
        return self._header.dims[2]
    def writeChunk(self,X,i1=-1,i2=-1,i3=-1):
        print("Writing chunk {} {} {}".format(i1,i2,i3),X.shape)
        if (len(X.shape)>=2):
            N1=X.shape[0]
        else:
            N1=1
        if (len(X.shape)>=2):
            N2=X.shape[1]
        else:
            N2=1
        if (len(X.shape)>=3):
            N3=X.shape[2]
        else:
            N3=1
        if (i2<0):
            return self._write_chunk_1d(X,i1)
        elif (i3<0):
            if N1 != self._header.dims[0]:
                print ("Unable to support DiskWriteMda N1 {} != {}".format(N1,self._header.dims[0]))
                return None
            return self._write_chunk_1d(X.ravel(order='F'),i1+N1*i2)
        else:
            if N1 != self._header.dims[0]:
                print ("Unable to support DiskWriteMda N1 {} != {}".format(N1,self._header.dims[0]))
                return None
            if N2 != self._header.dims[1]:
                print ("Unable to support DiskWriteMda N2 {} != {}".format(N2,self._header.dims[1]))
                return None
            return self._write_chunk_1d(X.ravel(order='F'),i1+N1*i2+N1*N2*i3)
    def _write_chunk_1d(self,X,i):
        print('_write_chunk_1d',X.shape,i)
        N=X.size
        f=open(self._path,"ab")
        try:
            print(self._header.num_bytes_per_entry)
            print('Writing data to file at position {} values: {}'.format(self._header.header_size+self._header.num_bytes_per_entry*i,X))
            f.seek(self._header.header_size+self._header.num_bytes_per_entry*i)
            X.astype(self._header.dt).tofile(f)
            f.close()
            return True
        except Exception as e: # catch *all* exceptions
            print (e)
            f.close()
            return False

def _dt_from_dt_code(dt_code):
    if dt_code == -2:
        dt='uint8'
    elif dt_code == -3:
        dt='float32'
    elif dt_code == -4:
        dt='int16'
    elif dt_code == -5:
        dt='int32'
    elif dt_code == -6:
        dt='uint16'
    elif dt_code == -7:
        dt='float64'
    elif dt_code == -8:
        dt='uint32'
    else:
        dt=None
    return dt

def _dt_code_from_dt(dt):
    if dt == 'uint8':
        return -2
    if dt == 'float32':
        return -3
    if dt == 'int16':
        return -4
    if dt == 'int32':
        return -5
    if dt == 'uint16':
        return -6
    if dt == 'float64':
        return -7
    if dt == 'uint32':
        return -8
    return None

def get_num_bytes_per_entry_from_dt(dt):
    if dt == 'uint8':
        return 1
    if dt == 'float32':
        return 4
    if dt == 'int16':
        return 2
    if dt == 'int32':
        return 4
    if dt == 'uint16':
        return 2
    if dt == 'float64':
        return 8
    if dt == 'uint32':
        return 4
    return None

def _read_header(path):
    f=open(path,"rb")
    try:
        dt_code=_read_int32(f)
        num_bytes_per_entry=_read_int32(f)
        num_dims=_read_int32(f)
        uses64bitdims=False
        if (num_dims<0):
            uses64bitdims=True
            num_dims=-num_dims
        if (num_dims<1) or (num_dims>6): # allow single dimension as of 12/6/17
            print ("Invalid number of dimensions: {}".format(num_dims))
            f.close()
            return None
        dims=[]
        dimprod=1
        if uses64bitdims:
            for j in range(0,num_dims):
                tmp0=_read_int64(f)
                dimprod=dimprod*tmp0
                dims.append(tmp0)
        else:
            for j in range(0,num_dims):
                tmp0=_read_int32(f)
                dimprod=dimprod*tmp0
                dims.append(tmp0)
        dt=_dt_from_dt_code(dt_code)
        if dt is None:
            print ("Invalid data type code: {}".format(dt_code))
            f.close()
            return None
        H=MdaHeader(dt,dims)
        if (uses64bitdims):
            H.uses64bitdims=True
            H.header_size=3*4+H.num_dims*8
        f.close()
        return H
    except Exception as e: # catch *all* exceptions
        print (e)
        f.close()
        return None

def _write_header(path,H,rewrite=False):
    if rewrite:
        f=open(path,"r+b")
    else:
        f=open(path,"wb")
    try:
        _write_int32(f,H.dt_code)
        _write_int32(f,H.num_bytes_per_entry)
        if H.uses64bitdims:
            _write_int32(f,-H.num_dims)
            for j in range(0,H.num_dims):
                _write_int64(f,H.dims[j])
        else:
            _write_int32(f,H.num_dims)
            for j in range(0,H.num_dims):
                _write_int32(f,H.dims[j])
        f.close()
        return True
    except Exception as e: # catch *all* exceptions
        print (e)
        f.close()
        return False

def readmda(path):
    H=_read_header(path)
    if (H is None):
        print ("Problem reading header of: {}".format(path))
        return None
    ret=np.array([])
    f=open(path,"rb")
    try:
        f.seek(H.header_size)
        #This is how I do the column-major order
        ret=np.fromfile(f,dtype=H.dt,count=H.dimprod)
        ret=np.reshape(ret,H.dims,order='F')
        f.close()
        return ret
    except Exception as e: # catch *all* exceptions
        print (e)
        f.close()
        return None

def writemda32(X,fname):
    return _writemda(X,fname,'float32')

def writemda64(X,fname):
    return _writemda(X,fname,'float64')

def writemda8(X,fname):
    return _writemda(X,fname,'uint8')

def writemda32i(X,fname):
    return _writemda(X,fname,'int32')

def writemda32ui(X,fname):
    return _writemda(X,fname,'uint32')    

def writemda16i(X,fname):
    return _writemda(X,fname,'int16')    

def writemda16ui(X,fname):
    return _writemda(X,fname,'uint16')    

def _writemda(X,fname,dt):
    dt_code=0
    num_bytes_per_entry=get_num_bytes_per_entry_from_dt(dt)
    dt_code=_dt_code_from_dt(dt)
    if dt_code is None:
        print ("Unexpected data type: {}".format(dt))
        return False

    f=open(fname,'wb')
    try:
        _write_int32(f,dt_code)
        _write_int32(f,num_bytes_per_entry)
        _write_int32(f,X.ndim)
        for j in range(0,X.ndim):
            _write_int32(f,X.shape[j])
        #This is how I do column-major order
        A=np.reshape(X,X.size,order='F').astype(dt)
        A.tofile(f)
        f.close()
        return True
    except Exception as e: # catch *all* exceptions
        print (e)
        f.close()
        return False

def appendmda(X,path):
    H=_read_header(path)
    if (H is None):
        print ("Problem reading header of: {}".format(path))
        return None
    if (len(H.dims) != len(X.shape)):
        print ("Incompatible number of dimensions in appendmda",H.dims,X.shape)
        return None
    num_entries_old=np.product(H.dims)
    num_dims=len(H.dims)
    for j in range(num_dims-1):
        if (X.shape[j]!=X.shape[j]):
            print ("Incompatible dimensions in appendmda",H.dims,X.shape)
            return None
    H.dims[num_dims-1]=H.dims[num_dims-1]+X.shape[num_dims-1]
    try:
        _write_header(path,H,rewrite=True)
        f=open(path,"r+b")
        f.seek(H.header_size+H.num_bytes_per_entry*num_entries_old)
        A=np.reshape(X,X.size,order='F').astype(H.dt)
        A.tofile(f)
        f.close()
    except Exception as e: # catch *all* exceptions
        print (e)
        f.close()
        return False
    
def _read_int32(f):
    return struct.unpack('<i',f.read(4))[0]
    
def _read_int64(f):
    return struct.unpack('<q',f.read(8))[0]

def _write_int32(f,val):
    f.write(struct.pack('<i',val))
    
def _write_int64(f,val):
    f.write(struct.pack('<q',val))

def mdaio_test():
    M=4
    N=12
    X=np.ndarray((M,N))
    for n in range(0,N):
        for m in range(0,M):
            X[m,n]=n*10+m
    writemda32(X,'tmp1.mda')
    Y=readmda('tmp1.mda')
    print (Y)
    print (np.absolute(X-Y).max())
    Z=DiskReadMda('tmp1.mda')
    print (Z.readChunk(i1=0,i2=4,N1=M,N2=N-4))

    A=DiskWriteMda('tmpA.mda',(M,N))
    A.writeChunk(Y,i1=0,i2=0)
    B=readmda('tmpA.mda')
    print (B.shape)
    print (B)

#mdaio_test()