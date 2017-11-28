from mlpy import DiskReadMda
import time

class TimeseriesChunkInfo:
    def __init__(self):
        self.t1a=0
        self.t2a=0
        self.t1=0
        self.t2=0
        self.size=0

class TimeseriesChunkReader:
    def __init__(self, chunk_size=0, chunk_size_mb=0, overlap_size=0, t1=-1, t2=-1, verbose=True):
        # Note that the actual chunk size will be the maximum of chunk_size,overlap_size and chunk_size_mb*1e6/(M*4)
        self._chunk_size=chunk_size
        self._chunk_size_mb=chunk_size_mb
        self._overlap_size=overlap_size
        self._t1=t1
        self._t2=t2
        self._elapsed_reading=0
        self._elapsed_running=0
        self._verbose=verbose
    def run(self, mdafile_path_or_diskreadmda, func):
        if (type(mdafile_path_or_diskreadmda)==str):
            X=DiskReadMda(mdafile_path_or_diskreadmda)
        else:
            X=mdafile_path_or_diskreadmda
        M,N = X.N1(),X.N2()
        cs=max([self._chunk_size,int(self._chunk_size_mb*1e6/(M*4)),M])        
        if self._t1<0:
            self._t1=0
        if self._t2<0:
            self._t2=N-1
        t=self._t1
        while t <= self._t2:
            t1=t
            t2=min(self._t2,t+cs-1)
            s1=max(0,t1-self._overlap_size)
            s2=min(N-1,t2+self._overlap_size)
            
            timer=time.time()
            chunk=X.readChunk(i1=0, N1=M, i2=s1, N2=s2-s1+1)
            self._elapsed_reading+=time.time()-timer
            
            info=TimeseriesChunkInfo()
            info.t1=t1
            info.t2=t2
            info.t1a=t1-s1
            info.t2a=t2-s1
            info.size=t2-t1+1
            
            timer=time.time()
            if not func(chunk, info):
                return False
            self._elapsed_running+=time.time()-timer                
                
            t=t+cs
        if self._verbose:
            print('Elapsed for TimeseriesChunkReader: %g sec reading, %g sec running' % (self._elapsed_reading,self._elapsed_running))
        return True
    def elapsedReading(self):
        return self._elapsed_reading
    def elapsedRunning(self):
        return self._elapsed_running