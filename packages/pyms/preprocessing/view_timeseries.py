import matplotlib.pyplot as plt
import numpy as np

def view_timeseries(timeseries,trange=None,channels=None,samplerate=30000,title='',fig_size=[18,6]):
    #timeseries=mls.loadMdaFile(timeseries)
    
    set_fig_size(fig_size[0],fig_size[1])
    
    M=timeseries.shape[0]
    N=timeseries.shape[1]
    
    channel_colors=_get_channel_colors(M)
    
    if not trange:
        trange=[0,np.minimum(1000,N)]
    if not channels:
        channels=np.arange(M).tolist()
        
    X=timeseries[channels][:,int(trange[0]):int(trange[1])]
    
    spacing_between_channels=np.max(np.abs(X.ravel()))
    
    y_offset=0
    for m in range(len(channels)):
        A=X[m,:]
        plt.plot(np.arange(X.shape[1]),A+y_offset,color=channel_colors[channels[m]])
        y_offset-=spacing_between_channels
    
    ax=plt.gca()
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    
    if title:
        plt.title(title,fontsize=title_fontsize)
    
    plt.show()
    return ax;

def set_fig_size(W,H):
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = W
    fig_size[1] = H
    plt.rcParams["figure.figsize"] = fig_size
    
def _get_channel_colors(M):
    cm = plt.get_cmap('gist_ncar')
    channel_colors=[]
    for m in range(M):
        channel_colors.append(cm(1.0*(m+0.5)/M))
    return channel_colors