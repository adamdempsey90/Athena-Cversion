import matplotlib.pyplot as plt
import numpy as np

def loadit(i,ny=104,nx=40):
    fname = 'conv.{:04d}.tab'.format(i)
    dat = np.loadtxt(fname)
    p = dat[:,-1].reshape(ny,nx)
    d = dat[:,4].reshape(ny,nx)
    y = dat[:,3].reshape(ny,nx)
    T = p/(d*.4)
    s = np.log(T*p**(-.4))
    s -= s[0,:].mean()
    return y[:,0],T,d,p,s

def plot2d(q,fig=None,ax=None):
    if ax is None:
        fig,ax = plt.subplots()
    ax.imshow(q,origin='lower')

def plotit(y,T,d,p,s,j=10,avg=True,fig=None,axes=None):
    if axes is None:
        fig,axes = plt.subplots(2,2,sharex=True)

    if avg:
        axes[0,0].plot(y,T.mean(axis=1),'.-')
        axes[0,1].plot(y,p.mean(axis=1),'.-')
        axes[1,0].plot(y,d.mean(axis=1),'.-')
        axes[1,1].plot(y,s.mean(axis=1),'.-')
    else:
        axes[0,0].plot(y,T[:,j],'.-')
        axes[0,1].plot(y,p[:,j],'.-')
        axes[1,0].plot(y,d[:,j],'.-')
        axes[1,1].plot(y,s[:,j],'.-')

    for ax in axes.flatten():
        ax.minorticks_on()
    axes[0,0].set_ylabel('$T$')
    axes[0,1].set_ylabel('$P$')
    axes[1,0].set_ylabel('$\\rho$')
    axes[1,1].set_ylabel('$s$')
    for ax in axes[-1,:]:
        ax.set_xlabel('$z$')
    fig.tight_layout()
    return fig,axes


def summary(i):
    fig,axes = plotit(*loadit(0))
    plotit(*loadit(i),fig=fig,axes=axes)
    plot2d(loadit(i)[-1])
