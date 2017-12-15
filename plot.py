import matplotlib.pyplot as plt
import numpy as np

def loadit(i,ny=200,nx=72):
    fname = 'conv.{:04d}.tab'.format(i)
    dat = np.loadtxt(fname)
    p = dat[:,-1].reshape(ny,nx)
    d = dat[:,4].reshape(ny,nx)
    vz = dat[:,6].reshape(ny,nx)
    y = dat[:,3].reshape(ny,nx)
    T = p/(d*.4)
    s = np.log(T*p**(-.4))
    s -= s[0,:].mean()
    return y[:,0],T,d,p,vz,s

def plot2d(q,fig=None,ax=None,symmetric=False):
    if ax is None:
        fig,ax = plt.subplots()
    kargs={}
    if symmetric:
        vmax =  abs(q).max()
        vmin = - vmax
        kargs = {'vmin':vmin,'vmax':vmax}

    ax.imshow(q,origin='lower',extent=(-.5,.5,-1,2),cmap='coolwarm',**kargs)

def plotit(y,T,d,p,vz,s,j=10,avg=True,fig=None,axes=None,**kargs):
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


def summary(ivals,figsize=(15,10),savefig=None,**kargs):
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(2,4)
    fig = plt.figure(figsize=figsize)

    axes = []
    for i in range(2):
        for j in range(2):
            axes.append(fig.add_subplot(gs[i,j]))
    axes = np.array(axes).reshape(2,2)
    ax2 = fig.add_subplot(gs[:,-2])
    ax3 = fig.add_subplot(gs[:,-1])

    plotit(*loadit(0),fig=fig,axes=axes,**kargs)
    try:
        nt = len(ivals)
    except TypeError:
        nt = 1
        ivals = [ivals]
    for i in ivals:
        y,T,d,p,vz,s = loadit(i)
        plotit(y,T,d,p,vz,s,fig=fig,axes=axes)

    plot2d(s-s.mean(axis=1)[:,np.newaxis],ax=ax2,fig=fig)
    plot2d(vz,ax=ax3,fig=fig,symmetric=True)
    ax2.set_title('$s-\\langle s \\rangle$')
    ax3.set_title('$v_z$')
    fig.tight_layout()
    if savefig is not None:
        fig.savefig(savefig,bbox_inches='tight')
    return fig,axes
