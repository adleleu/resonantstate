import numpy as np
from resonantstate.SFMR  import *

def plt_area_paper2(ax,Delta_lim,X_lim,alpha_area=0.1,legend=False):


    bound=(max(np.abs(X_lim))**2+max(np.abs(Delta_lim)))/3

    delta_grid=-np.arange(-bound,bound,.01)
    fixp_grid=np.zeros(delta_grid.size)
    sepp_grid=np.zeros(delta_grid.size)*np.nan
    sepm_grid=np.zeros(delta_grid.size)*np.nan
    intc_grid=np.zeros(delta_grid.size)
    lib_boundary=np.zeros(delta_grid.size)
    intlib_boundary=np.zeros(delta_grid.size)
    hyp_grid=np.zeros(delta_grid.size)

    for k in range(delta_grid.size):
        Xmin, Xmax, Xres, Xint, Xhyp=topology(np.array([delta_grid[k]]))
        fixp_grid[k]=Xres
        sepp_grid[k]=Xmax
        sepm_grid[k]=Xmin
        intc_grid[k]=Xint
        hyp_grid[k]=Xhyp
        if delta_grid[k]<0:
            lib_boundary[k]=max(X1X2(0, 0, delta_grid[k]))
        else:
            lib_boundary[k]=Xmax
            intlib_b=min(X1X2(0, 0, delta_grid[k]))
            if intlib_b<-0.001:
                intlib_boundary[k]=intlib_b
            else:
                intlib_boundary[k]=Xhyp



    #formally resonant
    Iddelta1=np.abs(delta_grid+2).argmin()
    Idsepp=np.where(sepp_grid)[0]
    ysepp=sepp_grid[Idsepp]
    xsepp=sepp_grid[Idsepp]**2-3*(delta_grid[Idsepp]+1)

    Id_trap=np.where(np.abs(delta_grid+1)<1)[0]

    ysepp=np.append(ysepp,fixp_grid[Id_trap]+1)
    xsepp=np.append(xsepp,(fixp_grid[Id_trap]+1)**2-3*(delta_grid[Id_trap]+1))

    yfixp=fixp_grid[:Iddelta1]
    xfixp=yfixp**2-3*(delta_grid[:Iddelta1]+1)

    fyd1=np.arange(yfixp[-1],ysepp[-1],.01)
    fxd1=fyd1**2+3*(1)

    xfill=np.concatenate((xfixp,fxd1,xsepp[::-1]))
    yfill=np.concatenate((yfixp,fyd1,ysepp[::-1]))


    if legend:
        ax.fill(xfill, yfill, 'red',alpha=alpha_area,label='resonance-trapped',zorder=-200)
    else:
        ax.fill(xfill, yfill, 'red',alpha=alpha_area,zorder=-200)

    #external circulation 
    Idext=np.where(delta_grid<-2)[0]
    xfill=np.concatenate((xsepp,fixp_grid[Idext]**2-3*(delta_grid[Idext]+1),np.array([Delta_lim[1]])))
    yfill=np.concatenate((ysepp,fixp_grid[Idext],np.array([X_lim[1]])))

    if legend:
        ax.fill(xfill,yfill,alpha=alpha_area,color='orange',label='upper area',zorder=-200)
    else:
        ax.fill(xfill,yfill,alpha=alpha_area,color='orange',zorder=-200)

    #internal circulation

    Id_delp1=np.where(delta_grid>0)
    xfixp=intc_grid[Id_delp1]**2-3*(delta_grid[Id_delp1]+1)
    intc_grid[Id_delp1][::-1]

    yhyp=hyp_grid[Id_delp1]
    xhyp=hyp_grid[Id_delp1]**2-3*(delta_grid[Id_delp1]+1)
    xfill=np.concatenate((xfixp,xhyp[::-1],np.array([Delta_lim[0]])))
    yfill=np.concatenate((intc_grid[Id_delp1],yhyp[::-1],np.array([X_lim[0]])))

    if legend:
        ax.fill(xfill,yfill,alpha=alpha_area,color='yellow',label='lower area',zorder=-200)
    else:
        ax.fill(xfill,yfill,alpha=alpha_area,color='yellow',zorder=-200)
