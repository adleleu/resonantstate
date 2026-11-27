import numpy as np
from resonantstate.SFMR  import *
import matplotlib.pyplot as plt


def abs_err(yi,yerre1,yerre2):
    """_summary_

    Parameters
    ----------
    yi :array of floats
        median values
    yerre1 :array of floats
        .16 quantile
    yerre2 :array of floats
        .84 quantile

    Returns
    -------
    yi :array of floats
        median values
    yerre1 :array of floats
        .16 quantile
    yerre2 :array of floats
        .84 quantile

        redefine errors if the median is negative 
    """
    if yi<0:
        yi=-yi
        syerre1=yerre1
        yerre1=abs(yerre2)
        yerre2=abs(syerre1)
    return yi,yerre1,yerre2


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



def theory_background(ax,DMMR_lim,er_lim):
    """colour the background of the SFM to describe which physcial process produce pairs of planet in a given area.

    Parameters
    ----------
    ax : matplotlib.axes
        _description_
    DMMR_lim : tuple of floats (DMMR_min,DMMR_max)
        _description_
    er_lim : tuple of floats (er_min,er_max)
        _description_
    """

  
    width_lines=.3
    alpha_area=.4
    alpha_line=.4
    bound=(max(np.abs(er_lim))**2+max(np.abs(DMMR_lim)))/3

    delta_grid=-np.arange(-bound,bound,.01)
    fixp_grid=np.zeros(delta_grid.size)
    sepp_grid=np.zeros(delta_grid.size)
    intc_grid=np.zeros(delta_grid.size)
    lib_boundary=np.zeros(delta_grid.size)
    hyp_grid=np.zeros(delta_grid.size)
    intlib_boundary=np.zeros(delta_grid.size)

    for k in range(delta_grid.size):
        Xmin, Xmax, Xres, Xint, Xhyp=topology(np.array([delta_grid[k]]))
        fixp_grid[k]=Xres
        sepp_grid[k]=Xmax
        intc_grid[k]=Xint
        hyp_grid[k]=Xhyp
        if delta_grid[k]<1:
            lib_boundary[k]=max(X1X2(0, 0, delta_grid[k]))
        else:
            lib_boundary[k]=Xmax
            intlib_b=min(X1X2(0, 0, delta_grid[k]))
            if intlib_b<-0.001:
                intlib_boundary[k]=intlib_b
            else:
                intlib_boundary[k]=Xhyp

    plt.rcParams["text.usetex"] = True
    label ="$\left\{  \\begin{array}{l} \\textnormal{Overstable libr.} \\\ \\textnormal{Mild instabilities} \end{array}  \\right.$"#   \\textnormal{Turbulence} \\\


    #Red - resonant
    Id=np.where(delta_grid>-2)[0]
    DMMR_grid=fixp_grid**2-3*(delta_grid+1)
    y_rentrance_1=np.arange(fixp_grid[Id][-1],fixp_grid[Id][-1]+1,.01)
    x_rentrance_1=y_rentrance_1**2+3
    Id3=np.where((delta_grid<=1)&(delta_grid>-2))[0][::-1]
    x_rentrance_2=(fixp_grid[Id3]+1)**2-3*(1+delta_grid[Id3])#DMMR_grid[Id3]
    y_rentrance_2=fixp_grid[Id3]+1
    Iddeltap=np.where(delta_grid>=1)[0][::-1]
    xfill=np.concatenate((DMMR_grid[Id],x_rentrance_1,x_rentrance_2,sepp_grid[Iddeltap]**2-3*(1+delta_grid[Iddeltap])))
    yfill=np.concatenate((fixp_grid[Id],y_rentrance_1,y_rentrance_2,sepp_grid[Iddeltap]))
    ax.fill(xfill,yfill,alpha=alpha_area,color='C3',label=label,zorder=-1e4,linewidth=0.0)


    #Blue fix point resonant
    label='$\left\{  \\begin{array}{l}\\textnormal{Laminar} \\\ \\textnormal{Stable libration}\end{array}  \\right.$'#'$\\textnormal{Stable libration}$'#
    Id = np.where(delta_grid>-2.83/3)[0]
    ax.fill_betweenx(fixp_grid[Id],fixp_grid[Id]**2-3*(1+delta_grid[Id]),fixp_grid[Id]**2-3*(1+delta_grid[Id])+width_lines*5,color='white',alpha=1,zorder=-1e2-1)
    ax.fill_betweenx(fixp_grid[Id],fixp_grid[Id]**2-3*(1+delta_grid[Id]),fixp_grid[Id]**2-3*(1+delta_grid[Id])+width_lines*5,color='C0',label=label,alpha=alpha_line,zorder=-1e2)
    Id = np.where((delta_grid<-2.83/3)&(delta_grid>-2))[0]
    print(Id.size,'ID-size')
    ax.fill_between(fixp_grid[Id]**2-3*(1+delta_grid[Id]),fixp_grid[Id],fixp_grid[Id]+width_lines*2,color='white',alpha=1,zorder=-1e2-1)
    ax.fill_between(fixp_grid[Id]**2-3*(1+delta_grid[Id]),fixp_grid[Id],fixp_grid[Id]+width_lines*2,color='C0',alpha=alpha_line,zorder=-1e2)
    
    #pink
    #pink="#CB87A9"
    Id_blue_upper = np.where(delta_grid<-2)[0]
    label ="$\left\{  \\begin{array}{l} \\textnormal{Tides} \\\ \\textnormal{Planetesimal disk} \end{array}  \\right.$"#  \\ \\


    x_blue_upper=fixp_grid[Id_blue_upper]**2-3*(1+delta_grid[Id_blue_upper])
    y_blue_upper=fixp_grid[Id_blue_upper]
    ax.fill_between(x_blue_upper,y_blue_upper,y_blue_upper+width_lines,color='white',alpha=1,zorder=-1e2-1)
    
    ax.fill_between(x_blue_upper,y_blue_upper,y_blue_upper+width_lines,color='C6',label=label,alpha=alpha_line,zorder=-1e2)

    Id = np.where(np.abs(delta_grid+1.35)<.15)[0]
    ax.fill_between(fixp_grid[Id]**2-3*(1+delta_grid[Id]),fixp_grid[Id],fixp_grid[Id]+width_lines,color='C6',alpha=alpha_line,zorder=-1e2)
    Id = np.where(np.abs(delta_grid+.75)<.15)[0]
    ax.fill_between(fixp_grid[Id]**2-3*(1+delta_grid[Id]),fixp_grid[Id],fixp_grid[Id]+width_lines,color='C6',alpha=alpha_line,zorder=-1e2)

    Id = np.where(delta_grid>2.7)[0]
    ax.fill_between(intc_grid[Id]**2-3*(1+delta_grid[Id]),intc_grid[Id],intc_grid[Id]-width_lines,color='white',alpha=1,zorder=-1e2-1)
    ax.fill_between(intc_grid[Id]**2-3*(1+delta_grid[Id]),intc_grid[Id],intc_grid[Id]-width_lines,color='C6',alpha=alpha_line,zorder=-1e2)


    Id = np.where(np.abs(delta_grid-2.2)<.15)[0]
    ax.fill_between(intc_grid[Id]**2-3*(1+delta_grid[Id]),intc_grid[Id],intc_grid[Id]-width_lines,color='C6',alpha=alpha_line,zorder=-1e2)
    Id = np.where(np.abs(delta_grid-1.6)<.15)[0]
    ax.fill_between(intc_grid[Id]**2-3*(1+delta_grid[Id]),intc_grid[Id],intc_grid[Id]-width_lines,color='C6',alpha=alpha_line,zorder=-1e2)


    #Gold lower
    Id_delp1=np.where(delta_grid>0)

    yfixpint=intc_grid[Id_delp1]
    xfixpint=intc_grid[Id_delp1]**2-3*(1+delta_grid[Id_delp1])

    yhyp=hyp_grid[Id_delp1]
    xhyp=hyp_grid[Id_delp1]**2-3*(1+delta_grid[Id_delp1])

    xfill=np.concatenate((xfixpint,xhyp[::-1],np.array([DMMR_lim[0]])))
    yfill=np.concatenate((yfixpint,yhyp[::-1],np.array([er_lim[0]])))

    ax.fill(xfill,yfill,alpha=alpha_area,color='gold',label='$\\textnormal{Post-disc instabilities}$',zorder=-1e3)

    #Gold upper
    Idsepp=np.where(sepp_grid)[0]
    ysepp=sepp_grid[Idsepp]
    xsepp=sepp_grid[Idsepp]**2-3*(1+delta_grid[Idsepp])


    xfill=np.concatenate((xsepp,x_rentrance_2[::-1],x_rentrance_1[::-1],x_blue_upper,np.array([DMMR_lim[1]])))
    yfill=np.concatenate((ysepp,y_rentrance_2[::-1],y_rentrance_1[::-1],y_blue_upper,np.array([er_lim[1]])))

    ax.fill(xfill,yfill,alpha=alpha_area,color='gold',zorder=-1e3)


    #Green dashes

    Id = np.where(delta_grid<.5)[0]
    px1=fixp_grid[Id]**2-3*(1+delta_grid[Id])
    py1=fixp_grid[Id]
    px2=(fixp_grid[Id]+1.2)**2-3*(1+delta_grid[Id])
    px2=px2[::-1]
    py2=py1[::-1]+2

    px=np.append(px1,px2)
    py=np.append(py1,py2)

    plt.rcParams["hatch.linewidth"] = 1.0
    ax.fill(px, py, facecolor="none", edgecolor='white', hatch="//", linewidth=0.0,zorder=-11,alpha=1)   
    ax.fill(px, py, facecolor="none", edgecolor='C2', hatch="//", linewidth=0.0,zorder=-10,label='$\\textnormal{Turbulent disk}$',alpha=1)   


    #Green dashes lower
    Id = np.where(delta_grid>0)[0]
    px1=intc_grid[Id]**2-3*(1+delta_grid[Id])
    py1=intc_grid[Id]
    px2=px1[::-1]#
    py2=py1[::-1]-1.5

    px=np.append(px1,px2)
    py=np.append(py1,py2)
    plt.rcParams["hatch.linewidth"] = 1.0
    ax.fill(px, py, facecolor="none", edgecolor='white', hatch="//", linewidth=0.0,zorder=-11,alpha=1)  
    ax.fill(px, py, facecolor="none", edgecolor='C2', hatch="//", linewidth=0.0,zorder=-10,alpha=1)   

    
    return ax



def plt_rob(ax,x,y,xerr1,xerr2,yerr1,yerr2,starname,keep_reason,edgecolor=np.array([]),alpha_text=0,yabs=False,coef_ms=1.8,alpha_grey=.5,color_rob='C1'):

    if edgecolor.size==0:
        edgecolor=np.full([x.size], "black", dtype=(str, 5))

    for xi, yi, xerre1,xerre2,yerre1,yerre2, name,krea,edge in zip(x,y,xerr1,xerr2,yerr1,yerr2,starname,keep_reason,edgecolor):
    
        if krea>0:
            color=color_rob
            alpha=1
            ms=3*coef_ms
            zorder=1e100
        else:

            color='grey'
            alpha=alpha_grey
            ms=1*coef_ms      
            zorder=-1  

        if yabs and yi<0:
            yi,yerre1,yerre2=abs_err(yi,yerre1,yerre2)

        if alpha>0:
            markers, caps, bars=ax.errorbar(xi, yi, xerr=np.array([xi-xerre1,xerre2-xi]).reshape(2,1),yerr=np.abs(np.array([yi-yerre1,yerre2-yi])).reshape(2,1), fmt='o', color=color, ecolor=color, capsize=1,alpha=alpha,ms=ms,zorder=zorder,markeredgecolor=edge)
            if color=='grey':
                for b in bars:  
                    b.set_alpha(0.3)
        if (alpha_text>0) :
            ax.text(xi, (yi),name,alpha=alpha_text)

    return ax
