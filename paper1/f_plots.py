import numpy as np
import matplotlib.pyplot as plt
from resonantstate.ell2SFM import *

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


def plt_rob(ax,x,y,xerr1,xerr2,yerr1,yerr2,starname,keep_reason,edgecolor=np.array([]),alpha_text=0,yabs=False,coef_ms=1.7):
    """_summary_

    Parameters
    ----------
    ax : _type_
        _description_
    x : array of floats
        _description_
    y : array of floats
        _description_
    xerr1 : array of floats
        _description_
    xerr2 : array of floats
        _description_
    yerr1 : array of floats
        _description_
    yerr2 : array of floats
        _description_
    starname : array of strs
        _description_
    keep_reason : array of floats
        _description_
    edgecolor : str, optional
        _description_, by default np.array([])
    alpha_text : int, optional
        _description_, by default 0
    yabs : bool, optional
        _description_, by default False
    coef_ms : float, optional
        _description_, by default 1.7

    Returns
    -------
    ax : ax 
        
    """
    if edgecolor.size==0:
        edgecolor=np.full([x.size], "black", dtype=(str, 5))

    for xi, yi, xerre1,xerre2,yerre1,yerre2, name,krea,edge in zip(x,y,xerr1,xerr2,yerr1,yerr2,starname,keep_reason,edgecolor):

        if krea<= -1:
            color='grey'
            alpha=0.3
            ms=3*coef_ms
        elif krea==0:
            color='C0'
            alpha=.7
            ms=3*coef_ms
        elif krea>0:
            color='C1'
            alpha=1
            ms=3*coef_ms


        if yabs and yi<0:
            yi,yerre1,yerre2=abs_err(yi,yerre1,yerre2)

        if alpha>0:
            ax.errorbar(xi, yi, xerr=np.array([xi-xerre1,xerre2-xi]).reshape(2,1),yerr=np.abs(np.array([yi-yerre1,yerre2-yi])).reshape(2,1), fmt='o', color=color, ecolor=color, capsize=1,alpha=alpha,ms=ms,zorder=200+krea,markeredgecolor=edge)
        
        if alpha_text>0:
            ax.text(xi, (yi),name,alpha=alpha_text)

    return ax


def plt_color(ax,x,y,xerr1,xerr2,yerr1,yerr2,col,starname,Tobsname,keep_reason,cname,sm,cmap,norm,alpha_text=0,yabs=False,plt_cbar=True):

    for xi, yi, xerre1,xerre2,yerre1,yerre2, ci, starname, name, krea in zip(x,y,xerr1,xerr2,yerr1,yerr2,col,starname,Tobsname,keep_reason):

        color = cmap(norm(ci))
    
        if krea<0:
            alpha=0.1
        elif krea==0:
            alpha=.5
        else:
            alpha=1

        if yabs and yi<0:
            yi,yerre1,yerre2=abs_err(yi,yerre1,yerre2)

        if alpha>0:
            ax.errorbar(xi, yi, xerr=np.array([xi-xerre1,xerre2-xi]).reshape(2,1),yerr=np.array([yi-yerre1,yerre2-yi]).reshape(2,1), fmt='o', color=color, ecolor='grey', capsize=2,alpha=alpha,zorder=200+krea,markeredgecolor='black')

        if alpha_text>0 and krea>=0:    
            ax.text(xi, (yi),starname,alpha=alpha_text)

    # Add colorbar to the correct axis
    if plt_cbar: 
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label(cname)

        return ax,cbar
    
    else:
        return ax
    


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

  
    width_lines=.2
    alpha_area=.2
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
    label ="$\left\{  \\begin{array}{l} \\textnormal{Overstable libration} \\\ \\textnormal{Turbulence} \\\ \\textnormal{Mild instabilities} \end{array}  \\right.$"#  \\ \\


    #cyan
    Id=np.where(delta_grid>-2)[0]
    DMMR_grid=fixp_grid**2-3*delta_grid

    Id2=np.where(delta_grid>0)[0]


    y_rentrance_1=np.arange(fixp_grid[Id][-1],fixp_grid[Id][-1]+1,.01)
    x_rentrance_1=y_rentrance_1**2+3*2


    Id3=np.where((delta_grid<=1)&(delta_grid>-2))[0][::-1]
    x_rentrance_2=(fixp_grid[Id3]+1)**2-3*delta_grid[Id3]#DMMR_grid[Id3]
    y_rentrance_2=fixp_grid[Id3]+1


    Iddeltap=np.where(delta_grid>=1)[0][::-1]


    xfill=np.concatenate((DMMR_grid[Id],x_rentrance_1,x_rentrance_2,sepp_grid[Iddeltap]**2-3*delta_grid[Iddeltap]))
    yfill=np.concatenate((fixp_grid[Id],y_rentrance_1,y_rentrance_2,sepp_grid[Iddeltap]))

    ax.fill(xfill,yfill,alpha=alpha_area,color='cyan',label=label)


    #green
    Id = np.where(delta_grid>-2.83/3)[0]
    ax.fill_betweenx(fixp_grid[Id],fixp_grid[Id]**2-3*delta_grid[Id],fixp_grid[Id]**2-3*delta_grid[Id]+width_lines*5,color='green',label='Laminar/Stable libration',alpha=alpha_line)
    Id = np.where((delta_grid<-2.83/3)&(delta_grid>-2))[0]
    print(Id.size,'ID-size')
    ax.fill_between(fixp_grid[Id]**2-3*delta_grid[Id],fixp_grid[Id],fixp_grid[Id]+width_lines,color='green',alpha=alpha_line)




    #blue
    Id_blue_upper = np.where(delta_grid<-2)[0]
    label ="$\left\{  \\begin{array}{l} \\textnormal{Tides} \\\ \\textnormal{Planetesimal disk} \end{array}  \\right.$"#  \\ \\


    x_blue_upper=fixp_grid[Id_blue_upper]**2-3*delta_grid[Id_blue_upper]
    y_blue_upper=fixp_grid[Id_blue_upper]
    ax.fill_between(x_blue_upper,y_blue_upper,y_blue_upper+width_lines,color='blue',label=label,alpha=alpha_line)

    Id = np.where(np.abs(delta_grid+1.35)<.15)[0]
    ax.fill_between(fixp_grid[Id]**2-3*delta_grid[Id],fixp_grid[Id],fixp_grid[Id]+width_lines,color='blue',alpha=alpha_line)
    Id = np.where(np.abs(delta_grid+.75)<.15)[0]
    ax.fill_between(fixp_grid[Id]**2-3*delta_grid[Id],fixp_grid[Id],fixp_grid[Id]+width_lines,color='blue',alpha=alpha_line)

    Id = np.where(delta_grid>2.7)[0]
    ax.fill_between(intc_grid[Id]**2-3*delta_grid[Id],intc_grid[Id],intc_grid[Id]-width_lines,color='blue',alpha=alpha_line)

    Id = np.where(np.abs(delta_grid-2.2)<.15)[0]
    ax.fill_between(intc_grid[Id]**2-3*delta_grid[Id],intc_grid[Id],intc_grid[Id]-width_lines,color='blue',alpha=alpha_line)
    Id = np.where(np.abs(delta_grid-1.6)<.15)[0]
    ax.fill_between(intc_grid[Id]**2-3*delta_grid[Id],intc_grid[Id],intc_grid[Id]-width_lines,color='blue',alpha=alpha_line)


    #yellow lower
    Id_delp1=np.where(delta_grid>1)

    yfixpint=intc_grid[Id_delp1]
    xfixpint=intc_grid[Id_delp1]**2-3*delta_grid[Id_delp1]

    yhyp=hyp_grid[Id_delp1]
    xhyp=hyp_grid[Id_delp1]**2-3*delta_grid[Id_delp1]

    xfill=np.concatenate((xfixpint,xhyp[::-1],np.array([DMMR_lim[0]])))
    yfill=np.concatenate((yfixpint,yhyp[::-1],np.array([er_lim[0]])))

    ax.fill(xfill,yfill,alpha=alpha_area,color='yellow',label='Post-disc instabilities')

    #yellow upper
    Idsepp=np.where(sepp_grid)[0]
    ysepp=sepp_grid[Idsepp]
    xsepp=sepp_grid[Idsepp]**2-3*delta_grid[Idsepp]


    xfill=np.concatenate((xsepp,x_rentrance_2[::-1],x_rentrance_1[::-1],x_blue_upper,np.array([DMMR_lim[1]])))
    yfill=np.concatenate((ysepp,y_rentrance_2[::-1],y_rentrance_1[::-1],y_blue_upper,np.array([er_lim[1]])))

    ax.fill(xfill,yfill,alpha=alpha_area,color='yellow')

    return ax