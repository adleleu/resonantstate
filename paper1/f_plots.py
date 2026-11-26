import numpy as np


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


    

