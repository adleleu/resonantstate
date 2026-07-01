import numpy as np
from resonantstate.constants import *
from resonantstate.ell2SFM import *

Tp=np.arange(1,20)


def _dist_Dpp1p(Pr):


    if Pr>2.5 or Pr<1.1:
        return np.inf,-1,False
    
    else :
        dist_H= -2*(f1s[Tp-1]/3/Tp**(1/2))**(-2/3)*(1-(Pr*Tp/(Tp+1))**(1/3))
        Id_closest=np.abs(dist_H).argmin()

        dist_H_best,p_best=dist_H[Id_closest],Tp[Id_closest]

        if dist_H_best>0:
            Pr_test=(2-(p_best/(1+p_best)*Pr)**(1/3))**3*(p_best+1)/p_best
            dist_H_test= -2*(f1s[Tp-1]/3/Tp**(1/2))**(-2/3)*(1-(Pr_test*Tp/(Tp+1))**(1/3))
            Id_closest_test=np.abs(dist_H_test).argmin()
            
            if Tp[Id_closest_test]==p_best:
                sym=True
            else:
                sym=False
        else:
            sym=True

        
        return dist_H[Id_closest],Tp[Id_closest],sym


def dist_Dpp1p(Pr):
    """deal with different entry shape for Pr

    Parameters
    ----------
    Pr : float or array
        period ratios

    Returns
    -------
    dist : float

    p : integer

    sym: boolean
        
    """

    if isinstance(Pr,np.ndarray):
        if len(Pr.shape)==1:
            dist_H,Tp,sym=np.zeros(Pr.size),np.zeros(Pr.size),np.zeros(Pr.size,dtype=bool)
            for k,pr in enumerate(Pr):
                dist_H[k],Tp[k],sym[k]=_dist_Dpp1p(pr)
            return dist_H,Tp,sym
    else:
        return  _dist_Dpp1p(Pr)


def get_PRM(samples,df_ana):
    nb_planets=len(df_ana['planets_list'])

    TP=np.zeros(nb_planets)
    TR=np.zeros(nb_planets)
    TM=np.zeros(nb_planets)

    for k in range(nb_planets):
        TP[k]=samples['period_days_'+str(k)].mean()
        TR[k]=(samples['radius_planet_star_ratio_'+str(k)]*samples['radius_star_r_sun']/(Rearth/Rsun)).mean()
        TM[k]=(samples['mass_planet_star_ratio_'+str(k)].values*samples['mass_star_m_sun'].values/(Mearth/Msun)).mean()

    Idsort=TP.argsort()

    return TP,TR,TM,nb_planets,Idsort



def get_SFM_quantities(samples,pair,p=None,sampling=1):
    """call functions of the package to extract the quantities used in the first paper : delta, DMMR, er, Der, es, IR, p

    Parameters
    ----------
    samples : pandas dataframe
        posterior of a dynamical analysis
    pair : tuple of ints
        identifiant of the pair to consider
    p : int, optional
        p+1:p first order MMR to consider
    sampling : int, optional
        thinning of the samples

    Returns
    -------
    delta : np.array of floats of the lengh of the thinned samples
        delta parameter of the second fundamental model for resonance (Henrard and Lemaitre 1983)
    DMMR : np.array of floats of the lengh of the thinned samples
        er^2-3*delta
    er : np.array of floats of the lengh of the thinned samples
        rer parameter of the second fundamental model for resonance (see GRSW paper 1 - Leleu et al 2026) also called resonant eccentricities
    Der : np.array of floats of the lengh of the thinned samples
        Der parameter of the second fundamental model for resonance (see GRSW paper 1 - Leleu et al 2026) 
    es : np.array of floats of the lengh of the thinned samples
        secular eccentricities (see GRSW paper 1 - Leleu et al 2026) 
    IR : np.array of floats of the lengh of the thinned samples
        position of the sample in the second fundamental model:
            2 : upper libration
            1 : formally resonant
            0 : upper circulation
            -1 : lower circulation
            -2 : lower libration
    p : int
        pair attributed to the p+1:p MMR
    
    """

    P1=np.mean(samples['period_days_'+str(pair[0])].values)
    P2=np.mean(samples['period_days_'+str(pair[1])].values)

    if isinstance(samples, pd.DataFrame):
            samples=np.vstack([samples[col] for col in samples.columns])
    

    if p==None:
        dist,p,sym=dist_Dpp1p(P2/P1)
    
    [X, Y, X2, Y2, delta]=samples2SFM( samples[:,::sampling], pair, p)
    [sig, Sig, sig2, Sig2, x1, x2, IR] = samples2usefull( samples[:,::sampling], pair,p)
    
    es=np.sqrt(2*Sig2)

    Der=np.zeros(delta.size)
    er=np.zeros(delta.size)


    for kl in range(delta.size):
        [fixp,o1,o2]=topology_light(delta[kl])
        Vx=np.array([x1[kl], x2[kl]])
        Idmaxx=abs(Vx).argmax()
        Idxlib=abs(Vx[Idmaxx]-np.array([fixp,o1,o2])).argmin()
        Der[kl]=Vx[Idmaxx]-np.array([fixp,o1,o2])[Idxlib]
        er[kl]=Vx[Idmaxx]



    DMMR=er**2-3*delta


    return delta, DMMR, er, Der, es, IR, p




