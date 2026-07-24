import rebound 
import numpy as np
from resonantstate.constants import *
def get_rebound_sim_coplanar(df_samples,Idp,Id_sample):
    """
    Get a rebound simulation for one of the samples of the posterior

    Parameters
    ----------
    df_samples : pandas dataframe
        samples of dynamical analysis
    Idp : list or tuple of int
        planet indicies
    Id_sample : int
        identifiant of the sample to consider

    Returns
    -------
    sim: rebound.Simulation
        Simulation object initialized based on posterior sample
    """    
    m00=df_samples['mass_star_m_sun'].values[Id_sample]
    

    sim = rebound.Simulation()
    sim.G = 4. * pi * pi

    sim.add(m=m00)

    #Create the system
    for k in Idp:
        sim.add(l=df_samples['mean_longitude_deg_'+str(k)].values[Id_sample]*pi/180,
                P=df_samples['period_days_'+str(k)].values[Id_sample]/365.25,
                k=df_samples['k_'+str(k)].values[Id_sample],
                h=df_samples['h_'+str(k)].values[Id_sample],
                m=df_samples['mass_planet_star_ratio_'+str(k)].values[Id_sample]*m00)

    sim.move_to_com()
    return sim

def rebound_sim_coplanar(df_samples,Idp,Id_sample,duration,Noutputs=1000):
    """n-body integration of one of the samples of the posterior

    Parameters
    ----------
    df_samples : pandas dataframe
        samples of dynamical analysis
    Idp : list or tuple of int
        planet indicies
    Id_sample : int
        identifiant of the sample to consider
    duration : _type_
        duration of the integration in years
    Noutputs : int, optional
        number of outputs

    Returns
    -------
    times : numpy array of floats
        times of the outputs
    la : numpy array of floats of dim len(Idp)*Noutputs
        mean longitude (deg)
    a : numpy array of floats of dim len(Idp)*Noutputs
        semi-major axis (au)
    P : numpy array of floats of dim len(Idp)*Noutputs
        period (day)
    e : numpy array of floats of dim len(Idp)*Noutputs
        eccentricity
    pomega : numpy array of floats of dim len(Idp)*Noutputs
        longitude of periastron

    """
    nb_planets=len(Idp)
    sim = get_rebound_sim_coplanar(df_samples,Idp,Id_sample)
    times = np.linspace(0.,duration, Noutputs)
    la = np.zeros((nb_planets,Noutputs))
    e = np.zeros((nb_planets,Noutputs))
    a = np.zeros((nb_planets,Noutputs))
    P = np.zeros((nb_planets,Noutputs))
    pomega = np.zeros((nb_planets,Noutputs))
    sim.integrator = "ias15" # IAS15 is the default integrator, so we actually don't need this line
    ps = sim.particles       # ps is now an array of pointers and will change as the simulation runs

    #run the integration
    for i,time in enumerate(times):
        sim.integrate(time)
        
        for j in range(nb_planets):
            la[j][i] = ps[j+1].l   # This stores the data which allows us to plot it later
            e[j][i] = ps[j+1].e
            a[j][i] = ps[j+1].a
            P[j][i] = ps[j+1].P*365.25
            pomega[j][i] = ps[j+1].pomega
            if ps[j+1].x>1 or ps[j+1].y>1:
                break
    
    return times,la,a,P,e,pomega
