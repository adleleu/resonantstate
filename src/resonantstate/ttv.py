
import ttvfast
import numpy as np
pi=np.pi
AUpdaytomps=149597870700/(24*3600)

def forward_modeling(sample,nb_planets,t_ini, dt, t_end,input_flag=1,rv_times=[]):
    '''
    input_flag: input coordinate system (0, 1 or 2)

    The `input_flag` argument corresponds to:
        0 = Jacobi
        1 = astrocentric elements
        2 = astrocentric cartesian
    '''

    planets=[]
    for k in range(nb_planets):
        Omega=sample[6+k*8]
        long_peri=np.arctan2(sample[4+k*8],sample[3+k*8])*180/pi
        mean_anomaly=np.mod(sample[1+k*8]-long_peri,360)
        eccentricity=np.sqrt(sample[3+k*8]**2+sample[4+k*8]**2) 

        planets.append(
            ttvfast.models.Planet(
                        mass=sample[7+k*8]*sample[8*nb_planets+1],  # M_sun
                        period=sample[2+k*8],     # days
                        eccentricity=eccentricity,
                        inclination=sample[5+k*8],         # degrees
                        longnode=Omega,           # degrees
                        argument=long_peri-Omega,            # degrees
                        mean_anomaly=mean_anomaly,      # degrees
                        )
                )

   
    stellar_mass=sample[8*nb_planets+1]
    
    if len(rv_times)>0:
        results = ttvfast.ttvfast(planets,stellar_mass,t_ini, dt, t_end,input_flag=input_flag,rv_times=rv_times)
    else:
        results = ttvfast.ttvfast(planets,stellar_mass,t_ini, dt, t_end,input_flag=input_flag)

    resu=np.array(results['positions'])
    
    Id_rel=np.where((resu[2])>-2)[0]
    Id_planetes=resu[0][Id_rel]
    tt=resu[2][Id_rel]

    transits=[]
    for k in range(nb_planets):
        transits.append(tt[Id_planetes==k])

    if len(rv_times)>0:
        return transits,np.array(results['rv'])*AUpdaytomps
    else:
        return transits


def fit_lin(y,T_vect,sig_y):
    #linear fit of y by the vectors,  sig_y is the error on y
    nb_col=len(T_vect)
    nb_line=y.size
    Mt=np.zeros((nb_col,nb_line))
    for k in range(len(T_vect)):
        Mt[k][:] = T_vect[k]
    M= Mt.transpose()
    W= np.diag(1/sig_y**2)
    Mresh=M.reshape((nb_line,nb_col))
    WM=np.dot(W,Mresh)
    MtWM=np.dot(Mt,WM)
    MtWMm1=np.linalg.inv(MtWM)
    Wy=np.dot(W,y)
    MtWy=np.dot(Mt,Wy)
    theta=np.dot(MtWMm1,MtWy)
    signal=np.dot(M,theta)
    sum_weight=sum((signal-y)**2/sig_y**2) 
    return sum_weight,theta,signal


def fit_trend(epoch,TTV,VW):
    return fit_lin(TTV,[1.0+epoch*0,epoch],VW)

def get_TTV_detrend(Vt0,per):
    TTV=(Vt0)%per
    TTV=TTV-TTV[0]

    for k in range(TTV.size-1):
        if(TTV[k+1]>(TTV[k]+per/2)):
            TTV[k+1:]=TTV[k+1:]-per
        if(TTV[k+1]<(TTV[k]-per/2)):
            TTV[k+1:]=TTV[k+1:]+per
    sum_weight,theta,signal=fit_trend(Vt0,TTV,TTV*0+1/48)
    TTV=TTV-signal
    return TTV



def TTV_amplitudes(sample,nb_planets,t_ini, dt, t_end,input_flag=1):
    transits=forward_modeling(sample,nb_planets,t_ini, dt, t_end,input_flag=input_flag)
    amplitudes=np.zeros(nb_planets)
    for k in range(nb_planets):
        if transits[k].size>2:
            TTV=get_TTV_detrend(transits[k],sample[2+k*8])
            amplitudes[k]=TTV.max()-TTV.min()
        else:
            amplitudes[k]=-1
    return amplitudes