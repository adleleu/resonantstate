
import ttvfast
import numpy as np
pi=np.pi

def forward_modeling(sample,nb_planets,t_ini, dt, t_end,input_flag=1):
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
    
    results = ttvfast.ttvfast(planets,stellar_mass,t_ini, dt, t_end,input_flag=input_flag)

    resu=np.array(results['positions'])
    
    Id_rel=np.where((resu[2])>-2)[0]
    Id_planetes=resu[0][Id_rel]
    tt=resu[2][Id_rel]

    transits=[]
    for k in range(nb_planets):
        transits.append(tt[Id_planetes==k])

    return transits