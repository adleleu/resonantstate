import numpy as np
from resonantstate.data_download  import get_metadata_observations, download_observations_samples
from resonantstate.analyse_samples import *
from resonantstate.ell2SFM import *
from resonantstate.simulations_resonance_analysis import *
from resonantstate.utils import *
from tqdm import tqdm


    
def choose_obs(case,
               flag_lowe=False, #take the low eccentricity prior when available
                Pmax=50, #max period of the outer planet of the pair
                delta_max=200, # max median of the delta parameter to avoid solution overly degenerated.
                Pmin_maxvalue=40,  #max period of the inner planet of the resonant group
                xintersect_spread_max=1, # max std of the xintersect parameter for class 0
                prec_max_rob=.07, #precision threshold to consider a pair robust, based on Hadden & Lithwick 2017 and Leleu et al 2023, we consider robust solutions for which the error on the mass is better than 7%
                download_destination_path=None,
                local_path=None):

    dataframe_observations = get_metadata_observations(local_path=None)
    dataframe_observations.head()

    #cleaning up the database
    #analysed individually and found unconstrained
    dataframe_observations=dataframe_observations[dataframe_observations["star_name"]!='Kepler-444'] 
    #synthetic dataset
    dataframe_observations=dataframe_observations[dataframe_observations["star_name"]!='SOL'] 

    #issue with these
    dataframe_observations=dataframe_observations[dataframe_observations["author_name"]!='Ragozzine'] 
    dataframe_observations=dataframe_observations[dataframe_observations["author_name"]!='MacDonald'] 

    #The author didn't attribute default solution to any analysis
    dataframe_observations.loc[(dataframe_observations['author_name']=='Almenara') , 'default'] = 1  
    dataframe_observations.loc[(dataframe_observations['author_name']=='Almenara')&(dataframe_observations['star_name']=='Kepler-138')&(dataframe_observations['analysis_id']==0), 'default']=0

    #get the list of unique stars to compare the solution.
    dataframe_observations.append(dataframe_observations)
    star_names=np.unique(dataframe_observations['star_name'].values)

    nb_star=star_names.size
    nb_star_done=0

    #downsampling to speed up
    sampling=10

    T_keep=[]

    print('looping over stars')
    for star_name in tqdm(star_names):

        nb_star_done+=1
        
        dataframe_system = dataframe_observations[(dataframe_observations["star_name"]==star_name)]
        if flag_lowe:
            favored=dataframe_system['eccentricity_prior']=='log-uniform'
            target_default=0
        else:
            favored=dataframe_system['robustness']>=0
            target_default=1

        if sum(favored)>0: 
        #favor analysis that checked robustness
            dataframe_system_favored=dataframe_system[favored]
            if sum(dataframe_system_favored['methods'].str.contains('photo-dynamics'))>0 or sum(dataframe_system_favored['methods'].str.contains('photodynamics'))>0:
            #favor photodynamical analysis
                dataframe_system_dl=dataframe_system_favored[(dataframe_system_favored['methods'].str.contains('photo-dynamics')|dataframe_system_favored['methods'].str.contains('photodynamics'))&(dataframe_system_favored['default']==target_default)]
            else:
                dataframe_system_dl=dataframe_system_favored[((dataframe_system_favored['methods'].str.contains('extracted times'))|(dataframe_system_favored['methods'].str.contains('TTV')))&(dataframe_system_favored['default']==target_default)]
        else:
            if sum(dataframe_system['methods'].str.contains('photo-dynamics'))>0 or sum(dataframe_system['methods'].str.contains('photodynamics'))>0:
            #favor photodynamical analysis
                dataframe_system_dl=dataframe_system[(dataframe_system['methods'].str.contains('photo-dynamics')|dataframe_system['methods'].str.contains('photodynamics'))&(dataframe_system['default']==1)]
            else:
                dataframe_system_dl=dataframe_system[((dataframe_system['methods'].str.contains('extracted times'))|(dataframe_system['methods'].str.contains('TTV')))&(dataframe_system['default']==1)]

        if dataframe_system_dl['author_name'].values.size==0:
            input('no'+star_name+' nb favored='+str(sum(favored)) )
        
   
        df_list = download_observations_samples(dataframe_system_dl, local_path=local_path,download_destination=download_destination_path)
       
        TP1,TP2,T_star_name=[],[],[] #store periods and star name to check duplicates

        for l in range(len(df_list)):

            df_ana=df_list[l]
            
            robustness=df_ana['robustness']
            samples=df_ana['samples']


            TP,TR,TM,nb_planets,Idsort=get_PRM(samples,df_ana)
            
            for k in range(nb_planets-1):
                Id1=Idsort[k]
                Id2=Idsort[k+1]
                pair=(Id1,Id2)


                if TP[Id2]<Pmax and TP[Id1]<Pmax  :

                    samples_array=np.vstack([samples[col] for col in samples.columns])
                    near_res=get_nearest_resonance(TP[Id2]/TP[Id1], second_order = False, kmax=12, difference_order = 0.2)

                    dist,p,selection_criterion=dist_Dpp1p(TP[Id2]/TP[Id1])
                    
                    if selection_criterion:

                        [X, Y, X2, Y2, delta]=samples2SFM( samples_array[:,::sampling], pair, near_res[0])
                        [sig, Sig, sig2, Sig2, x1, x2, IR] = samples2usefull( samples_array[:,::sampling], pair,near_res[0])

                        #find the max value of the intersection with the delta,x_res plane of the SFM
                        xintersect=np.maximum(x1,x2)
                        
                        #find the period of the innermost planet of the 'resonant' group of planets.
                        Idk=np.where(TP[Idsort]==TP[Id1])[0][0]
                        Idmin=Idsort[Idk]
                        Pmin=TP[Idmin]
                        while Idk>=0 and (TP[Idsort[Idk+1]]/TP[Idsort[Idk]])<2.1:
                            Pmin=TP[Idsort[Idk]]
                            Idk-=1 

                        #based on Hadden & Lithwick 2017 and Leleu et al 2023, we consider robust solutions for which the error on the mass is better than prec_max_rob
                        mtotmstar=samples['mass_planet_star_ratio_'+str(Id1)].values+samples['mass_planet_star_ratio_'+str(Id2)].values

                        #criteria for solutions to consider
                        qdel=np.quantile(delta,(.16,.5,.84))
                        consider=abs(qdel[1])<delta_max and Pmin<Pmin_maxvalue


                        if consider:
                            #check if the pair is new
                            Id_already=np.where((abs(TP[Id1]-np.array(TP1))<.2) & (abs(TP[Id2]-np.array(TP2))<.2) & (np.array(T_star_name)==df_ana['star_name']))[0]

                            dic_sol={}
                            dic_sol['P1']=TP[Id1]
                            dic_sol['P2']=TP[Id2]
                            dic_sol['R1']=TR[Id1]
                            dic_sol['R2']=TR[Id2]
                            dic_sol['M1']=TM[Id1]
                            dic_sol['M2']=TM[Id2]
                            dic_sol['Pmin']=Pmin
                            dic_sol['Id1']=Id1
                            dic_sol['Id2']=Id2
                            dic_sol['analysis_id']=df_ana['analysis_id']
                            dic_sol['star_name']=df_ana['star_name']

                            if Id_already.size==0:
                                #if the pair is new
                                if robustness[Id1]==1 and robustness[Id2]==1 :
                                    dic_sol['reason']=1
                                elif mtotmstar.std()/mtotmstar.mean() < prec_max_rob: 
                                    dic_sol['reason']=.5
                                elif xintersect.std()<xintersect_spread_max:
                                    dic_sol['reason']=0
                                elif df_ana['eccentricity_prior']=='log-uniform': #favor low eccentricity solution
                                    dic_sol['reason']=-1
                                else :
                                    dic_sol['reason']=-2

                                TP1.append(TP[Id1])
                                TP2.append(TP[Id2])
                                T_star_name.append(df_ana['star_name'])

                                T_keep.append(dic_sol)

                            else:
                                #if the pair is not new
                                print('already have '+df_ana['star_name']+' '+str(TP[Id1])+' '+str(TP[Id2])+' '+str(Id_already.size))

                                new_reason=-2
                                if robustness[Id1]==1 and robustness[Id2]==1 and T_keep[Id_already[0]]['reason']<1:
                                    new_reason=1
                                elif mtotmstar.std()/mtotmstar.mean() < prec_max_rob and T_keep[Id_already[0]]['reason']<0.5:
                                    new_reason=.5
                                elif  mtotmstar.std()/mtotmstar.mean()<.1 and T_keep[Id_already[0]]['reason']<0:
                                    new_reason=0         
                                elif df_ana['eccentricity_prior']=='log-uniform' and T_keep[Id_already[0]]['reason']<-1:
                                    new_reason=-1      

                                if new_reason>-2:
                                    dic_sol['reason']=new_reason
                                    T_keep[Id_already[0]]=dic_sol
                                    print('updated')
                    

    df_keep = pd.DataFrame(T_keep)                            
    df_keep.to_parquet('data/keep_'+case)