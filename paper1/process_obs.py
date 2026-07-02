import pandas as pd
import numpy as np
from resonantstate.data_download  import get_metadata_observations, download_observations_samples
from resonantstate.analyse_samples import *
from resonantstate.ell2SFM import *
from resonantstate.simulations_resonance_analysis import *
from resonantstate.utils import *
from resonantstate import ttv
from astropy.timeseries import LombScargle
from tqdm import tqdm


def process_obs(case,
                download_destination_path = None,
                ):


    dataframe_observations = get_metadata_observations()
    dataframe_observations.head()
   

    df_keep=pd.read_parquet('data/keep_'+case)


    #We pre allocate the new quantities we compute in this function
    T_additional_columns=['xres','xres_e1','xres_e2',
                        'delta','delta_e1','delta_e2',  
                        'xsec','xsec_e1','xsec_e2',
                        'Prdist','Prdist_e1','Prdist_e2',
                        'ein','ein_e1','ein_e2','eout','eout_e1','eout_e2',
                        'xlib','xlib_e1','xlib_e2',
                        'mutinc','mutinc_e1','mutinc_e2',
                        'ttv_amp_1','ttv_amp_2',
                        'ttv_phase_1','ttv_phase_2',
                        'author_name','prior_e',
                        'nb_planets',
                        'nb_nearby_planets','resonantstate','bibtex']


    for col in T_additional_columns:
        df_keep[col] = None


    #We also store 300 samples of some chosen quantities. this allow properly explore potential orrelations between parameters
    T_samples_columns=['author','reason','Pmin','delta','xres','xsec','xlib','mutinc','rhoin','rhoout']
    df_samples = pd.DataFrame(None, index=range(df_keep.shape[0]*300), columns=T_samples_columns)




    nb_pair=df_keep.shape[0]
    nb_pair_done=0


    print('loop over near resonant pairs')
    for index, row in tqdm(df_keep.iterrows(), total=df_keep.shape[0]):
        
        nb_pair_done+=1

        dataframe_system = dataframe_observations[(dataframe_observations["star_name"]==row['star_name'])&(dataframe_observations["analysis_id"]==row['analysis_id'])]
        
        #get the infos from the analysis. this should be unique if the pre-selection was done correctly
        df_ana = download_observations_samples(dataframe_system, download_destination_path)[0]

        robustness=df_ana['robustness']
        samples=df_ana['samples']
        nb_planets=len(df_ana['planets_list'])
        nbsamples=df_ana['samples'].shape[0]
        Id1=row['Id1']
        Id2=row['Id2']
        pair=(Id1,Id2)


        delta, DMMR, xres, xlib, Sig2, IR , p = get_SFM_quantities(samples,pair)
        
        TP=np.zeros(nb_planets)
        TR=np.zeros(nb_planets)
        for k in range(nb_planets):
            TP[k]=samples['period_days_'+str(k)].mean()
            TR[k]=(samples['radius_planet_star_ratio_'+str(k)]*samples['radius_star_r_sun']/0.00916794).mean()

        cond_group=np.logical_and(TP<3.5*TP[Id2],TP>TP[Id1]/3.5)

        Prdist,p,sym=dist_Dpp1p(samples['period_days_'+str(Id2)].values/samples['period_days_'+str(Id1)].values)

        
        IRv, IRc = np.unique(IR, return_counts=True)

        df_keep.loc[index,'nb_planets']=nb_planets
        df_keep.loc[index,'nb_nearby_planets']=sum(cond_group)

        df_keep.loc[index,['resonantstate']]=IRv[np.argmax(IRc)]
        df_keep.loc[index,['er','er_e1','er_e2']]=np.quantile(xres,(0.5,.16,.84))
        df_keep.loc[index,['delta','delta_e1','delta_e2']]=np.quantile(delta,(0.5,.16,.84))
        df_keep.loc[index,['DMMR','DMMR_e1','DMMR_e2']]=np.quantile(xres**2-3*(delta+1),(0.5,.16,.84))
        df_keep.loc[index,['es','es_e1','es_e2']]=np.quantile(np.sqrt(2*Sig2),(0.5,.16,.84))
        df_keep.loc[index,['Der','Der_e1','Der_e2']]=np.quantile(xlib,(0.5,.16,.84))
        df_keep.loc[index,['ein','ein_e1','ein_e2']]=np.quantile(np.sqrt(samples['h_'+str(Id1)]**2+samples['k_'+str(Id1)]**2),(0.5,.16,.84))
        df_keep.loc[index,['eout','eout_e1','eout_e2']]=np.quantile(np.sqrt(samples['h_'+str(Id2)]**2+samples['k_'+str(Id2)]**2),(0.5,.16,.84))
        df_keep.loc[index,['Prdist','Prdist_e1','Prdist_e2']]= np.quantile(Prdist,(0.5,.16,.84))

        df_keep.loc[index,'bibtex']=df_ana['bibtex']
        df_keep.loc[index,'author_name']=df_ana['author_name']
        df_keep.loc[index,'prior_e']=df_ana['eccentricity_prior']


    
    df_keep.to_parquet('data/'+case+'_rs')
    df_samples.to_parquet('data/'+case+'_rs_samples')