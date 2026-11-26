from resonantstate.data_download  import get_metadata_simulations,download_simulations_run_table, download_simulations
from resonantstate.analyse_samples import *
from resonantstate.SFMR import *
import numpy as np
from tqdm import tqdm
import requests
from pathlib import Path
import json

pi=np.pi

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



def get_infos(author_name,run_id,simulation_type):

    simulation_metadata_table = get_metadata_simulations()


    simulations = simulation_metadata_table[(simulation_metadata_table["author_name"] == author_name)&(simulation_metadata_table["run_id"]==run_id)&(simulation_metadata_table["simulation_type"]==simulation_type)]

    simulation=simulations.iloc[[0]]

    print(simulation['author_name'].values[0]+'_'+simulation['simulation_type'].values[0]+'_'+str(simulation['run_id'].values[0]))

    table = download_simulations_run_table(simulation)

    readme_url=table.readme_file.values[0]
    metadata_url=table.metadata_file.values[0]

    if (readme_url is not None  ) and (readme_url != 'None'):
        file_readme = requests.get(readme_url, verify=False)
        if not file_readme.ok: 
            raise Exception(f"URL {readme_url} responded with status code: {file_readme.status_code}")
        readme = file_readme.text
    else:
        readme = None


    file_metadata = requests.get(metadata_url, verify=False)
    if not file_metadata.ok: 
        raise Exception(f"URL {metadata_url} responded with status code: {file_metadata.status_code}")
    metadata = json.loads(file_metadata.text)
    
        
    return readme,metadata

def process_sim(author_name,run_id,simulation_type,Pmax=50,take_all=False,line_to_take=-2):

    simulation_metadata_table = get_metadata_simulations()


    simulations = simulation_metadata_table[(simulation_metadata_table["author_name"] == author_name)&(simulation_metadata_table["run_id"]==run_id)&(simulation_metadata_table["simulation_type"]==simulation_type)]


    for sim in range(simulations.shape[0]):
        simulation=simulations.iloc[[sim]]

        print('starting '+simulation['author_name'].values[0]+'_'+simulation['simulation_type'].values[0]+'_'+str(simulation['run_id'].values[0]))
            

        table = download_simulations_run_table(simulation)

        T_keep=[]
        
        for sim_Id in tqdm(np.unique(table['simulation_id'].values)):
        
            run_1 = table[table["simulation_id"] == sim_Id]
            
            TP=run_1['period_days'].values
            TR=run_1['radius_planet_r_earth'].values

            if take_all:
                Id=np.arange(TP.size)
            else:
                cond=np.logical_and(TP<Pmax , TR>1)
                cond2=np.logical_and(cond , TP>1.5)
                Id=np.where(cond2)[0]

            TPs=TP[Id] #list of period of selected planets
            Ids=TP[Id].argsort() # Ids sorting the period of selected planets
            TPs.sort()


            #The first check is done on the metadata table values, to see if it is worth loading the simluation
            if TPs.size>1:
                
                if (TPs[1:]/TPs[:-1]).min()<2.5: #we only load simluations for which a pair of planet is has a period ratio less than 2.1
                    
                    run_load_ok=False
                    try:
                        run1_data = download_simulations(run_1)
                        run_load_ok=True
                    except:
                        print('failed run ',str(run_1))
                    
                    if run_load_ok:
                        
                        df_evol_full=run1_data[0]['simulation']
                        df_evol=run1_data[0]['simulation'].iloc[[line_to_take]]
                        
                        TP,TR,TM,nb_planets,Idsort=get_PRM(df_evol,run1_data[0])
                        
                        if take_all:
                            Id=np.arange(TP.size)
                        else:
                            cond=np.logical_and(TP<Pmax , TR>1)
                            cond2=np.logical_and(cond , TP>1.5)
                            Id=np.where(cond2)[0] #indice of the planets to keep from the simluation

                        Ids=TP[Id].argsort() #how they are sorted
                        Id=Id[Ids]

                        evol_array=np.vstack([df_evol[col] for col in df_evol.columns])

                        evol=evol_array[:nb_planets*8+3,-1]

                        for k in range(Id.size-1):
                            Id1=Id[k]
                            Id2=Id[k+1]
                            pair=(Id1,Id2)

                            dist,p,sym=dist_Dpp1p(TP[Id2]/TP[Id1])
                            
                            if sym:

                                
                                SFM = get_SFMR_quantities(df_evol,pair)
                                

                                # count number of meighbouring planets
                                cond=np.logical_and(TP<3.5*TP[Id2],TP>TP[Id1]/3.5)
                                Idnp=np.where(cond)[0]
                                
                        
                                #find the period of the innermost planet of the (near-)resonant group
                                Idk=np.where(TP[Idsort]==TP[Id1])[0][0]
                                Idmin=Idsort[Idk]
                                Pmin=TP[Idmin]
                                while Idk>=0 and (TP[Idsort[Idk+1]]/TP[Idsort[Idk]])<2.1:
                                    Pmin=TP[Idsort[Idk]]
                                    Idk-=1 


                                dic_pair={}
                                
                                dic_pair['simulation_type']=simulation['simulation_type'].values[0]
                                dic_pair['run_id']=simulation['run_id'].values[0]
                                dic_pair['author_name']=simulation['author_name'].values[0]
                                dic_pair['simulation_id']=sim_Id
                                dic_pair['P1']=TP[Id1]
                                dic_pair['P2']=TP[Id2]
                                dic_pair['R1']=TR[Id1]
                                dic_pair['R2']=TR[Id2]
                                dic_pair['M1']=TM[Id1]
                                dic_pair['M2']=TM[Id2]
                                dic_pair['er']=SFM.er[0]
                                dic_pair['es']=SFM.es[0]
                                dic_pair['Der']=SFM.Der[0]
                                dic_pair['delta']=SFM.delta[0]
                                dic_pair['DMMR']=SFM.DMMR[0]
                                dic_pair['Pmin']=Pmin
                                dic_pair['Prdist']=dist
                                dic_pair['nb_nearby_planets']=Idnp.size
                                dic_pair['resonantstate']=SFM.RS[0]
                                dic_pair['Id1']=Id1
                                dic_pair['Id2']=Id2
                                dic_pair['p']=p


                                T_keep.append(dic_pair)
                            


        df_keep = pd.DataFrame(T_keep)                            
        df_keep.to_parquet('data/'+simulation['author_name'].values[0]+'_'+simulation['simulation_type'].values[0]+'_'+str(simulation['run_id'].values[0]))
            
