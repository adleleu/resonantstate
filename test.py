import numpy as np
import pandas as pd
import requests
import io
import matplotlib.pyplot as plt
URL_metadata_csv='https://dace.unige.ch/downloads/resonant_state/observations/tables/metadata_table.csv'
file = requests.get(URL_metadata_csv, verify=False)
df = pd.read_csv(io.StringIO(file.text))

Kepl51=df[df['star_name']=='Kepler-51'].drop_duplicates(subset=['url_of_samples'])




#contacts=np.unique(Kepl51['contact_mail'].values)
#system_url=np.unique(Kepl51['url_of_samples'].values)


for index, row in Kepl51.iterrows():
    file = requests.get(row['url_of_samples'], verify=False)
    df = pd.read_csv(io.StringIO(file.text))
    for p in range(1,row['nb_planets']+1):
        e=np.sqrt(df['k_j_'+str(p)].values**2+df['h_j_'+str(p)].values**2)
        m=df['mass_planet_star_ratio_'+str(p)].values
        P=df['period_days_'+str(p)].values
        plt.figure(p)
        plt.scatter(m,e,label=row['star_name']+'_'+str(row['analysis_id']),alpha=.1)
        plt.legend()
