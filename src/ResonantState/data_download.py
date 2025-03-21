import requests
import pandas as pd

def get_metadata_observations():
    URL_metadata_parquet='https://dace.unige.ch/downloads/resonant_state/observations/tables/metadata_table.parquet'
    metadata_parquet = requests.get(URL_metadata_parquet,verify = False) #check stream option , stream=True
    return pd.read_parquet(metadata_parquet)