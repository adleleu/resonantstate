import logging
import pandas as pd
import requests
import io
from pathlib import Path

log_fmt = '%(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt)
logger = logging.getLogger(__name__)        
    
# Disable warnings due to the use of the requests library with unverified HTTPS requests
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


URL_OBSERVATIONS_METADATA_TABLE = "https://dace.unige.ch/downloads/resonant_state/observations/tables/metadata_table.csv"
URL_SIMULATIONS_METADATA_TABLE = "https://dace.unige.ch/downloads/resonant_state/simulations/tables/metadata_table.csv"

KEY_URL_OF_SAMPLES                  = "url_of_samples"                      # str

            
def get_metadata_observations(local_metadata_path=None):
    """This method retrieves the metadata table for the observations and return it in a pandas dataframe

    Returns
    -------
    pd.DataFrame
        The observations metadata table in a pandas dataframe
    
    """
    
    if local_metadata_path:
        dataframe = pd.read_csv(local_metadata_path)
            
    else:
        file = requests.get(URL_OBSERVATIONS_METADATA_TABLE, verify=False)
        if not file.status_code == 200: # check if the request was successful
            raise Exception(f"URL {URL_OBSERVATIONS_METADATA_TABLE} responded with status code: {file.status_code}")
        dataframe = pd.read_csv(io.StringIO(file.text))
        
    return dataframe
        

def get_metadata_simulations(local_metadata_path=None):
    """This method retrieves the metadata table for the simulations and return it in a pandas dataframe
    
    Returns
    -------
    pd.DataFrame
        The simulations metadata table in a pandas dataframe
    
    """
    
    if local_metadata_path:
        dataframe = pd.read_csv(local_metadata_path)
            
    else:
        file = requests.get(URL_SIMULATIONS_METADATA_TABLE, verify=False)
        if not file.status_code == 200: # check if the request was successful
            raise Exception(f"URL {URL_SIMULATIONS_METADATA_TABLE} responded with status code: {file.status_code}")
        dataframe = pd.read_csv(io.StringIO(file.text))
        
    return dataframe        



def download_samples(dataframe, download_destination=None):
    """This method retrieves the samples tables form the systems in the given dataframe, and returns a dictionnary of pandas dataframes with the samples. If a download destination is given, it saves the samples in the given directory.

    Parameters
    ----------
    dataframe : pd.DataFrame
        The dataframe containing the metadata of the planets samples that we want to download
    download_destination : str or None
        If given, save the samples in the given directory
    
    Returns
    -------
    dict
        A dictionnary of pandas dataframes with the samples names as keys and the samples dataframes as values
    """
    
    if isinstance(dataframe, pd.Series):
        raise ValueError("The dataframe should be a pandas DataFrame, not a pandas Series")
        
    if not isinstance(dataframe, pd.DataFrame):
        raise ValueError("The dataframe should be a pandas DataFrame")
    
    
    if download_destination is not None:
        if not isinstance(download_destination, Path):
            download_destination = Path(download_destination)
        download_destination.mkdir(parents=True, exist_ok=True)
    
    
    unique_urls = []
    for row in dataframe.iterrows():
        index_row, planet_metadata = row
        url_of_sample = planet_metadata[KEY_URL_OF_SAMPLES]
        
        if url_of_sample not in unique_urls:
            unique_urls.append(url_of_sample)
    
    
    dataframes = {}
    for url in unique_urls:
        file = requests.get(url, verify=False)
        
        dataframe = pd.read_parquet(io.BytesIO(file.content))
        
        if download_destination is not None:
            file_name = download_destination / Path(url).name
            dataframe.to_parquet(file_name)
        
        key_name = Path(url).name
        for potential_endfile_names in [".csv", ".parquet", "_samples", "_sample", "_evol", "_evols"]:
            key_name = key_name.replace(potential_endfile_names, "")
        
        dataframes[key_name] = dataframe
    
    return dataframes