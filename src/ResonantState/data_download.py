import logging
import pandas as pd
import requests
import io
from pathlib import Path
from src.ResonantState.table_entries_catalogs import MetaDataTableEntries as mte
from src.ResonantState.table_entries_catalogs import ObservationsSamplesDictionnary as osd
from src.ResonantState.table_entries_catalogs import URLS as urls
import json


log_fmt = '%(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt)
logger = logging.getLogger(__name__)        
    
# Disable warnings due to the use of the requests library with unverified HTTPS requests
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


            
def get_metadata_observations():
    """This method retrieves the metadata table for the observations and return it in a pandas dataframe

    Returns
    -------
    pd.DataFrame
        The observations metadata table in a pandas dataframe
    
    """
    url = urls.OBSERVATIONS_METADATA_TABLE.value

    file = requests.get(url, verify=False)
    if not file.status_code == 200: # check if the request was successful
        raise Exception(f"URL {url} responded with status code: {file.status_code}")
    dataframe = pd.read_parquet(io.BytesIO(file.content))
    
    return dataframe


def get_metadata_summary_simulations():
    """This method retrieves the metadata table for the simulations and return it in a pandas dataframe
    
    Returns
    -------
    pd.DataFrame
        The simulations metadata table in a pandas dataframe
    
    """
    url = urls.SIMULATIONS_SUMMARY_METADATA_TABLE.value
  
    
    file = requests.get(url, verify=False)
    if not file.status_code == 200: # check if the request was successful
        raise Exception(f"URL {url} responded with status code: {file.status_code}")
    dataframe = pd.read_parquet(io.BytesIO(file.content))
        
    return dataframe  



def _get_metadata_complet_simulations(local_metadata_path=None):
    """This method retrieves the metadata table for the simulations and return it in a pandas dataframe
    
    Returns
    -------
    pd.DataFrame
        The simulations metadata table in a pandas dataframe
    
    """
    url = urls.SIMULATIONS_METADATA_TABLE.value
    if local_metadata_path:
        dataframe = pd.read_csv(local_metadata_path)
            
    else:
        file = requests.get(url, verify=False)
        if not file.status_code == 200: # check if the request was successful
            raise Exception(f"URL {url} responded with status code: {file.status_code}")
        dataframe = pd.read_parquet(io.BytesIO(file.content))
        
    return dataframe        



def download_observation_samples(dataframe, download_destination=None):
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
    
    
    unique_samples_urls = []
    metadata_urls = []
    for row in dataframe.iterrows():
        index_row, planet_metadata = row
        
        url_of_sample   = planet_metadata[mte.URL_OF_SAMPLES_PARQUET.value]
        url_of_metadata = planet_metadata[mte.METADATA_FILE.value]
        
        if url_of_sample not in unique_samples_urls:
            
            unique_samples_urls.append(url_of_sample)
            metadata_urls.append(url_of_metadata)
           
    
    
    return_samples = []
    for url_index in range(len(unique_samples_urls)):
        sample_url = unique_samples_urls[url_index]
        metadata_url = metadata_urls[url_index]
        
        file_sample = requests.get(sample_url, verify=False)
        dataframe_sample = pd.read_parquet(io.BytesIO(file_sample.content))
        file_metadata = requests.get(metadata_url, verify=False)
        metadata = json.loads(file_metadata.text)
    
    
        name = f"{metadata[mte.STAR_NAME.value]}_{metadata[mte.ANALYSIS_ID.value]}"
        
        
        samples_dict = {
            osd.SAMPLE_NAME.value:      name,
            osd.PLANETS_LIST.value:     metadata[mte.PLANET_LIST.value],
            osd.CODE.value:            metadata[mte.CODE_USED.value],
            osd.BIBTEX.value:          metadata[mte.BIBTEX.value],
            osd.CONTACT_EMAIL.value:   metadata[mte.CONTACT_EMAIL.value],
            osd.SAMPLE.value:          dataframe_sample
        }
        return_samples.append(samples_dict)
    
        if download_destination is not None:
            file_name = download_destination / f"{name}.parquet"
            dataframe.to_parquet(file_name)
        
    return return_samples



# def download_simulation_samples(dataframe, download_destination=None, local=False):
#     """This method retrieves the samples tables form the systems in the given dataframe, and returns a dictionnary of pandas dataframes with the samples. If a download destination is given, it saves the samples in the given directory.

#     Parameters
#     ----------
#     dataframe : pd.DataFrame
#         The dataframe containing the metadata of the planets samples that we want to download
#     download_destination : str or None
#         If given, save the samples in the given directory
    
#     Returns
#     -------
#     dict
#         A dictionnary of pandas dataframes with the samples names as keys and the samples dataframes as values
#     """
    
#     if isinstance(dataframe, pd.Series):
#         raise ValueError("The dataframe should be a pandas DataFrame, not a pandas Series")
        
#     if not isinstance(dataframe, pd.DataFrame):
#         raise ValueError("The dataframe should be a pandas DataFrame")
    
    
#     if download_destination is not None:
#         if not isinstance(download_destination, Path):
#             download_destination = Path(download_destination)
#         download_destination.mkdir(parents=True, exist_ok=True)
    
    
#     unique_samples_urls = []
#     metadata_urls = []
#     for row in dataframe.iterrows():
#         index_row, metadata_line = row
        
#         url_of_sample   = metadata_line.get(mte.URL_OF_SAMPLES_PARQUET.value, None)
#         url_of_metadata = metadata_line.get(mte.METADATA_FILE.value, None)

#         if url_of_sample not in unique_samples_urls:
#             unique_samples_urls.append(url_of_sample)
#             metadata_urls.append(url_of_metadata)
        

    
#     return_samples = []
#     for url_index in range(len(unique_samples_urls)):
#         sample_url = unique_samples_urls[url_index]
#         metadata_url = metadata_urls[url_index]
        
#         if not local:
#             file_sample = requests.get(sample_url, verify=False)
#             dataframe_sample = pd.read_parquet(io.BytesIO(file_sample.content))
#             file_metadata = requests.get(metadata_url, verify=False)
#             metadata = json.loads(file_metadata.text)
#         else:
#             file_sample = Path(sample_url)
#             relative_path_sample =file_sample.relative_to(file_sample.parent.parent.parent.parent)
#             logger.info(f"Relative path sample: {relative_path_sample}")
#             return 0    
        
#         if metadata.get(mte.STAR_NAME.value, False) and metadata.get(mte.ANALYSIS_ID.value, False):
#             name = f"{metadata[mte.STAR_NAME.value]}_{metadata[mte.ANALYSIS_ID.value]}"
#         elif metadata.get(mte.SYSTEM_ID.value, False) and metadata.get(mte.RUN_ID.value, False):
#             name = f"{metadata[mte.SYSTEM_ID.value]}_{metadata[mte.RUN_ID.value]}"
#         else:
#             logger.error(f"Metadata does not contain the required keys for activity A or B")
#             name = Path(sample_url).stem
        
#         samples_dict = {
#             "name":  name,
#             "planets": f"{metadata[mte.PLANET_LIST.value]}",
#             "code": metadata[mte.CODE_USED.value],
#             "bibtex": metadata[mte.BIBTEX.value],
#             "samples": dataframe_sample
#         }
#         return_samples.append(samples_dict)
    
#         if download_destination is not None:
#             file_name = download_destination / name
#             dataframe.to_parquet(file_name)
        
#     return return_samples


# def download_simulation_run_tables(dataframe, download_destination=None, local=False):
#     """This method retrieves the samples tables form the systems in the given dataframe, and returns a dictionnary of pandas dataframes with the samples. If a download destination is given, it saves the samples in the given directory.

#     Parameters
#     ----------
#     dataframe : pd.DataFrame
#         The dataframe containing the metadata of the planets samples that we want to download
#     download_destination : str or None
#         If given, save the samples in the given directory
    
#     Returns
#     -------
#     dict
#         A dictionnary of pandas dataframes with the samples names as keys and the samples dataframes as values
#     """
    
#     if isinstance(dataframe, pd.Series):
#         raise ValueError("The dataframe should be a pandas DataFrame, not a pandas Series")
        
#     if not isinstance(dataframe, pd.DataFrame):
#         raise ValueError("The dataframe should be a pandas DataFrame")
    
    
#     if download_destination is not None:
#         if not isinstance(download_destination, Path):
#             download_destination = Path(download_destination)
#         download_destination.mkdir(parents=True, exist_ok=True)
    
    
#     unique_samples_urls = []
#     metadata_urls = []
#     for row in dataframe.iterrows():
#         index_row, metadata_line = row
        
#         url_of_sample   = metadata_line.get(mte.URL_OF_SAMPLES_PARQUET.value, None)
#         url_of_metadata = metadata_line.get(mte.METADATA_FILE.value, None)

#         if url_of_sample not in unique_samples_urls:
#             unique_samples_urls.append(url_of_sample)
#             metadata_urls.append(url_of_metadata)
        

    
#     return_samples = []
#     for url_index in range(len(unique_samples_urls)):
#         sample_url = unique_samples_urls[url_index]
#         metadata_url = metadata_urls[url_index]
        
#         if not local:
#             file_sample = requests.get(sample_url, verify=False)
#             dataframe_sample = pd.read_parquet(io.BytesIO(file_sample.content))
#             file_metadata = requests.get(metadata_url, verify=False)
#             metadata = json.loads(file_metadata.text)
#         else:
#             file_sample = Path(sample_url)
#             relative_path_sample =file_sample.relative_to(file_sample.parent.parent.parent.parent)
#             logger.info(f"Relative path sample: {relative_path_sample}")
#             return 0    
        
#         if metadata.get(mte.STAR_NAME.value, False) and metadata.get(mte.ANALYSIS_ID.value, False):
#             name = f"{metadata[mte.STAR_NAME.value]}_{metadata[mte.ANALYSIS_ID.value]}"
#         elif metadata.get(mte.SYSTEM_ID.value, False) and metadata.get(mte.RUN_ID.value, False):
#             name = f"{metadata[mte.SYSTEM_ID.value]}_{metadata[mte.RUN_ID.value]}"
#         else:
#             logger.error(f"Metadata does not contain the required keys for activity A or B")
#             name = Path(sample_url).stem
        
#         samples_dict = {
#             "name":  name,
#             "planets": f"{metadata[mte.PLANET_LIST.value]}",
#             "code": metadata[mte.CODE_USED.value],
#             "bibtex": metadata[mte.BIBTEX.value],
#             "samples": dataframe_sample
#         }
#         return_samples.append(samples_dict)
    
#         if download_destination is not None:
#             file_name = download_destination / name
#             dataframe.to_parquet(file_name)
        
#     return return_samples