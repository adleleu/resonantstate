import pandas as pd
import requests
import io
from pathlib import Path
from resonantstate.table_entries_catalogs import MetaDataTableEntriesObservations as mteObs
from resonantstate.table_entries_catalogs import MetaDataTableEntriesSimulations as mteSim
from resonantstate.table_entries_catalogs import MetaDataTableEntriesSimulationsSummary as mteSimSum
from resonantstate.table_entries_catalogs import AuthorsMetaDataTableEntries as mteAuthors
from resonantstate.table_entries_catalogs import ObservationsSamplesDictionnary as obsDict
from resonantstate.table_entries_catalogs import SimulationsDictionnary as simDict
from resonantstate.table_entries_catalogs import URLS as urls
import json


          
# Disable warnings due to the use of the requests library with unverified HTTPS requests
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning) 


# #######################################################################################################
# OBSERVATIONS FUNCTION #################################################################################
# #######################################################################################################

def get_metadata_observations():
    """This method retrieves the metadata table for the observations and returns it in a pandas dataframe

    Returns
    -------
    pd.DataFrame
        The observations metadata table in a pandas dataframe
    
    """
    
    url = urls.OBSERVATIONS_METADATA_TABLE.value

    # we query the table from DACE and make sure the status code is sucessful (http 200 OK = sucessful http request)
    file = requests.get(url, verify=False)
    if not file.ok: # check if the request was successful
        raise Exception(f"URL {url} responded with status code: {file.status_code}")
    dataframe = pd.read_parquet(io.BytesIO(file.content)) # translate the content of the parquet file into a pandas dataframe
    
    return dataframe


def download_observations_samples(dataframe, download_destination=None):
    """This method retrieves the samples for the systems in the given dataframe, and returns a dictionnary that contains informations about the samples and the samples themselves. 
    
    If a download destination is given, it saves the samples and the coresponding metadata in the given directory.

    Parameters
    ----------
    dataframe : pd.DataFrame
        The dataframe containing the metadata of the planets samples that we want to download
    download_destination : str or None
        If given, save each sample and metadata in the given directory 
    
    Returns
    -------
    List[Dict]
        A list of dictionnaries with the samples informations and the samples themselves
        keys of the dictionnary are:
        
        - samples_name:     the name of the samples (str)
        - planets_list:     the list of planets in the sample (list of str)
        - samples:          the samples themselves (pandas.DataFrame) 
        - readme:           the readme of the samples (str|None)
        - additional_infos: the additional infos of the samples (str|None)
        - author_name:      the name of the author of the samples (str)
        - star_name:        the name of the star of the sampled system (str)
        - analysis_id:      the analysis id of the samples (str)
        - contact_email:    the email of the contact person for the sample (str)
        - default:          the default value of the samples (str)
        - robustness:       the robustness of the samples (str)
        - initial_condition_date_bjd: the initial condition date of the samples (str)
        - nb_planets:       the number of planets in the sample (int)
        - gaia_id:          the gaia id of the sample (str)
        - mass_priors:      the mass priors of the sample (str)
        - eccentricity_priors:  the eccentricity priors of the sample (str)
        - transit_definition:   the transit definition of the sample (str|None)
        - methods:          the methods used to generate the sample (str)
        - instruments:      the instruments used to generate the sample (str)
        - bibtex:           the bibtex reference of the sample (str)
        - code:             the code used to generate the sample (str)
        - other_remarks:    the other remarks of the sample (str)
        
    """
    # common mistake when manipulating pandas dataframe and does not work with the package so it has its own check
    if isinstance(dataframe, pd.Series):
        raise ValueError("The dataframe should be a pandas DataFrame, not a pandas Series")
        
    if not isinstance(dataframe, pd.DataFrame):
        raise ValueError("The dataframe should be a pandas DataFrame")
    
    # if the download destination is given, we create the directory if it does not exist
    if download_destination is not None:
        if not isinstance(download_destination, Path):
            download_destination = Path(download_destination)
        download_destination.mkdir(parents=True, exist_ok=True)
    
    # we create a list of unique urls for the samples and the metadata
    unique_samples_urls = []
    metadata_urls = []
    readme_urls = []
    additional_infos_urls = []
    planet_metadatas = []

    for row in dataframe.iterrows():
        index_row, planet_metadata = row
        
        url_of_sample           = planet_metadata[mteObs.URL_OF_SAMPLES.value]
        url_of_metadata         = planet_metadata[mteObs.METADATA_FILE.value]
        url_of_readme           = planet_metadata.get(mteObs.README_FILE.value, None) # may be None
        if url_of_readme == "":
            url_of_readme = None
        url_of_additional_info  = planet_metadata.get(mteObs.CONFIG_FILE.value, None) # may be None
        if url_of_additional_info == "":
            url_of_additional_info = None

        if url_of_sample not in unique_samples_urls:
            
            unique_samples_urls.append(url_of_sample)
            metadata_urls.append(url_of_metadata)
            readme_urls.append(url_of_readme)
            additional_infos_urls.append(url_of_additional_info)
            planet_metadatas.append(planet_metadata)
           
    
    # iterate over the unique urls and download the samples and metadata and add the information dictionnary to the list
    return_samples = []
    for url_index in range(len(unique_samples_urls)):
        sample_url          = unique_samples_urls[url_index]
        metadata_url        = metadata_urls[url_index]
        readme_url          = readme_urls[url_index]
        additional_info_url = additional_infos_urls[url_index]
        
        planet_metadata = planet_metadatas[url_index]
        file_sample = requests.get(sample_url, verify=False)
        if not file_sample.ok: 
            raise Exception(f"URL {sample_url} responded with status code: {file_sample.status_code}")
        dataframe_sample = pd.read_parquet(io.BytesIO(file_sample.content))
        
        file_metadata = requests.get(metadata_url, verify=False)
        if not file_metadata.ok: 
            raise Exception(f"URL {metadata_url} responded with status code: {file_metadata.status_code}")
        metadata = json.loads(file_metadata.text)

        if readme_url is not None:
            file_readme = requests.get(readme_url, verify=False)
            if not file_readme.ok: 
                raise Exception(f"URL {readme_url} responded with status code: {file_readme.status_code}")
            readme = file_readme.text
        else:
            readme = None
        
        if additional_info_url is not None:   
            file_additional_info = requests.get(additional_info_url, verify=False)
            if not file_additional_info.ok: 
                raise Exception(f"URL {additional_info_url} responded with status code: {file_additional_info.status_code}")
            additional_infos = json.loads(file_additional_info.text) 
        else:
            additional_infos = None
        
        name = f"{metadata[mteObs.STAR_NAME.value]}_{metadata[mteObs.ANALYSIS_ID.value]}"
        
        # name definitions from the table_entries_catalogs.py file
        samples_dict = {
            obsDict.SAMPLE_NAME.value:                  name,
            obsDict.PLANETS_LIST.value:                 metadata[mteAuthors.PLANETS_LIST.value],
           
            obsDict.SAMPLES.value:                      dataframe_sample,
            obsDict.README.value:                       readme,
            obsDict.ADDITIONAL_INFO.value:              additional_infos,
            mteObs.AUTHOR_NAME.value:                  planet_metadata.get(mteObs.AUTHOR_NAME.value, None),
            mteObs.STAR_NAME.value:                     planet_metadata[mteObs.STAR_NAME.value],                     
            mteObs.ANALYSIS_ID.value:                   planet_metadata[mteObs.ANALYSIS_ID.value],  
            obsDict.CONTACT_EMAIL.value:                 planet_metadata[mteObs.CONTACT_EMAIL.value],  
            mteObs.DEFAULT.value:                       planet_metadata[mteObs.DEFAULT.value],  
            mteObs.ROBUSTNESS.value:                    planet_metadata[mteObs.ROBUSTNESS.value],  
            mteObs.INITIAL_CONDITION_DATE_BJD.value:    planet_metadata[mteObs.INITIAL_CONDITION_DATE_BJD.value],  
            mteObs.NB_PLANETS.value:                    planet_metadata[mteObs.NB_PLANETS.value],  
            mteObs.GAIA_ID.value:                       planet_metadata[mteObs.GAIA_ID.value],  
            mteObs.MASS_PRIORS.value:                   planet_metadata[mteObs.MASS_PRIORS.value],  
            mteObs.ECCENTRICITY_PRIORS.value:           planet_metadata[mteObs.ECCENTRICITY_PRIORS.value],  
            mteObs.TRANSIT_DEFINITION.value:            planet_metadata.get(mteObs.TRANSIT_DEFINITION.value, None),
            mteObs.METHODS.value:                       planet_metadata[mteObs.METHODS.value],  
            mteObs.INSTRUMENTS.value:                   planet_metadata[mteObs.INSTRUMENTS.value],  
            mteObs.BIBTEX.value:                        planet_metadata[mteObs.BIBTEX.value],  
            mteObs.CODE_USED.value:                     planet_metadata[mteObs.CODE_USED.value],  
            mteObs.OTHER_REMARKS.value:                 planet_metadata[mteObs.OTHER_REMARKS.value],  
        }
        
        return_samples.append(samples_dict)
    
        # save the files if the download destination is given
        if download_destination is not None:
            sample_file_name = download_destination / f"{name}.parquet"
            metadata_file_name = download_destination / f"{name}_metadata.json"
            readme_file_name = download_destination / f"{name}_readme.txt"
            dataframe.to_parquet(sample_file_name)
            with open(metadata_file_name, 'w') as f:
                json.dump(metadata, f, indent=4)
            if readme is not None:
                with open(readme_file_name, 'w') as f:
                    f.write(readme)
            if additional_infos is not None:
                additional_infos_file_name = download_destination / f"{name}_additional_infos.json"
                with open(additional_infos_file_name, 'w') as f:
                    json.dump(additional_infos, f, indent=4)
        
    return return_samples




# #######################################################################################################
# SIMULATIONS FUNCTION ##################################################################################
# #######################################################################################################

def get_metadata_simulations():
    """This method retrieves the summary metadata table for the simulations and return it in a pandas dataframe
    
    Returns
    -------
    pd.DataFrame
        The simulations metadata table in a pandas dataframe
    
    """
    url = urls.SIMULATIONS_METADATA_TABLE.value
  
    
    # we query the table from DACE and make sure the status code is sucessful (OK = sucessful http request)
    file = requests.get(url, verify=False)
    if not file.ok: # check if the request was successful
        raise Exception(f"URL {url} responded with status code: {file.status_code}")
    dataframe = pd.read_parquet(io.BytesIO(file.content)) # translate the content of the parquet file into a pandas dataframe
        
    return dataframe  



def download_simulations_run_table(dataframe):
    """
    This method retrieves the runs tables from the selected line of the summary metadata datafame, and returns another dataframe with more details about the chosen run. For now, only allows to retrieve a single run table at a time.
    Parameters
    ----------
    dataframe : pd.DataFrame
        The dataframe containing the metadata of the simulation run that we want to download
    
    Returns
    -------
    pd.DataFrame
        A dataframe with the detailed metadata of the run
        
    """
    
    if isinstance(dataframe, pd.Series):
        raise ValueError("The dataframe should be a pandas DataFrame, not a pandas Series")
        
    if not isinstance(dataframe, pd.DataFrame):
        raise ValueError("The dataframe should be a pandas DataFrame")
    
    
    if dataframe.shape[0] > 1:
        raise ValueError(f"The dataframe should only contain one line (header non included). Given shape: {dataframe.shape}")
    
    url_of_table = dataframe[mteSimSum.URL_OF_TABLE.value].values[0]
    
    file = requests.get(url_of_table, verify=False)
    if not file.ok: 
        raise Exception(f"URL {url_of_table} responded with status code: {file.status_code}")
    dataframe = pd.read_parquet(io.BytesIO(file.content))

    return dataframe




def download_simulations(dataframe, download_destination=None):
    """This method retrieves the simulations present in the given dataframe, and returns a list of dictionnary that contains information about the evolutions and the evolutions themselves. If a download destination is given, it saves the evolutions in the given directory.

    Parameters
    ----------
    dataframe : pd.DataFrame
        The dataframe containing the metadata of the simulations that we want to download.
    download_destination : str or None
        If given, save the simulations in the given directory with their corresponding metadatas.
    
    Returns
    -------
    List[Dict]
        A list of dictionnaries with the simulations informations and the simulations themselves
        keys in of dictionnary are:
        
        - simulation_name:  the name of the simulation (str)
        - additional_infos: the additional infos of the simulation (dict)
        - simulation:       the simulation itself (pandas.DataFrame) 
        - readme:           the readme of the simulation (str|None)
        - contact_email:    the email of the contact person for the simulation (str)
        - author_name:      the name of the author of the simulation (str)
        - planets_list:     the list of planets in the simulation (list of str)
        - star_name:        the name of the star of the simulation (str)
        - run_id:          the run id of the simulation (str)
        - simulation_id:    the simulation id of the simulation (str)
        - run_name:         the run name of the simulation (str)
        - simulation_type:  the type of the simulation (str)
        - physics_implemented: the physics implemented in the simulation (str)
        - nb_planets:       the number of planets in the simulation (int)
        - code_used:        the code used for the simulation (str)
        - bibtex:           the bibtex reference of the simulation (str)
        - other_remarks:    the other remarks of the simulation (str)
    """
    
    # common mistake when manipulating pandas dataframe and does not work with the package so it has its own check
    if isinstance(dataframe, pd.Series):
        raise ValueError("The dataframe should be a pandas DataFrame, not a pandas Series")
        
    if not isinstance(dataframe, pd.DataFrame):
        raise ValueError("The dataframe should be a pandas DataFrame")
    
    
    if download_destination is not None:
        if not isinstance(download_destination, Path):
            download_destination = Path(download_destination)
        download_destination.mkdir(parents=True, exist_ok=True)
    
    
    unique_simulations_urls = []
    metadata_urls = []
    additional_infos_urls = []
    readme_files_urls = []
    planet_metadatas = []
    for row in dataframe.iterrows():
        index_row, metadata_line = row
        
        url_of_sample           = metadata_line.get(mteSim.URL_OF_SIMULATION.value, None)
        url_of_metadata         = metadata_line.get(mteSim.METADATA_FILE.value, None)
        url_of_additional_info  = metadata_line.get(mteSim.CONFIG_FILE.value, None) 
        url_of_readme           = metadata_line.get(mteSim.README_FILE.value, None) # may be None
        if url_of_additional_info == "":
            url_of_additional_info = None
        if url_of_readme == "":
            url_of_readme = None
            pass # no condition on the readme as it may not be present
        if url_of_sample is None or url_of_sample == "":
            raise ValueError(f"Missing URL_OF_SIMULATION in the metadata line: {metadata_line}")
        if url_of_metadata is None or url_of_metadata == "":
            raise ValueError(f"Missing METADATA_FILE in the metadata line: {metadata_line}")
        if url_of_sample not in unique_simulations_urls:
            unique_simulations_urls.append(url_of_sample)
            metadata_urls.append(url_of_metadata)
            additional_infos_urls.append(url_of_additional_info)
            readme_files_urls.append(url_of_readme)
            planet_metadatas.append(metadata_line)

    
    return_samples = []
    for url_index in range(len(unique_simulations_urls)):
        sample_url = unique_simulations_urls[url_index]
        metadata_url = metadata_urls[url_index]
        aditional_info_url = additional_infos_urls[url_index]
        readme_file_rul = readme_files_urls[url_index]
        planet_metadata = planet_metadatas[url_index]
        file_sample = requests.get(sample_url, verify=False)
        if not file_sample.ok:
            raise Exception(f"URL {sample_url} responded with status code: {file_sample.status_code}")
        dataframe_sample = pd.read_parquet(io.BytesIO(file_sample.content))
        
        file_metadata = requests.get(metadata_url, verify=False)
        if not file_metadata.ok:
            raise Exception(f"URL {metadata_url} responded with status code: {file_metadata.status_code}")
        metadata = json.loads(file_metadata.text)
        
        if aditional_info_url is not None:
            file_additional_info = requests.get(aditional_info_url, verify=False)
            if not file_additional_info.ok:
                raise Exception(f"URL {aditional_info_url} responded with status code: {file_additional_info.status_code}")
            additional_infos = json.loads(file_additional_info.text)
        else:
            additional_infos = None
        
        if readme_file_rul is not None:
            file_readme = requests.get(readme_file_rul, verify=False)
            if not file_readme.ok:
                raise Exception(f"URL {readme_file_rul} responded with status code: {file_readme.status_code}")
            readme = file_readme.text
        else:
            readme = None
        
        name = f"{metadata[mteSim.SIMULATION_TYPE.value]}_{metadata[mteSim.RUN_ID.value]}_{metadata[mteSim.SIMULATION_ID.value]}"
        
        samples_dict = {
            simDict.SIMULATION_NAME.value:              name,
            simDict.ADDITIONAL_INFO.value:              additional_infos,
            simDict.SIMULATION.value:                   dataframe_sample,
            simDict.README.value:                       readme,
            
            simDict.CONTACT_EMAIL.value:                planet_metadata[mteSim.CONTACT_EMAIL.value],
            mteSim.AUTHOR_NAME.value:                   planet_metadata.get(mteSim.AUTHOR_NAME.value,None),
            simDict.PLANETS_LIST.value:                 metadata[mteAuthors.PLANETS_LIST.value],
            mteSim.STAR_NAME.value:                     planet_metadata[mteSim.STAR_NAME.value],
            mteSim.RUN_ID.value:                        planet_metadata[mteSim.RUN_ID.value],
            mteSim.SIMULATION_ID.value:                 planet_metadata[mteSim.SIMULATION_ID.value],
            mteSim.RUN_NAME.value:                      planet_metadata[mteSim.RUN_NAME.value],
            mteSim.SIMULATION_TYPE.value:               planet_metadata[mteSim.SIMULATION_TYPE.value],
            mteSim.PHYSICS_IMPLEMENTED.value:           planet_metadata[mteSim.PHYSICS_IMPLEMENTED.value],
            mteSim.NB_PLANETS.value:                    planet_metadata[mteSim.NB_PLANETS.value],  
            mteSim.CODE_USED.value:                     planet_metadata[mteSim.CODE_USED.value],
            mteSim.BIBTEX.value:                        planet_metadata[mteSim.BIBTEX.value],
            mteSim.OTHER_REMARKS.value:                 planet_metadata[mteSim.OTHER_REMARKS.value],
            
            
        }
        samples_dict.update(metadata)    
        
        return_samples.append(samples_dict)
    
        if download_destination is not None:
            simulation_file_name = download_destination / f"{name}.parquet"
            metadata_file_name = download_destination / f"{name}_metadata.json"
            dataframe.to_parquet(simulation_file_name)
            with open(metadata_file_name, 'w') as f:
                json.dump(metadata, f, indent=4) 
            if additional_infos is not None:
                additional_infos_file_name = download_destination / f"{name}_additional_infos.json"
                with open(additional_infos_file_name, 'w') as f:
                    json.dump(additional_infos, f, indent=4)
            if readme is not None:
                readme_file_name = download_destination / f"{name}_readme.txt"
                with open(readme_file_name, 'w') as f:
                    f.write(readme)
        
    return return_samples