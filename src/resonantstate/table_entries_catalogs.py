from enum import Enum


class URLS(Enum):
    """
    Enum class for the urls of the tables.
    """
    OBSERVATIONS_METADATA_TABLE         = "https://dace.unige.ch/downloads/resonant_state/observations/tables/metadata_table.parquet"
    SIMULATIONS_METADATA_TABLE          = "https://dace.unige.ch/downloads/resonant_state/simulations/tables/metadata_table.parquet"
    #SIMULATIONS_SUMMARY_METADATA_TABLE  = "https://dace.unige.ch/downloads/resonant_state/simulations/tables/metadata_table_summary.parquet"
    

class ObservationsSamplesDictionnary(Enum):
    SAMPLE_NAME     = "sample_name"  
    PLANETS_LIST    = "planets_list" 
    CODE            = "code"        
    BIBTEX          = "bibtex"
    CONTACT_EMAIL   = "contact_email"
    SAMPLE          = "sample"
   

# Name of the column for the metadata table created from the author metadata file
class MetaDataTableEntries(Enum):
    
    # common for both activities:
    PLANET_NAME                     = "planet_name"                         # TOI-178 c
    STAR_NAME                       = "star_name"                           # TOI-178
    CONTACT_EMAIL                   = "contact_mail"                        # adrien.leleu@unige.ch
    PERIOD_DAYS                     = "period_days"                         # in days (float)
    RADIUS_PLANET_R_EARTH           = "radius_planet_r_earth"               # in Earth radius (float)
    MASS_PLANET_M_EARTH             = "mass_planet_m_earth"                 # in Earth mass (float)
    
    MASS_STAR_M_SUN                 = "mass_star_m_sun"                     # in solar mass (float)
    RADIUS_STAR_R_SUN               = "radius_star_r_sun"                   # in solar radius (float)
    
    PLANET_STAR_MASS_RATIO           = "planet_star_mass_ratio"             # float
    PLANET_STAR_RADIUS_RATIO         = "planet_star_radius_ratio"           # float
    SEMI_MAJOR_AXIS_AU              = "semi_major_axis_au"                  # in AU (float)
    MEAN_LONGITUDE_DEG              = "mean_longitude_deg"                  # in deg (float)
    ECCENTRICITY                    = "eccentricity"                        # (float)
    LONGITUDE_OF_PERIASTRON_DEG     = "longitude_of_periastron_deg"         # in deg (float )
    INCLINATION_DEG                 = "inclination_deg"                     # in deg (float)
    LONGITUDE_OF_ASCENDING_NODE_DEG = "longitude_of_ascending_node_deg"     # in deg (float)
    INITIAL_CONDITION_DATE_BJD      = "initial_condition_date_BJD"          # 283223.2838382 (float)
    NB_PLANETS                      = "nb_planets"                          # int
    
    
    BIBTEX                          = "bibtex"                              # bibtex 2024A&a...6888A.211L
    CODE_USED                       = "code_used"                           # like bibtex
    OTHER_REMARKS                   = "other_remarks"                       # str
    
    METADATA_FILE                   = "metadata_file"                       # metadata file path (str)
    CONFIG_FILE                     = "additional_info"                     # additional info file path (str)
    #NOTEBOOK_FILE                   = "notebook_file"                      # notebook file path (str)
    README_FILE                     = "readme_file"                         # readme file path (str)
    #URL_OF_SAMPLES                  = "url_of_samples"                     # str
    URL_OF_SAMPLES_PARQUET          = "url_of_samples"                      # str


    
    
    
    ############################### only in activity A ######################
    ANALYSIS_ID                     = "analysis_id"                         # 0,1,2,..
    ROBUSTNESS                      = "robustness"                          # between -1, 0 and 1
    DEFAULTS                        = "default"                             # 0 or 1 (int)
    GAIA_ID                         = "GAIA_ID"                             # str 
    MASS_PRIORS                     = "mass_prior"                          # uniform, log-uniform (str)
    ECCENTRICITY_PRIORS             = "eccentricity_prior"                  # same
    
    METHODS                         = "methods"                             # ["rv", "transit"] (list of str)
    INSTRUMENTS                     = "instruments"                         # ["Tess", "NGTS, "CHEOPS", "ESPRESSO"] (list of str)
    
    
    
    ############################### only in activity B ######################
    RUN_ID                          = "run_ID"                              # 0,1,2,..
    SYSTEM_ID                       = "System_ID"
    SIMULATION_TYPE                 = "simulation_type"
    
    
    ############################### only in author files ####################
    PLANET_LIST                     = "planet_list"

