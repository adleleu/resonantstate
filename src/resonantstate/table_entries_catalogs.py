from enum import Enum

MISSING_VALUE = ""
DEFAULT_VALUE = "DEFAULT_VALUE"

class URLS(Enum):
    OBSERVATIONS_METADATA_TABLE         = "https://dace.unige.ch/downloads/resonant_state/observations/tables/metadata_table.parquet"
    SIMULATIONS_METADATA_TABLE          = "https://dace.unige.ch/downloads/resonant_state/simulations/tables/metadata_table.parquet"
    

class ObservationsSamplesDictionnary(Enum):
    SAMPLE_NAME                     = "sample_name"  
    SAMPLE                          = "sample"
    README                          = "readme"
    NB_PLANETS                      = "nb_planets"
   
    PLANETS_LIST                    = "planets_list" 
    CODE                            = "code"        
    STAR_NAME                       = "star_name"                       
    ANALYSIS_ID                     = "analysis_id"
    CONTACT_EMAIL                   = "contact_email"


class SimulationsDictionnary(Enum):
    SIMULATION_NAME = "simulation_name"  
    PLANETS_LIST    = "planets_list" 
    CODE            = "code"        
    BIBTEX          = "bibtex"
    CONTACT_EMAIL   = "contact_email"
    SIMULATION      = "simulation"
    ADDITIONAL_INFO = "additional_info"
    README          = "readme"


# Name of the column for the metadata table created from the author metadata file
class MetaDataTableEntriesObservations(Enum):
    
    ############################### common for both activities ##############
    PLANET_NAME                     = "planet_name"                         # TOI-178 c
    STAR_NAME                       = "star_name"                           # TOI-178
    ANALYSIS_ID                     = "analysis_id"                         # 0,1,2,..
    CONTACT_EMAIL                   = "contact_mail"                        # adrien.leleu@unige.ch
    AUTHOR_NAME                     = "author_name"                         
    DEFAULT                         = "default"                             # 0 or 1 (int)
    ROBUSTNESS                      = "robustness"                          # between -1, 0 and 1
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
    GAIA_ID                         = "gaia_id"                             # str 
    
    MASS_PRIORS                     = "mass_prior"                          # uniform, log-uniform (str)
    ECCENTRICITY_PRIORS             = "eccentricity_prior"                  # same
    
    TRANSIT_DEFINITION              = "transit_definition"                   # "TESS", "NGTS", "CHEOPS" (str)
    METHODS                         = "methods"                             # ["rv", "transit"] (list of str)
    INSTRUMENTS                     = "instruments"                         # ["Tess", "NGTS, "CHEOPS", "ESPRESSO"] (list of str)
    BIBTEX                          = "bibtex"                              # bibtex 2024A&a...6888A.211L
    CODE_USED                       = "code_used"                           # like bibtex
    OTHER_REMARKS                   = "other_remarks"                       # str
    
    METADATA_FILE                   = "metadata_file"                       # metadata file path (str)
    CONFIG_FILE                     = "additional_info"                     # additional info file path (str)
    README_FILE                     = "readme_file"                         # readme file path (str)
    URL_OF_SAMPLES                  = "url_of_samples"                      # str


class MetaDataTableEntriesSimulations(Enum):
    
    ############################### common for both activities ##############
    SIMULATION_TYPE                 = "simulation_type"
    RUN_ID                          = "run_id"                              # 0,1,2,..
    RUN_NAME                        = "run_name"                             # str
    SIMULATION_ID                   = "simulation_id"
    CONTACT_EMAIL                   = "contact_mail"                        # adrien.leleu@unige.ch
    AUTHOR_NAME                     = "author_name"                         

    PLANET_NAME                     = "planet_name"                         # TOI-178 c
    PERIOD_DAYS                     = "period_days"                         # in days (float)
    RADIUS_PLANET_R_EARTH           = "radius_planet_r_earth"               # in Earth radius (float)
    MASS_PLANET_M_EARTH             = "mass_planet_m_earth"                 # in Earth mass (float)
    
    MASS_STAR_M_SUN                 = "mass_star_m_sun"                     # in solar mass (float)
    RADIUS_STAR_R_SUN               = "radius_star_r_sun"                   # in solar radius (float)
    
    PLANET_STAR_MASS_RATIO          = "planet_star_mass_ratio"             # float
    PLANET_STAR_RADIUS_RATIO        = "planet_star_radius_ratio"           # float
   
    SEMI_MAJOR_AXIS_AU              = "semi_major_axis_au"                  # in AU (float)
    MEAN_LONGITUDE_DEG              = "mean_longitude_deg"                  # in deg (float)
    ECCENTRICITY                    = "eccentricity"                        # (float)
    LONGITUDE_OF_PERIASTRON_DEG     = "longitude_of_periastron_deg"         # in deg (float )
    INCLINATION_DEG                 = "inclination_deg"                     # in deg (float)
    LONGITUDE_OF_ASCENDING_NODE_DEG = "longitude_of_ascending_node_deg"     # in deg (float)

    NB_PLANETS                      = "nb_planets"                          # int
    
    STAR_NAME                       = "star_name"                           # TOI-178
    PHYSICS_IMPLEMENTED             = "physics_implemented"                 # str
    
    BIBTEX                          = "bibtex"                              # bibtex 2024A&a...6888A.211L
    CODE_USED                       = "code_used"                           # like bibtex
    OTHER_REMARKS                   = "other_remarks"                       # str
    
    METADATA_FILE                   = "metadata_file"                       # metadata file path (str)
    CONFIG_FILE                     = "additional_info"                     # additional info file path (str)
    README_FILE                     = "readme_file"                         # readme file path (str)

    
    URL_OF_SIMULATION               = "url_of_simulation"                   # str
    
    
    
class MetaDataTableEntriesSimulationsSummary(Enum):

    RUN_ID                          = "run_id"                         # 0,1,2,..
    RUN_NAME                        = "run_name"                             # str
    SIMULATION_ID                   = "simulation_id"

    CONTACT_EMAIL                   = "contact_mail"                        # adrien.leleu@unige.ch
    
    SIMULATION_TYPE                 = "simulation_type"
    
    BIBTEX                          = "bibtex"                              # bibtex 2024A&a...6888A.211L
    CODE_USED                       = "code_used"                           # like bibtex
    
    NUMBER_OF_RUNS                  = "number_of_runs"                       # int
    MAXIMUM_NUMBER_OF_PLANETS       = "maximum_number_of_planets"             # int
    URL_OF_TABLE                    = "url_of_table"                       # str
    


# Name of the column for the metadata table recieved by the authors
class AuthorsMetaDataTableEntries(Enum):
    
    ############################### common for both activities ##############
    CONTACT_EMAIL                   = "contact_mail"
    AUTHOR_NAME                     = "author_name"
    BIBTEX                          = "bibtex"
    CODE_USED                       = "code_used"
    STAR_NAME                       = "star_name"
    PLANET_LIST                     = "planet_list"
    OTHER_REMARKS                   = "other_remarks"
    
    ############################### only in activity A ######################
    GAIA_ID                         = "gaia_id"                             # str
    INSTRUMENTS                     = "instruments"
    METHODS                         = "methods"
    ANALYSIS_ID                     = "analysis_id"
    ECCENTRICITY_PRIORS             = "eccentricity_prior"
    MASS_PRIORS                     = "mass_prior"
    ROBUSTNESS                      = "robustness"
    INITIAL_CONDITION_DATE_BJD      = "initial_condition_date_BJD"
    DEFAULT                         = "default"
    
    
    ############################### only in activity B ######################
    RUN_ID                          = "run_id"
    SIMULATION_ID                   = "system_id"
    SIMULATION_TYPE                 = "simulation_type"
    PHYSICS_IMPLEMENTED             = "physics_implemented"
    
    
    

        

    
    
    