# ResonantState database [draft]

This package contains basic functions that helps navigate the resonant state databases, please check the notebooks https://github.com/adleleu/resonantstate_nb to see how the package can be used.

The resonant state database is split into two : 



- **Resonant State Observations Catalog** : https://dace.unige.ch/resonantStateObservationsCatalog/ this database stores samples of MCMC posterior of the Transit Timing Variations (TTVs) or photo-dynamical analysis of observed (near-)resonant planetary systems.

- **Resonant State Simulations Catalog** : https://dace.unige.ch/resonantStateSimulationsCatalog/ this database stores simlulated planetary systems and their temporal evolution. The simulated physics can vary : 1D proto-planetary discs, hydro-dynamical simulations, tidal evolution, etc.


## Resonant State Observations Catalog 

### The observation metadata
The entries to this catalog are sumarized in a metadata table that can be queried, and visualised on the DACE website : [*link*](https://dace.unige.ch/resonantStateObservationsCatalog/). Each line describe a planet that was characterised in an analysis. The couple star_name + analysis_id is an unique identifiant. 



- **planet_name** (*char*):  Planet name on the Nasa Exoplanet Archive 
- **star_name** (*char*):  Star name on the Nasa Exoplanet Archive 
- **analysis_id** (*int*): Unique analysis indentifiant per star 
- **contact_mail** (*char*): Email adress of the contact person
- **default** (*int*): 1 if the author consider that this is the analysis to use.
- **robustness** (*int*): was the planetary masses tested for robustness to degeneracies (see for example : https://ui.adsabs.harvard.edu/abs/2017AJ....154....5H/abstract, https://ui.adsabs.harvard.edu/abs/2023A%26A...669A.117L/abstract): -1 if the robustness was not tested, 0 if the mass was found not robust, 1 if the mass was found robust
- **period_days** (*float*): orbital period of the planet in days
- **radius_planet_r_earth** (*float*): median of the radius of the planet in Earth radii
- **mass_planet_m_earth** (*float*): median of the mass of the planet in Earth masses
- **mass_star_m_sun** (*float*): median of the mass of the planet in Sun masses
- **radius_star_r_sun** (*float*): median of the radius of the planet in Sun radii
- **planet_star_mass_ratio** (*float*): median of the mass of the planet in Stellar masses
- **planet_star_radius_ratio** (*float*): median of the radius of the planet in Stellar radii
- **semi_major_axis_au** (*float*): median of the semi-major axis  of the planet in AU
- **mean_longitude_deg** (*float*): median of the mean longitude of the planet in dregree
- **eccentricity** (*float*): eccentricity of the planet
- **longitude_of_periastron_deg** (*float*): median of the longitude of periastron of the planet in dregree
- **inclination_deg** (*float*): median of inclination of the planet in dregree
- **longitude_of_ascending_node_deg** (*float*): median of the longitude of ascending node of the planet in dregree
- **initial_condition_date_BJD** (*float*): date at which the orbital elements are determined in Barycentric Julian Date.
- **nb_planets** (*int*): number of planets analysed
- **GAIA_ID** (*str*): indentifiant of a Gaia data release
- **mass_prior** (*str*): type of prior used for the planetary mass (e.g. uniform, log-uniform, $\beta$ distribution, etc.)
- **eccentricity_prior** (*str*): type of prior used for the planetary eccentricity (e.g. uniform, log-uniform, $\beta$ distribution, etc.)
- **methods** (list of *str*): list of method used (e.g. TTVs, photo-dynamic, RVs)
- **instruments** (list of *str*): list of instruments whose data was used in the analysis (e.g. Kepler, TESS, CHEOPS, ESPRESSO)
- **bibtex** (*str*): reference for the analysis
- **code_used** (*str*): code used for the analysis
- **other_remarks** (*str*): other remarks provided by the contact person
- **metadata_file** (*str*): link to the metadata.json that was provided by the contact person for this analysis
- **additional_info** (*str*): link to the additional_info.json file that was provided by the contact person for this analysis
- **readme_file** (*str*): link to the readme_file.txt file that was provided by the contact person for this analysis
- **url_of_samples** (*str*): link to the samples for this analysis


### The samples 

The samples are stored in a parquet file. For an analysis with N planets, the columns are as follow:

- **Column  0** :   sample Identifier 
- **Columns 8*j+1** : Mean longitude $\lambda_j$ of planet j, in degrees 
- **Columns 8*j+2** : Period of planet j, in days 
- **Columns 8*j+3** : $k_j = e_j*\cos(\varpi_j)$, where $e_j$ is the eccentricity and $\varpi_j$ the longitude of the periapsis of planet j 
- **Columns 8*j+4** : $h_j = e_j*\sin(\varpi_j)$ , where $e_j$ is the eccentricity and $\varpi_j$ the longitude of the periapsis of planet j 
- **Columns 8*j+5** : Inclination $I_j$ of planet j in degrees (set to 90 for coplanar analysis)
- **Columns 8*j+6** : Longitude of the ascending node $\Omega_j$ of planet j, in degrees (set to 0 for the outermost planet, set to 0 for all planets for coplanar analysis) 
- **Columns 8*j+7** : $m_j/M_*$, where $m_j$ is the mass of planet j and $M_*$ is the stellar mass 
- **Columns 8*j+8** : $R_j/R_*$ Radii of the planets divided by stellar radius (set to literature value if not fitted)
- **Column  8*N+1** : Mass of the star in Sun masses (set to literature value if not fitted) 
- **Column  8*N+2** : Radii of the star in Sun radii  (set to literature value if not fitted) 

The number of planets N is the same as the length of the list in the field **planet_list** of the metadata.json file. The ordering is the same as well, so that list can be used to retrieve the name of each planet.

Additional columns can be added by the contact person. If so, their name and unit in should be defined in the Readme.txt / additional_info.json associated to that analysis.



## The Resonant State Simulations Catalog 

### The simulation metadata


The entries to this catalog are sumarized in a metadata table that can be queried, and visualised on the DACE website : *[link](https://dace.unige.ch/resonantStateSimulationsCatalog/)*. each line correspond to a run of simulations. A run is typically a setup that was run for one or more time with different initial conditions. Each line contain the link to a parquet file that describe the simulated planetary systems similarly to the Observations Catalog metadata table (see above).



- **contact_mail** (*char*): Email adress of the contact person.
- **simulation_type** (*str*): type of simulation (population synthesis, tidal evolution, etc.).
- **run_id** (*int*) unique identifier for a given simulation type.
- **run_name** (*str*) name given at the run by the author.


- **bibtex** (*str*): reference for the analysis
- **code_used** (*str*): code used for the analysis

- **maximum_number_of_planets** (*int*) average number of planet per simulated system


[to be completed]


### The run metadata


- **contact_mail** (*char*): Email adress of the contact person.
- **simulation_type** (*str*): type of simulation (population synthesis, tidal evolution, etc.).
- **run_id** (*int*) unique identifier for a given simulation type.
- **run_name** (*str*) name given to the run by the author.

- **simulation_id** (*int*) unique identifier for a given simulation within a run.
- **planet_name** (*char*):  planet name.

- **period_days** (*float*): orbital period of the planet in days.
- **radius_planet_r_earth** (*float*): median of the radius of the planet in Earth radii.
- **mass_planet_m_earth** (*float*): median of the mass of the planet in Earth masses.
- **mass_star_m_sun** (*float*): median of the mass of the planet in Sun masses.
- **radius_star_r_sun** (*float*): median of the radius of the planet in Sun radii
- **planet_star_mass_ratio** (*float*): median of the mass of the planet in Stellar masses
- **planet_star_radius_ratio** (*float*): median of the radius of the planet in Stellar radii
- **semi_major_axis_au** (*float*): median of the semi-major axis of the planet in AU
- **mean_longitude_deg** (*float*): median of the mean longitude of the planet in dregree
- **eccentricity** (*float*): median of eccentricity of the planet
- **longitude_of_periastron_deg** (*float*): median of the longitude of periastron of the planet in dregree
- **inclination_deg** (*float*): median of inclination of the planet in dregree
- **longitude_of_ascending_node_deg** (*float*): median of the longitude of ascending node of the planet in dregree

- **nb_planets** (*int*): number of surviving planets at the end of the simulation.

- **star_name** (*char*):  Star name 

- **physics_implemented** (list of *str*): list of physical processes that are implemented (e.g. protoplanetary disk, tides, etc.)

- **bibtex** (*str*): reference for the simulation
- **code_used** (*str*): code used for the simulation
- **other_remarks** (*str*): other remarks provided by the contact person
- **metadata_file** (*str*): link to the metadata.json that was provided by the contact person for this simulation
- **additional_info** (*str*): link to the additional_info.json file that was provided by the contact person for this simulation
- **readme_file** (*str*): link to the readme_file.txt file that was provided by the contact person for this simulation
- **url_of_simulation** (*str*): link to the temporal evolution of the simulation


[to be completed]
    

### The simulation files 


The simulation are stored in a parquet file. For an analysis with N planets, the 8*N+2 first columns are as follow:

- **Column  0** :   time stamp 
- **Columns 8*j+1** : Mean longitude lambda_j of planet j, in degrees 
- **Columns 8*j+2** : Period of planet j, in days 
- **Columns 8*j+3** : $k_j = e_j*\cos(\varpi_j)$, where $e_j$ is the eccentricity and $\varpi_j$ the longitude of the periapsis of planet j 
- **Columns 8*j+4** : $h_j = e_j*\sin(\varpi_j)$ , where $e_j$ is the eccentricity and $\varpi_j$ the longitude of the periapsis of planet j 
- **Columns 8*j+5** : Inclination $I_j$ of planet j in degrees (set to 90 for coplanar analysis)
- **Columns 8*j+6** : Longitude of the ascending node $\Omega_j$ of planet j, in degrees (set to 0 for the outermost planet, set to 0 for all planets for coplanar analysis) 
- **Columns 8*j+7** : $m_j/M_*$, where m_j is the mass of planet j and M* is the stellar mass 
- **Columns 8*j+8** : $R_j/R_*$ Radii of the planets divided by stellar radius (set to literature value if not fitted)
- **Column  8*N+1** : Mass of the star in Sun masses (set to literature value if not fitted) 
- **Column  8*N+2** : Radii of the star in Sun radii  (set to literature value if not fitted) 

The number of planets N is the same as the length of the list in the field **planet_list** of the metadata.json file. The ordering is the same as well, so that list can be used to retrieve the name of each planet.

Additional columns can be added by the contact person. If so, their name and unit in should be defined in the Readme.txt / additional_info.json metadata associated to that analysis.


[to be completed]