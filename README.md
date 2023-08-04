## CISTAR-flowsheets
This folder contains folders, sub-folders, and files that are in active development for multiscale shale gas upgrading process optimization. The following sub-sections will highlight the contents of this folder.

### archived
This folder contains archived calculation files and folders. Validity of calculations using files and functions from this folder are not guaranteed.

### data
This folder contains parameter and feed composition data in .xlsx and .csv files.

### initialization_files
This folder contains solution .json files with process model state data that can be used to initialize the IDAES flowsheet.

### plots
This folder contains .pdf version of plots generated using matplotlib

### results
This folder contains .csv and .xlsx files with optimal solution results
* optimal_data_wrt_c_tax_rates.csv: Optimal solution and decision variables for optimization with different GHG emissions tax rates.
* optimal_data_wrt_region.csv: Optimal solution and decision variables for optimization with different NGL feed compositions.
* optimal_data_wrt_ROK_models.csv: Optimal solution and decision variables for optimization with different ROK models.
* solution_data.xlsx: Master file with stream data and heat integration data for all optimal solutions.

### src
This folder contains python files used to define functions for calculations.

### Jupyter notebooks
All the analysis is performed using Jupyter Notebooks. Open and refer to a notebook to learn about the calculations performed in it.
