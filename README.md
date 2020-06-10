System requirements: Windows OS

Decomposition of time series of hyperspectral data collected continuously with high resolution field spectrometers over vegetation targets at seasonal, daily and sub-diurnal time scales using Singular Spectrum Analysis method (Golyandina et al., 2001). The R package spectral.methods (Buttlar, 2015) is used for the SSA decomposition. All SSA calculations are done via the truncated and fast SSA algorithm of Korobeynikov et al, 2020 (R package Rssa).

The code is developed in R environment (R Core Team, 2020), and is specifically compiled to be used with the output of the FloX processing chain (https://www.jb-hyperspectral.com/products/flox/, GUI version>=13). The FloX (JB Hyperspectral Devices UG, Germany) is the commercially available instrument developed with intention to gain insights in short-to-long term vegetation processes by collecting long term time series of F and reflectance. 

The following fields must be included in the input .csv file:
•	doy.dayfract
•	datetime [UTC]
•	SZA
•	Lat
•	Lon
•	E_stability [%]
•	E_stability full [%]
•	sat value L
•	sat value E
•	sat value E2	
•	sat value L full	
•	sat value E full	
•	sat value E2 full
•	SIF_A_sfm [mW m-2nm-1sr-1]
•	SIF_B_sfm [mW m-2nm-1sr-1]
•	PAR [W m-2]
•	NDVI
•	PRI

The code is intended to be used for the quality check of the time series and application of SSA on solar-induced fluorescence (SIF) retrieved in O2A (687 nm) and O2B (760 nm) absorption bands (SIF A and SIF B, Mohammed et al., 2019) and Photochemical Reflectance Index (PRI, Gamon et al., 1992). Decomposition is implemented at long-term (seasonal), diurnal and sub-diurnal time scales using filterTSeriesSSA function of spectral.methods R package. 

The output consists of “ssa_decomposed_time_series.csv” and two plots showing the results of the decomposition “SSA_decomposition_output.png”  and “SSA_decomposition_scatterplots.png”. 

Along with the code, we provide a test dataset of FloX data acquired over winter wheat during 20 days from 30/05/2019 to 18/06/2019. For additional information on the dataset please contact Micol Rossini (micol.rossini@unimib.it) and Mirco Migliavacca (mmiglia@bgc-jena.mpg.de).

This project received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska‐Curie Grant Agreement No 721995.
