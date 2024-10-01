This readme file was generated on 2024-09-30 by Stephanie A. Cunningham
<help text in angle brackets should be deleted before finalizing your document>
<[text in square brackets should be changed for your specific dataset]>


GENERAL INFORMATION

Title of Dataset: 

<provide at least two contacts>
Author Information
Name: Stephanie A. Cunningham
ORCID:
Institution: 
Address: 
Email: steph.cunningham.929@gmail.com

Co-author Information
Name: Jacqueline L. Frair
ORCID:
Institution: SUNY College of Environmental Science and Forestry
Address: 1 Forestry Dr. Syracuse, NY 13210
Email: jfrair@esf.edu

Author/Alternate Contact Information
Name: 
ORCID:
Institution: 
Address: 
Email: 

Date of data collection: Samples were collected during legal trapping seasons for fishers in New York (variable by region), October-December of 2018, 2019, and 2020 

Geographic location of data collection: Across New York state, USA

Information about funding sources that supported the collection of the data:
Rodenticide screening was funded by the Grober Graduate Research Fellowship at SUNY ESF and the Wildlife Health Lab at Cornell University. S. Cunningham was supported by the New York Department of Environmental Conservation via Federal Aid in Wildlife Restoration Act Grant W-173-G.


SHARING/ACCESS INFORMATION

Licenses/restrictions placed on the data: 

Links to publications that cite or use the data: 
Silveira et al. 2024. Drivers of anticoagulant rodenticide exposure in fishers (Pekania pennanti) across the northeastern United States. Frontiers in Ecology and Evolution. DOI: 10.3389/fevo.2024.1304659 

Links to other publicly accessible locations of the data: 
Version control history of code available at: https://github.com/s-cunningham/NYFisherRodenticide

Recommended citation for this dataset: 
Cunningham SA, Frair JL, Jensen PG, and KL Schuler. 2024. Data from: Patterns of rodenticide exposure in a mammalian forest carnivore related to hardwood mast cycles and wildland-urban intermix. 


DATA & FILE OVERVIEW

File List: <list all files (or folders, as appropriate for dataset organization) contained in the dataset, with a brief description>
Data

Code
01_all_years_age_data.R
02_summarizing_AR_results.R
03_estimating_trapping_land_cover.R
04_gis_in_r_woods_woodywetlands.R
05_spatial_layers.R
06_checking_spatial_autocorrelation.R
06_combining_covariates.R
07_bromadiolone_iterations.R
07_brodifacoum_iterations.R
07_diphacinone_iterations.R
07_ncompounds_iterations.R
08_plotting_brodifacoum_model.R
08_plotting_bromadiolone_model.R
08_plotting_diphacinone_model.R
09_plotting_binary_models_together.R
09_plotting_ncompounds_predicted.R
10_plotting_AR_distributions.R

Relationship between files, if important:
Code is numbered in the order of processing and analysis.  

Additional related data collected that was not included in the current data package: 

Are there multiple versions of the dataset? No
If yes, name of file(s) that was updated: 
Why was the file updated? 
When was the file updated? 


METHODOLOGICAL INFORMATION

Description of methods used for collection/generation of data: <include links or references to publications or other documentation containing experimental design or protocols used in data collection>

Methods for processing the data: <describe how the submitted data were generated from the raw or collected data>
We used a multiple imputation approach to produce missing point locations of harvested fishers (Murray, 2018). Similar methods have been used successfully for imputing missing data in epidemiological (Sterne et al., 2009) and wildlife (Blanchong et al., 2006; Frair et al., 2004; McClintock, 2017) applications. We generated 10 random locations within each combination of town and WMU to represent potential locations for each sample. When generating random points, we removed areas of town/WMU polygons classified by NLCD as open water, barren land, emergent herbaceous wetlands, developed, or agriculture (hay/pasture and cultivated crops) because we expected trapping of fishers to be unlikely to potentially impossible in these areas. Matching fisher harvest locations from Jensen and Humphries (2019)—collected between 2005 and 2013—to the 2008 NLCD layer indicated that 84% of fishers harvested in northern NY were in forested areas, with an additional 11% trapped in woody wetlands; therefore, because forest and woody wetlands accounted for 95% of harvested fishers, we used these landcover classifications to assign our random points (Table S1). Likewise, we restricted random points to be generated ≤2 km from a road, as trappers generally arrange their traplines along roads (Hodgman et al., 1994; Wiebe et al., 2013).


Instrument- or software-specific information needed to interpret the data: <include full name and version of software, and any necessary packages or libraries needed to run scripts>





Describe any quality-assurance procedures performed on the data: 
Precise trapping locations were not recorded, so we used the town, county, and wildlife management unit (WMU) information recorded on the furbearer possession tag. Because WMU and town boundaries are often arbitrary rather than following roads or natural features, we inspected location information and corrected parts of the record attributes that did not align, i.e., we edited town or WMU if two of the three location attributes matched but not the third, presuming that the third was incorrectly recorded. We first determined the center point of the smallest spatial unit for each harvested fisher (i.e., usually the intersection of town and WMU).

People involved with sample collection, processing, analysis and/or submission: 
Fisher carcasses were submitted by licensed trappers across New York. New York State Department of Environmental Conservation staff extracted tissue samples.


 


DATA-SPECIFIC INFORMATION FOR: [FILENAME]
<repeat this section for each dataset, folder or file, as appropriate>

Number of variables: 

Number of cases/rows: 

Variable List: <list variable name(s), description(s), unit(s) and value labels as appropriate for each>

Missing data codes: <list code/symbol and definition>

Specialized formats or other abbreviations used: 





DATA-SPECIFIC INFORMATION FOR: [FILENAME]
<repeat this section for each dataset, folder or file, as appropriate>

Number of variables: 

Number of cases/rows: 

Variable List: <list variable name(s), description(s), unit(s) and value labels as appropriate for each>

Missing data codes: <list code/symbol and definition>

Specialized formats or other abbreviations used: 