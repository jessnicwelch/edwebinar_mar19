# Introduction to Geospatial Analysis in R  

### NASA Earthdata Webinar  

### *Presented by the ORNL DAAC*  <https://daac.ornl.gov>  

### *March 13, 2019*  

### **Keywords: R, Carbon Monitoring System, GIS**

***
**UPDATE:** The entire repo was updated 2021-12-03 and again on 2024-01-25. The recent update includes a revision of the primary tutorial to replace deprecated packages (e.g., rgdal).

# 1. Overview

The open-source software environment R is gaining popularity among many scientists, including geologists, biologists, and environmental scientists. R provides an alternative to traditional GIS software with numerous packages for geospatial analysis. This webinar will begin with a brief introduction to an example geospatial dataset from the ORNL DAAC and an overview of common geospatial operations in R. Next, we will demonstrate how to import files into R, overlay layers, reduce spatial extent, select and reclassify values, and make a map.

# 2. Dataset

## 2.1 CMS: Forest Carbon Stocks, Emissions, and Net Flux for the Conterminous US: 2005-2010

This dataset provides maps of estimated carbon in forests of the 48 conterminous states of the United States (CONUS) for the years 2005-2010. Committed carbon stocks were estimated for forest aboveground biomass, belowground biomass, standing dead stems, and litter for the year 2005. Carbon emissions were estimated from land use conversion to agriculture, fire damage, insect damage, logging, wind, and weather events in the forests for the years 2006-2010. Committed net carbon flux was estimated as the sum of carbon emissions and sequestration. The maps are provided at 100 m spatial resolution in GeoTIFF format. Average annual carbon estimates, by United States county, are provided in shapefile format for (1) emissions for the multiple disturbance sources, (2) sequestration, and (3) the committed net carbon flux.

Hagen, S., N. Harris, S.S. Saatchi, T. Pearson, C.W. Woodall, S. Ganguly, G.M. Domke, B.H. Braswell, B.F. Walters, J.C. Jenkins, S. Brown, W.A. Salas, A. Fore, Y. Yu, R.R. Nemani, C. Ipsan, and K.R. Brown. 2016. **CMS: Forest Carbon Stocks, Emissions, and Net Flux for the Conterminous US: 2005-2010.** ORNL DAAC, Oak Ridge, Tennessee, USA. <https://doi.org/10.3334/ORNLDAAC/1313>

## 2.2 Carbon Monitoring System (CMS)

The Carbon Monitoring System (CMS) is designed to make significant contributions in characterizing, quantifying, understanding, and predicting the evolution of global carbon sources and sinks through improved monitoring of carbon stocks and fluxes. CMS will use the full range of NASA satellite observations and modeling/analysis capabilities to establish the accuracy, quantitative uncertainties, and utility of products for supporting national and international policy, regulatory, and management activities. CMS will maintain a global emphasis while providing finer scale regional information, utilizing space-based and surface-based data and will rapidly initiate generation and distribution of products both for user evaluation and to inform near-term policy development and planning.

# 3. Prerequisites

Participants should have a basic understanding of R and some exposure to geospatial data and analysis, such as in ArcGIS.

## 3.1 R

1. [Download R](https://cran.r-project.org/)  
2. [Download RStudio](https://www.rstudio.com/products/rstudio/download/#download)  *Recommended*  
3. [Review R Manuals](https://cran.r-project.org/manuals.html)  *Recommended*  

## 3.2 Data

1. [Sign in to NASA Earthdata](https://urs.earthdata.nasa.gov/users/new)  
2. [Download GrossEmissions_v101_USA_Fire.tif](https://daac.ornl.gov/daacdata/cms/CMS_Forest_Carbon_Fluxes/data//GrossEmissions_v101_USA_Fire.tif)  
3. [Download GrossEmissions_v101_USA_Insect.tif](https://daac.ornl.gov/daacdata/cms/CMS_Forest_Carbon_Fluxes/data//GrossEmissions_v101_USA_Insect.tif)  

# 4. Procedure

## 4.1 Tutorial  

1. [R Markdown](https://github.com/jessnicwelch/edwebinar_mar19/blob/master/edwebinar_mar19_ornldaac_tutorial.Rmd)  
2. [Markdown](https://github.com/jessnicwelch/edwebinar_mar19/blob/master/edwebinar_mar19_ornldaac_tutorial.md)  

## 4.2 Supplemental  

1. [R Markdown](https://github.com/jessnicwelch/edwebinar_mar19/blob/master/edwebinar_mar19_ornldaac_supplemental.Rmd)  
2. [Markdown](https://github.com/jessnicwelch/edwebinar_mar19/blob/master/edwebinar_mar19_ornldaac_supplemental.md)  

## 4.3 Webinar  

1. [Presentation Slides](https://daac.ornl.gov/resources/tutorials/r-geospatial-webinar/Earthdata_R_Geospatial_webinar_Mar2019_v20190312.pdf)  
2. [Video](https://daac.ornl.gov/resources/tutorials/r-geospatial-webinar/)

# 5. Credits

* [R](https://www.r-project.org/) - 4.3.2 (2023-10-31) -- "Eye Holes"  
* [RStudio](https://www.rstudio.com/products/rstudio/) - IDE and notebook construction  
* [raster package](https://CRAN.R-project.org/package=raster) - geospatial data manipulations  
* [sf package](https://cran.r-project.org/web/packages/sf/index.html) - support for simple features
* [terra package](https://cran.r-project.org/package=terra) - spatial data analysis with raster (grid) data  
* [tigris](https://CRAN.R-project.org/package=tigris) - states function  
