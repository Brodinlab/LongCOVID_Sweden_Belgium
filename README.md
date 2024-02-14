# LongCOVID_Sweden_Belgium
Code and data to reproduce the figures of the maniscript entitled: "Restrained memory CD8+ T cell responses favors viral persistence and elevated IgG responses in patients with severe Long COVID" ![image](https://github.com/Brodinlab/LongCOVID_Sweden_Belgium/assets/1262298/0e93b5b3-8ffe-4ee8-9128-ad2b38dbeaa1)

## Table of contents
* [General info](#general-info)
* [Dependencies](#dependencies)
* [Repo description](#repo-description)

## General info
This project used multiple omics data:
- Plasma protein expression (Olink - NPX values)
- Cell frequencies (CyTOF - FlowSOM cell frequencies)
- single-cell TCR and mRNA sequencing (BD Rhapsody - TCR sequences and gene counts)
	
## Dependencies
Project is created with:
* RStudio version: 4.2.1
* Unix/Linux

## Repo description
- ```FlowSOM/``` contains scripts for preprocessing CyTOF data and FlowSOM clustering
- ```LC_OCanalysis_clean.R``` used for Olink, cell frequencies, and serological data analysis  
