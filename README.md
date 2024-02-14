# LongCOVID_Sweden_Belgium
Code to reproduce the figures of the manuscript entitled: "Restrained memory CD8+ T cell responses favors viral persistence and elevated IgG responses in patients with severe Long COVID" 

Currently available on [MedRxiv](https://dx.doi.org/10.1186%2Fs13059-014-0550-8)

## Table of contents
* [General info](#general-info)
* [Dependencies](#dependencies)
* [Repo description](#repo-description)

## General info
This project used multiple omics data:
- Plasma protein expression (Olink - NPX values)
- Cell frequencies (CyTOF - FlowSOM cell frequencies)
- Single-cell TCR and mRNA sequencing (BD Rhapsody - TCR sequences and gene counts)
	
## Dependencies
Project is created with:
* RStudio version: 4.2.1
* Unix/Linux

## Repo description
- ```FlowSOM/``` contains scripts for preprocessing CyTOF data and FlowSOM clustering
- ```LC_OCanalysis_clean.R``` used for Olink, cell frequencies, and serological data analysis  
