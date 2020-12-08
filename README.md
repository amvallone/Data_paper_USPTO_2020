# Data_paper_USPTO_2020

This repository contains the data and coding used in two under review papers

The data include a set of 13 tables, each one corresponding to one calendar year from 2007 to 2019. In the data repository, each table is name as a combination of the letter “t” and the granting year (i.e., the file “t2007.rdata” contains the 2007 data table).  Each table contains the fallowing variables:
- patent_id: Unique patent identification number in the USPTO system
- inventor_id: Unique inventor identification number in the USPTO system 
- pc_section: The patent IPC section 
- ipc_class: The patent IPC class
- ipc_subclass: The patent IPC subclass
- patent_year: The patent’s granted year 

The file Collect.R contains the codes and the generated function used for data collection.

The file Analyzed.R contains the codes which were used in the paper 
