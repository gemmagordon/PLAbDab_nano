# PLAbDab_nano

The database can be downloaded [here](https://opig.stats.ox.ac.uk/webapps/plabdab-nano/)

- For all entries: all_sequences.csv.gz 
- For separated VHH and VNAR entries: vhh_sequences.csv.gz, vnar_sequences.csv.gz

**still requires shark models for search - for VNAR annotation - may need to change before can make public (depends on ANARCII release)**

model at:
/vols/bitbucket/gordon/plabdabnano/pdn_reqs/shark_model.pt

needs to be in: 
PLAbDab-nano_project/test_repos/PLAbDab_nano/PLAbDab_nano/shark_number


---

<div align="center">    
 
# PLAbDab-nano: a database of camelid and shark nanobodies from patents and literature

---
    
by 
Gemma Gordon $^{1}$, Alexander Greenshields-Watson $^{1}$, Parth Agarwal $^{1}$ Ashley Wong $^{1}$, Alissa Hummer $^{1}$, Fergus Boyles $^{1}$, Ana Lujan $^{2}$, and Charlotte M. Deane $^{1,*}$

$^{1}$ Oxford Protein Informatics Group, Department of Statistics, University of Oxford, Oxford, United Kingdom  
$^{2}$ Twist Bioscience, South San Francisco, California, US  
$^{*}$ To whom correspondence should be addressed  
</div>

![PDN_logo_text](https://github.com/oxpig/PLAbDab_nano/assets/72207136/7631ffb5-251b-4494-979f-9e6527698f8c)

## Abstract
Nanobodies are essential proteins of the adaptive immune systems of camelids and cartilaginous fish, complementing conventional antibodies. Properties such as their smaller size, greater solubility and high thermostability make VHH and VNAR modalities a promising therapeutic format and a valuable resource for a wide range of biological applications. The volume of academic literature and patents related to nanobodies has risen significantly over the past few decades. Here, we present PLAbDab-nano, a nanobody accompaniment to the existing Patent and Literature Antibody Database (PLAbDab). PLAbDab-nano is a self-updating, searchable repository containing approximately 5000 annotated VHH and VNAR sequences. We describe the methods used to curate the entries in PLAbDab-nano, and highlight how PLAbDab-nano could be used to design diverse libraries, as well as compare and improve query sequences against those in the database. PLAbDab-nano is freely available as a searchable web server (https://opig.stats.ox.ac.uk/webapps/plabdab-nano/). 

## Usage

To use the PLAbDab-nano python package you must first download a copy of the database. This can be obtained from  <a href="https://opig.stats.ox.ac.uk/webapps/plabdab-nano/">here.</a> The file will come in a compressed format. The database columns are as follows: 

| Column name       | Contents                                                                                                                                                               |
|-------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| source            | Whether the entry was sourced from GenBank, SAbDab or TheraSAbDab                                                                                                      |
| model             | ID for the predicted structure of the sequence.  Where sequences are redundant they may share the same model.  VNARs cannot be modelled by NBB2 and so read 'FAILED'.  |
| type              | VHH or VNAR                                                                                                                                                            |
| cdr_lengths       | Lengths of CDR1,2,3 loops. VNARs do not have a CDR2 loop.                                                                                                              |
| cdr_sequences     | Dict of CDR sequences. VNARs do not have a CDR2 loop.                                                                                                                  |
| sequence          | The VHH or VNAR sequence.                                                                                                                                              |
| ID                | ID for the entry, taken from GenBank ID, the PDB ID or therapeutic name.                                                                                               |
| definition        | Title of individual sequence entry from GenBank, PDB or patent.                                                                                                        |
| reference_authors | Authors that published the sequence and/or structure.                                                                                                                  |
| reference_title   | Publication title associated with sequence.                                                                                                                            |
| organism          | Species from which nanobody was derived.                                                                                                                               |
| update_date       | Date that sequence was published.                                                                                                                                      |
| targets_mentioned | Possible antigen targets of nanobody.                                                                                                                                  |

Once the database file is extracted, you can point your version of PLAbDab-nano to it as follows:

```python
from PLAbDab_nano import PLAbDab_nano

plabdabnano = PLAbDab_nano(data_directory, n_jobs = 10)
```

This ''plabdabnano'' object can then be used to perform searches. PLAbDab-nano will output a pandas dataframe for each search. For more info on how to use the API, please check the notebooks. 

## Install

PLAbDab-nano requires NanoBodyBuilder2 from <a href="https://github.com/oxpig/ImmuneBuilder">ImmuneBuilder</a> to model queries when using the structure search. To install it please follow the instructions on the <a href="https://github.com/oxpig/ImmuneBuilder">ImmuneBuilder repo.</a>

<a href="https://github.com/oxpig/kasearch">KA-search</a>, <a href="https://github.com/oxpig/ANARCI">ANARCI</a>, and BLAST (blastp) are also required. Once these are installed, PLAbDab-nano can be installed by running the following command in the PLAbDab-nano directory:

```bash
$ pip install .
```

## Citing this work

The code and data in this package is based on the <a href="https://doi.org/10.1093/nar/gkae881">following paper.</a> If you use it, please cite:

```tex
@article{Gordon2024,
  title = {PLAbDab-nano: a database of camelid and shark nanobodies from patents and literature},
  volume = {53},
  ISSN = {1362-4962},
  url = {http://dx.doi.org/10.1093/nar/gkae881},
  DOI = {10.1093/nar/gkae881},
  number = {D1},
  journal = {Nucleic Acids Research},
  publisher = {Oxford University Press (OUP)},
  author = {Gordon,  Gemma L and Greenshields-Watson,  Alexander and Agarwal,  Parth and Wong,  Ashley and Boyles,  Fergus and Hummer,  Alissa and Lujan Hernandez,  Ana G and Deane,  Charlotte M},
  year = {2024},
  month = oct,
  pages = {D535–D542}
}
```
