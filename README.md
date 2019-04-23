# MedicalEmbeddingsByField 

Class project for CPSC 503 evaluating how the performance of clinical concept embedding vary by the field of medicine they're related to. 

The current implementation uses five different evaluation methods to judge embeddings from four different embedding sets. 

The project report is located [here](https://github.com/jjnunez11/MedicalEmbeddingsByField/blob/master/doc/med_embeddings_by_field.pdf)

Run all analyses by running [main.py](https://github.com/jjnunez11/MedicalEmbeddingsByField/blob/master/eval/main.py)

The code for the analyses are located within [`eval`](https://github.com/jjnunez11/MedicalEmbeddingsByField/tree/master/eval). The first three were written from scratch, the last two adapted from [Choi et al's Github](https://github.com/clinicalml/embeddings)
- `analysis_beam_boostrap.py` runs an analysis based on the method in Beam et al's paper 
- `analysis_new_sysvec.py`  runs a new analysis using a medical system's repersentative vector 
- `analysis_yu_umnsrs_cor.py` runs an analysis from Yu et al comparing vector similiarity vs human judgements. 
- `analysis_choi_mrp.py` calculates Medical Relatedness Property from Choi et al
- `analysis_choi_mcsp.py` calculates Medical Conceptual Similiarity Property from Choi et al

Additionally, I needed to write new functions shared by the difference analyses, in addition to adapting those from Choi et al's work. The functions I added or changed are alos located within `eval`:
- `embed_helpers.py` which contain functions for reading and processing the embedding files
- `cui_icd9_helpers.py` which contains function for flitering, converting, and other operations on the dictionaries between CUI and ICD9

Raw results files are located within `results` folder

The `data` folder conaints the various data used for the project, including embeddings, and the various dictionaries. Many come from Choi et al's implementation.

## Following data must be obtained prior to runnig: ##

The actual embeddings are not incloded, nor are the UMN SRS similairiy tables, nor is the MRCONSO.RRF file.
- The first three embeddings can be found at [Choi et al's Github](https://github.com/clinicalml/embeddings)
- Beam et al's embeddings can be found [here](https://figshare.com/s/00d69861786cd0156d81)
- The UMN SRS tables can be found [here](https://conservancy.umn.edu/handle/11299/196265)
- MRCONSO.RRF is not needed for this project so far, but is used for building some of the dictionaries in Choi et al's work. It can be downloaded from the NIH ULMS site after autoriziation is obtained. 
