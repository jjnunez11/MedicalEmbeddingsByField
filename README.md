# MedicalEmbeddingsByField 

Class project for CPSC 503 evaluating how the performance of clinical concept embedding vary by the field of medicine they're related to. 

The current implementation uses five different evaluation methods to judge embeddings from four different embedding sets. 

The project report is located [here](../doc/med_embeddings_by_field.pdf)

Run all analyses by running [main.py](../eval/main.py)

The code for the analyses are located within [../eval](../eval/). The analyses from Choi were adapted, all others were written by me. Auxillary functions are located within the two helper files. Some functions are reused from Choi et al's implementation. The ones with docstrings indicating JJN were written or adapted by me. 

Raw results files are located within [/results/][(../results)

The [/data/](../data) folder conaints the various data used for the project, including embeddings, and the various dictionaries. Many come from Choi et al's implementation.

## Following data must be obtained prior to runnig: ##

The actual embeddings are not incloded, nor are the UMN SRS similairiy tables, nor is the MRCONSO.RRF file.
- The first three embeddings can be found at [Choi et al's Github](https://github.com/clinicalml/embeddings)
- Beam et al's embeddings can be found [here](https://figshare.com/s/00d69861786cd0156d81)
- The UMN SRS tables can be found [here](https://conservancy.umn.edu/handle/11299/196265)
- MRCONSO.RRF is not needed for this project so far, but is used for building some of the dictionaries in Choi et al's work. It can be downloaded from the NIH ULMS site after autoriziation is obtained. 
