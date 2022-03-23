## MSc Data Science project
#### Student: Alexandra-Gabriela Gheorghe
#### S/N: 13162760
#### Supervisor: Dr. Cen Wan

### Project title: Exploring the use of protein interaction networks and functional similarity metrics for determining the longevity impact of ageing genes in *Caenorhabditis elegans*

This repository contains the source code and important files generated as part of the MSc Data Science project.

The structure of this repository is as follows:
1. The **source code** folder contains the obo2csv.py file used to convert the Gene Ontology to csv format. It also contains the node2vec.py implementation file for the node2vec algorithm. Credits and references are included in the README file associated with the folder. 
2. The **source datasets** folder contains the original go-basic.OBO file, the original GenAge ageing genes dataset, as well as references to the STRING PPI resources used. These files were extensively manipulated during the course of this project.
3. The **GO DAG generation** folder contains the .csv files generated from the Gene Ontology, the Jupyter notebook hosting the code for generating aspect-based GO DAGS, and the files that store GO DAGs as edgelists.
4. The **GenAge clean-up, gene products annotations and propagation** folder contains the Jupyter notebook with the GenAge data preparation steps, extraction of protein annotations, extraction of STRING IDs and propagation of all annotations. The output is comprised of a genage_full file of the clean, filtered GenAge dataset with associated gene products IDs, and a dictionary file of all BP-annotated STRING IDs as keys and their annotations (including propagations) as values. Please see report for full methodology, in particular the **Data acquisition, clean-up and preparation** subsection.
5. The **GO embedding, functional similarity, network scores** folder contains two Jupyter notebooks, one for the embedding of the GO DAG and computation of MHD for determining functional similarity scores for the PPI EXP subnetwork, the other for further processing of the interaction confidence score, generation of the aggregated score, and generation of PPI EXP edgelists weighed by each score respectively (IC, FSim, Aggregated score).
6. The **Classification with score and dimensionality tuning** folder contains the Jupyter notebook showing the code for embedding the PPI EXP subnetwork using different scores as a guide, at different output dimensionalities. It also contains the classification of ageing genes based on the extracted embeddings as predictors, using a range of classifiers and reporting different performance metrics. The folder also contains a subfolder of all the embeddings generated at this stage. See subsection 3.4 of the report.
7. The **Classification with further parametrization** folder contains the Jupyter notebook showing the experiments performed during the additional parametrization stage (see subsection 3.4.1 of the methodology), with the associated ageing gene classification steps. A subfolder of all the embeddings generated is included.

Please see the report submitted to Birkbeck college as part of this project. Please consult the README file included in each folder.
