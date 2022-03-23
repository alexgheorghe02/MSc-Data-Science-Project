This folder contains the implementation of:

1. obo2csv - The obo2csv.py file contains the workflow used to convert the Gene Ontology from obo to csv format. The reference code was imported from the public repository of Chengxin Zhang of University of Michigan. Please see https://github.com/kad-ecoli/python_scripts/blob/master/obo2csv.py. The original code was adapted to include all the strictly transitive relationships in the Gene Ontology (only 'is_a' was considered in the original), to remove blocks of code that were redundant and not needed for the scope of this project, and adjust the rest of the implementation for compatibility with these changes.

2. node2vec algorithm - The node2vec.py file contains the implementation of node2vec, based on the reference code by Grover and Leskovec. The hub attention implementation was modified from Schliski et al., for compatibility with the alias sampling proposed by Grover and Leskovec. Additionally, the code was adapted to read graphs with string format node labels, and for Skip-gram determinism. Please see: 

* A. Grover and J. Leskovec. node2vec: Scalable Feature Learning for Networks. KDD '16: Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, pages 855–864, 2016. Github reference: https://github.com/aditya-grover/node2vec
* F. Schliski, J. Schlötterer, and M. Granitzer. Influence of Random Walk Parametrization on Graph Embeddings. Advances in Information Retrieval. ECIR 2020. Lecture Notes in Computer Science, 12036, 2020. Github reference: https://github.com/Kombustor/submission-ecir2020-randomwalks/


