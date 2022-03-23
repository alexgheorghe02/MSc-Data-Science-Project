The GO DAG embedding, MHD.ipynb notebook contains the code for embedding the GO DAG (see GO embeddings.zip file), generating the PPI EXP subnetwork (see PPI_BP_edgelist_EXP.txt file), functions for calculating the MHD and computing fuctional similarity scores for all protein pairs in the PPI EXP subnetwork (See FSim_PPI_EXP.txt file).

The Aggregating scores.ipynb notebook contains the code for scaling the interaction confidence score to values [0, 1], aggregating the interaction confidence and functional similarity into a new score, and generating files representing the PPI EXP subnetwork weighed by different scores:

1. For IC score - PPI_EXP_scaled_IC.txt
2. For FSim score - FSim_PPI_EXP.txt
3. For Aggregated score - PPI_EXP_aggregated_scores.txt
