import os
import pandas as pd

from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names


if __name__ == '__main__':

    in_file  = 'mazie_endosprem_scRNA-seq_count.tsv'
    tf_file  = 'TFgenesallmaize.tsv'
    out_file = 'maize_endosperm_GRN.tsv'

    # ex_matrix is a DataFrame with gene names as column names
    ex_matrix = pd.read_csv(in_file, sep='\t')

    # tf_names is read using a utility function included in Arboreto
    tf_names = load_tf_names(tf_file)

    # compute the GRN
    network = grnboost2(expression_data=ex_matrix,
                        tf_names=tf_names)

    # write the GRN to file
    network.to_csv(out_file, sep='\t', index=False, header=False)
    ###################################repeat 50 times#####################################################
