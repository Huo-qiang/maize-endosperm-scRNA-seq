import os
import pandas as pd

from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names

maize_ex_path = 'maize_expression_data.tsv'
maize_tf_path = 'maize_transcription_factors.tsv'

ex_matrix = pd.read_csv(maize_ex_path, sep='\t')

tf_names = load_tf_names(net1_tf_path)

%%time
network = grnboost2(expression_data=ex_matrix,
                    tf_names=tf_names)

network.to_csv('maize_scRNA_network.tsv', sep='\t', header=False, index=False)
#####repeat 50 times##########
