BLAST_OPTIONS_COMMON = ['dbsize',
'word_size', 'gapopen', 'gapextend', 'qcov_hsp_perc', 
'xdrop_ungap', 'xdrop_gap', 'xdrop_gap_final', 'searchsp', 
'sum_stats', 'seg', 'soft_masking', 'matrix', 
'threshold', 'culling_limit', 
'window_size', 'num_threads', 
'comp_based_stats']

BLAST_FLAGS = ['lcase_masking', 'ungapped', 'use_sw_tback', 
'save_pssm_after_last_round', 'save_each_pssm']

BLASTP_OPTIONS = ['task']

PSIBLAST_OPTIONS = ['gap_trigger', 'num_iterations', 'out_pssm', 
'out_ascii_pssm', 'pseudocount', 'inclusion_ethresh']
