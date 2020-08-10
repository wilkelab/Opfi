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


def build_blastp_command(query, db, evalue, kwargs, output_fields, out):
    """
    Generate a python subprocess command (list) for
    executing blastp.

    kwargs should be valid blastp options. Note that
    at this time, only certain options are allowed.
    """

    cmd = ["blastp", "-query", query, "-db", db, "-evalue", str(evalue), "-out", out]
    for key, value in kwargs.items():
        if key in (BLAST_OPTIONS_COMMON + BLASTP_OPTIONS):
            cmd.append("-{}".format(key))
            cmd.append(str(value))
        elif key in (BLAST_FLAGS) and value:
            cmd.append("-{}".format(key))

    cmd.append("-outfmt")
    cmd.append("6 {}".format(output_fields))
    return cmd

def build_psiblast_command(query, db, evalue, kwargs, output_fields, out):
    """
    Generate a python subprocess command (list) for
    executing psiblast.

    kwargs should be valid blastp options. Note that
    at this time, only certain options are allowed.
    """

    cmd = ["psiblast", "-query", query, "-db", db, "-evalue", str(evalue), "-out", out]
    for key, value in kwargs.items():
        if key in (BLAST_OPTIONS_COMMON + PSIBLAST_OPTIONS):
            cmd.append("-{}".format(key))
            cmd.append(str(value))
        elif key in (BLAST_FLAGS) and value:
            cmd.append("-{}".format(key))

    cmd.append("-outfmt")
    cmd.append("6 {}".format(output_fields))
    return cmd
