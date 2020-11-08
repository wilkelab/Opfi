BLAST_OPTIONS_COMMON = ['dbsize',
'word_size', 'gapopen', 'gapextend', 'qcov_hsp_perc', 
'xdrop_ungap', 'xdrop_gap', 'xdrop_gap_final', 'searchsp', 
'sum_stats', 'seg', 'soft_masking', 'matrix', 
'threshold', 'culling_limit', 'window_size', 'num_threads', 
'comp_based_stats', 'gilist', 'seqidlist', 'negative_gilist',
'db_soft_mask', 'db_hard_mask', 'entrez_query', 'max_hsps',
'best_hit_overhang',  'best_hit_score_edge', 'max_target_seqs',
'import_search_strategy', 'export_search_strategy', 'num_alignments']

BLAST_FLAGS = ['lcase_masking', 'ungapped', 'use_sw_tback', 'remote']

PSIBLAST_FLAGS = ['lcase_masking', 'use_sw_tback', 
'save_pssm_after_last_round', 'save_each_pssm', 'remote']

BLASTP_OPTIONS = ['task']

PSIBLAST_OPTIONS = ['gap_trigger', 'num_iterations', 'out_pssm', 
'out_ascii_pssm', 'pseudocount', 'inclusion_ethresh']

BLASTN_OPTIONS = ['filtering_algorithm', 'sum_stats', 'window_masker_db', 'window_size', 'template_type', 'version', 'parse_deflines', 'min_raw_gapped_score', 'string', 'format', 'max_hsps', 'taxids', 'negative_taxids', 'num_alignments', 'strand', 'off_diagonal_range', 'subject_besthit', 'num_sequences', 'no_greedy', 'negative_taxidlist', 'culling_limit', 'xdrop_ungap', 'open_penalty', 'DUST_options', 'sorthits', 'xdrop_gap_final', 'negative_gilist', 'subject', 'use_index', 'bool_value', 'filename', 'seqidlist', 'task_name', 'sort_hits', 'database_name', 'lcase_masking', 'query_loc', 'subject_loc', 'sort_hsps', 'line_length', 'boolean', 'db_hard_mask', 'negative_seqidlist', 'template_length', 'filtering_db', 'filtering_database', 'penalty', 'searchsp', 'ungapped', 'type', 'gapextend', 'db_soft_mask', 'dbsize', 'qcov_hsp_perc', 'sorthsps', 'window_masker_taxid', 'index_name', 'export_search_strategy', 'float_value', 'soft_masking', 'gilist', 'entrez_query', 'show_gis', 'best_hit_score_edge', 'gapopen', 'subject_input_file', 'range', 'html', 'word_size', 'best_hit_overhang', 'perc_identity', 'input_file', 'num_descriptions', 'xdrop_gap', 'dust', 'taxidlist', 'max_target_seqs', 'num_threads', 'task', 'remote', 'int_value', 'extend_penalty', 'reward', 'import_search_strategy', 'num_letters']


def build_blastp_command(query, db, evalue, kwargs, output_fields, out, blastp_path):
    """
    Generate a python subprocess command (list) for
    executing blastp.

    kwargs should be valid blastp options. Note that
    at this time, only certain options are allowed.
    """

    cmd = [blastp_path, "-query", query, "-db", db, "-evalue", str(evalue), "-out", out]
    for key, value in kwargs.items():
        if key in (BLAST_OPTIONS_COMMON + BLASTP_OPTIONS):
            cmd.append("-{}".format(key))
            cmd.append(str(value))
        elif key in (BLAST_FLAGS) and value:
            cmd.append("-{}".format(key))

    cmd.append("-outfmt")
    cmd.append("6 {}".format(output_fields))
    return cmd

def build_psiblast_command(query, db, evalue, kwargs, output_fields, out, psiblast_path):
    """
    Generate a python subprocess command (list) for
    executing psiblast.

    kwargs should be valid blastp options. Note that
    at this time, only certain options are allowed.
    """

    cmd = [psiblast_path, "-query", query, "-db", db, "-evalue", str(evalue), "-out", out]
    for key, value in kwargs.items():
        if key in (BLAST_OPTIONS_COMMON + PSIBLAST_OPTIONS):
            cmd.append("-{}".format(key))
            cmd.append(str(value))
        elif key in (PSIBLAST_FLAGS) and value:
            cmd.append("-{}".format(key))

    cmd.append("-outfmt")
    cmd.append("6 {}".format(output_fields))
    return cmd

def build_blastn_command(query, db, evalue, kwargs, output_fields, out, blastn_path):
    """
    Generate a python subprocess command (list) for
    executing blastn.

    kwargs should be valid blastn options. 
    """

    cmd = [blastn_path, "-query", query, "-db", db, "-evalue", str(evalue), "-out", out]
    for key, value in kwargs.items():
        if key in BLASTN_OPTIONS:
            cmd.append("-{}".format(key))
            cmd.append(str(value))

    cmd.append("-outfmt")
    cmd.append("6 {}".format(output_fields))
    return cmd
