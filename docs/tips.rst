Inputs and Outputs
==================

.. _building-sequence-databases:

Building sequence databases
---------------------------

To search for gene clusters with Opfi, users must compile representative protein (or nucleic acid) sequences for any genes expected in target clusters (or for any non-essential accessory genes of interest). These may be from a pre-existing, private collection of sequences (perhaps from a previous bioinformatics analysis). Alternatively, users may download sequences from a publically available database such as `Uniprot <https://www.uniprot.org/>`_ (maintained by the `European Bioinformatics Institute <https://www.ebi.ac.uk/>`_ ) or one of the `databases <https://www.ncbi.nlm.nih.gov/>`_ provided by the National Center for Biotechnology Information. 

Once target sequences have been compiled, they must be converted to an application-specific database format. Opfi currently supports :program:`BLAST+`, :program:`mmseqs2`, and :program:`diamond` for homology searching:

* `Instructions for creating sequence databases for BLAST using makeblastdb <https://www.ncbi.nlm.nih.gov/books/NBK569841/>`_
* `Instructions for creating sequence databases for mmseqs2 using mmseqs createdb <https://github.com/soedinglab/mmseqs2/wiki#searching>`_
* `Diamond makedb command options <https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#makedb-options>`_

The FASTA file format
#####################

Both genomic input data and reference sequence data should be in `FASTA <https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp>`_ format. This is a simple flat text representation of biological sequence data, where individual sequences are delineated by the ``>`` greater than character. For example:

.. code-block:: 

    >UniRef50_Q02ML7 CRISPR-associated endonuclease Cas1 n=1700 RepID=CAS1_PSEAB
    MDDISPSELKTILHSKRANLYYLQHCRVLVNGGRVEYVTDEGRHSHYWNIPIANTTSLLL
    GTGTSITQAAMRELARAGVLVGFCGGGGTPLFSANEVDVEVSWLTPQSEYRPTEYLQRWV
    GFWFDEEKRLVAARHFQRARLERIRHSWLEDRVLRDAGFAVDATALAVAVEDSARALEQA
    PNHEHLLTEEARLSKRLFKLAAQATRYGEFVRAKRGSGGDPANRFLDHGNYLAYGLAATA
    TWVLGIPHGLAVLHGKTRRGGLVFDVADLIKDSLILPQAFLSAMRGDEEQDFRQACLDNL
    SRAQALDFMIDTLKDVAQRSTVSA
    >UniRef50_Q2RY21 CRISPR-associated endonuclease Cas1 1 n=1034 RepID=CAS1A_RHORT
    MADPAFVPLRPIAIKDRSSIVFLQRGQLDVVDGAFVLIDQEGVRVQIPVGGLACLMLEPG
    TRITHAAIVLCARVGCLVIWVGERGTRLYAAGQPGGARADRLLFQARNALDETARLNVVR
    EMYRRRFDDDPPARRSVDQLRGMEGVRVREIYRLLAKKYAVDWNARRYDHNDWDGADIPN
    RCLSAATACLYGLCEAAILAAGYAPAIGFLHRGKPQSFVYDVADLYKVETVVPTAFSIAA
    KIAAGKGDDSPPERQVRIACRDQFRKSGLLEKIIPDIEEILRAGGLEPPLDAPEAVDPVI
    PPEEPSGDDGHRG

The sequence definition (defline) comes directly after the ``>`` character, and should be on a separate line from the sequence (which can be on one or more subsequent lines). There is no specific defline format, however, Opfi requires that, for both genomic input and sequence data, each definition line contain a unique sequence identifer. This should be a single word/token immediately following the ``>`` character (i.e. spaces between the ``>`` character and the identifier are not allowed). Any additional text on the defline is parsed as a single string, and appears in the output CSV (see :ref:`opfi-output-format`).

.. tip::

    Biological sequences downloaded from most public databases will have an accession number/identifier by default.

.. _labeling-sequences:

Annotating sequence databases
#############################

To take full advantage of the rule-based filtering methods in :mod:`operon_analyzer.rules`, users are encouraged to annotate reference sequences with a name/label that is easily searched. Labels can be as broad or as specific as is necessary to provide meaningful annotation of target gene clusters.

Gene labels are parsed from sequence deflines; specifically, Opfi looks for the second word/token following the ``>`` character. For example, the following FASTA sequence has been annotated with the label "cas1":

.. code-block:: 

    >UniRef50_Q02ML7 cas1 CRISPR-associated endonuclease Cas1 n=1700 RepID=CAS1_PSEAB
    MDDISPSELKTILHSKRANLYYLQHCRVLVNGGRVEYVTDEGRHSHYWNIPIANTTSLLL
    GTGTSITQAAMRELARAGVLVGFCGGGGTPLFSANEVDVEVSWLTPQSEYRPTEYLQRWV
    GFWFDEEKRLVAARHFQRARLERIRHSWLEDRVLRDAGFAVDATALAVAVEDSARALEQA
    PNHEHLLTEEARLSKRLFKLAAQATRYGEFVRAKRGSGGDPANRFLDHGNYLAYGLAATA
    TWVLGIPHGLAVLHGKTRRGGLVFDVADLIKDSLILPQAFLSAMRGDEEQDFRQACLDNL
    SRAQALDFMIDTLKDVAQRSTVSA

After running :class:`gene_finder.pipeline.Pipeline`, users could select candidates with hits against this sequence using the following rule set:

.. code-block:: Python

    from operon_analyzer.rules import RuleSet

    rs = RuleSet.require("cas1")

In practice, a genomics search might use a reference database of hundreds (or even thousands) of representative protein sequences, in which case labeling each sequence individually would be tedious. It is recommended to organize sequences into groups of related proteins that can be given a single label. This script uses the Python package :program:`Biopython` to annotate sequences in a multi-sequence FASTA file:

.. code-block:: Python 

    from Bio import SeqIO
    import os, sys

    def annotate_reference(prot_ref_file, label):
        records = list(SeqIO.parse(ref_fasta, "fasta"))
            
        for record in records:
            des = record.description.split()
            prot_id = des.pop(0)
            des_with_label = "{} {} {}".format(prot_id, label, " ".join(des))
            record.description = des_with_label

        SeqIO.write(records, ref_fasta, "fasta")

    if __name__ == "__main__":
        ref_fasta = sys.argv[1]
        label = sys.argv[2]
        annotate_reference(ref_fasta, label)

It is possible to use the entire sequence description (i.e. all text following the sequence identifier) as the gene label. This is particularly useful when using a pre-built database like `nr <https://www.ncbi.nlm.nih.gov/refseq/about/nonredundantproteins/>`_, which contains representative protein sequences for many different protein families. When using sequence databases that haven't been annotated, users should set ``parse_descriptions=False`` for each :class:`gene_finder.pipeline.Pipeline` ``add_step()`` method call.

Converting sequence files to a sequence database
################################################

Once reference sequences have been compiled (and, optionally, labeled) they must be converted to a sequence database format that is specific to the homology search program used. Currently, Opfi supports :program:`BLAST`, :program:`mmseqs2`, and :program:`diamond`. Each software package is automatically installed with a companion utility program for generating sequence databases. The following example shows what a typical call to :program:`makeblastdb`, the BLAST+ database utility program, might look like:

.. code-block:: bash 

    makeblastdb -in "my_sequences.fasta" -out my_sequences/db -dbtype prot -title "my_sequences" -hash_index

The command takes a text/FASTA file ``my_sequences.fasta`` as input, and writes the resulting database files to the directory ``my_sequences``. Database files are prefixed with "db". ``-dbtype prot`` specifies that the input is amino acid sequences. We use ``-title`` to name the database (required by BLAST). ``-hash_index`` directs :program:`makeblastdb` to generate a hash index of protein sequences, which can speed up computation time.

.. tip::

    :program:`mmseqs2` and :program:`diamond` have similar database creation commands, see :ref:`building-sequence-databases`. 

BLAST advanced options
----------------------

BLAST+ programs have a number of tunable parameters that can, for example, be used to adjust the sensitivity of the search algorithm. We anticipate that application defaults will be sufficient for most users; nevertheless, it is possible to use non-default program options by passing them as keyword arguments to :class:`gene_finder.pipeline.Pipeline` ``add_step()`` methods. 

For example, when using :program:`blastp` on the command line, we could adjust the number of CPUs to four by passing the argument ``-num_threads 4`` to the program. When using Opfi, this would look like ``num_threads=4``. 

Flags (boolean arguments that generally do not precede additional data) are also possible. For example, the command line flag ``-use_sw_tback`` tells :program:`blastp` to compute locally optimal Smith-Waterman alignments. The correct way to specify this behavior via the :class:`gene_finder.pipeline.Pipeline` API would be to use the argument ``use_sw_tback=True``. 

Below is a list of options accepted by Opfi. Note that some BLAST+ options are not allowed, mainly those that modify BLAST output.

.. csv-table::
    :header: "Program", "Allowed Options"

    ":program:`blastp` and :program:`psiblast`", "dbsize word_size gapopen gapextend qcov_hsp_perc xdrop_ungap xdrop_gap xdrop_gap_final searchsp sum_stats seg soft_masking matrix threshold culling_limit window_size num_threads comp_based_stats gilist seqidlist negative_gilistdb_soft_mask db_hard_mask entrez_query max_hspsbest_hit_overhang best_hit_score_edge max_target_seqsimport_search_strategy export_search_strategy num_alignments"
    ":program:`blastp` only", "task"
    ":program:`psiblast` only", "gap_trigger num_iterations out_pssm out_ascii_pssm pseudocount inclusion_ethresh"
    ":program:`blastp` (flags)", "lcase_masking ungapped use_sw_tback remote"
    ":program:`psiblast` (flags)", "lcase_masking use_sw_tback save_pssm_after_last_round save_each_pssm remote"
    ":program:`blastn`", "filtering_algorithm sum_stats window_masker_db window_size template_type version parse_deflines min_raw_gapped_score string format max_hsps taxids negative_taxids num_alignments strand off_diagonal_range subject_besthit num_sequences no_greedy negative_taxidlist culling_limit xdrop_ungap open_penalty DUST_options sorthits xdrop_gap_final negative_gilist subject use_index bool_value filename seqidlist task_name sort_hits database_name lcase_masking query_loc subject_loc sort_hsps line_length boolean db_hard_mask negative_seqidlist template_length filtering_db filtering_database penalty searchsp ungapped type gapextend db_soft_mask dbsize qcov_hsp_perc sorthsps window_masker_taxid index_name export_search_strategy float_value soft_masking gilist entrez_query show_gis best_hit_score_edge gapopen subject_input_file range html word_size best_hit_overhang perc_identity input_file num_descriptions xdrop_gap dust taxidlist max_target_seqs num_threads task remote int_value extend_penalty reward import_search_strategy num_letters"

You can read more about BLAST+ options in the `BLAST+ appendices <https://www.ncbi.nlm.nih.gov/books/NBK279684/>`_. 

.. note::

    Using advanced options with :program:`mmseqs2` and :program:`diamond` is not supported at this time. 

.. _opfi-output-format:

Opfi output format
------------------

Results from :class:`gene_finder.pipeline.Pipeline` searches are written to a single CSV file. Below is an example from the tutorial (see :ref:`example-usage`):

.. csv-table::
    :file: csv/example_output.csv
    :header-rows: 0

The first two columns contain the input genome/contig sequence ID (sometimes called an accession number) and the coordinates of the candidate gene cluster, respectively. Since an input file can have multiple genomic sequences, these two fields together uniquely specify a candidate gene cluster. Each row represents a single annotated feature in the candidate locus. Features from the same candidate are always grouped together in the CSV. 

Descriptions of each output field are provided below. Alignment statistic naming conventions are from the BLAST documentation, see `BLAST+ appendices <https://www.ncbi.nlm.nih.gov/books/NBK279684/>`_ (specifically "outfmt" in table C1). This `glossary <https://www.ncbi.nlm.nih.gov/books/NBK62051/>`_ of common BLAST terms may also be useful in interpreting alignment statistic meaning. 

.. csv-table::
    :file: csv/fieldnames.csv
    :header-rows: 1
