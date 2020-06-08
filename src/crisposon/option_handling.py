from crisposon.blast_options import (BLAST_OPTIONS_COMMON,
                                        BLASTP_OPTIONS,
                                        PSIBLAST_OPTIONS,
                                        BLAST_FLAGS)

def build_blastp_command(query, db, evalue, kwargs, output_fields, out):
    """
    Generate a python subprocess command (list) for
    executing blastp.

    kwargs should be valid blastp options. Note that
    at this time, only certain options are allowed.
    """

    cmd = ["blastp", "-query", query, "-db", db, "-evalue", evalue, "-out", out]
    for key, value in kwargs.items():
        if key in (BLAST_OPTIONS_COMMON + BLASTP_OPTIONS):
            cmd.append("-{}".format(key))
            cmd.append(value)
        elif key in (BLAST_FLAGS) and value:
            cmd.append("--{}".format(key))

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

    cmd = ["psiblast", "-query", query, "-db", db, "-evalue", evalue, "-out", out]
    for key, value in kwargs.items():
        if key in (BLAST_OPTIONS_COMMON + BLASTP_OPTIONS):
            cmd.append("-{}".format(key))
            cmd.append(value)
        elif key in (BLAST_FLAGS) and value:
            cmd.append("--{}".format(key))

    cmd.append("-outfmt")
    cmd.append("6 {}".format(output_fields))
    return cmd
    