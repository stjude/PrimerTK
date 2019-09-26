#!/usr/bin/env python3

def primer3_input(header, sequence, args):
    """ Generates a primer3 input file with specified arguments.
    Args:
        header (string): The primer header string
        sequence (string): the sequence used to design primers
        args (Namespace): Arparse results.
    Returns:
        primer3_input_string (string): The string to write to the primer3 input file.
    """

    primer3_input_string = "SEQUENCE_ID="+header+'\n'\
    +"SEQUENCE_TEMPLATE="+sequence+'\n'\
    +"SEQUENCE_TARGET=%s" %args.sequence_target+'\n'\
    +"PRIMER_FIRST_BASE_INDEX=1"+'\n'\
    +"PRIMER_TASK=pick_detection_primers"+'\n'\
    +"PRIMER_MIN_THREE_PRIME_DISTANCE=3"+'\n'\
    +"PRIMER_MAX_LIBRARY_MISPRIMING=12.00"+'\n'\
    +"PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00"+'\n'\
    +"PRIMER_PRODUCT_SIZE_RANGE=%s" %args.product_size_range+'\n'\
    +"PRIMER_MAX_END_STABILITY=9.0"+'\n'\
    +"PRIMER_MAX_SELF_ANY_TH=45.00"+'\n'\
    +"PRIMER_MAX_SELF_END_TH=35.00"+'\n'\
    +"PRIMER_PAIR_MAX_COMPL_ANY_TH=45.00"+'\n'\
    +"PRIMER_PAIR_MAX_COMPL_END_TH=35.00"+'\n'\
    +"PRIMER_MAX_HAIRPIN_TH=24.00"+'\n'\
    +"PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00"+'\n'\
    +"PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00"+'\n'\
    +"PRIMER_TM_FORMULA=1"+'\n'\
    +"PRIMER_SALT_CORRECTIONS=1"+'\n'\
    +"PRIMER_SALT_MONOVALENT=50.0"+'\n'\
    +"PRIMER_INTERNAL_SALT_MONOVALENT=50.0"+'\n'\
    +"PRIMER_SALT_DIVALENT=1.5"+'\n'\
    +"PRIMER_INTERNAL_SALT_DIVALENT=1.5"+'\n'\
    +"PRIMER_DNTP_CONC=0.6"+'\n'\
    +"PRIMER_INTERNAL_DNTP_CONC=0.6"+'\n'\
    +"PRIMER_DNA_CONC=50.0"+'\n'\
    +"PRIMER_INTERNAL_DNA_CONC=50.0"+'\n'\
    +"PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1"+'\n'\
    +"PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=1"+'\n'\
    +"PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s" %args.thermopath+'\n'\
    +"PRIMER_PICK_LEFT_PRIMER=1"+'\n'\
    +"PRIMER_PICK_RIGHT_PRIMER=1"+'\n'\
    +"PRIMER_PICK_INTERNAL_OLIGO=1"+'\n'\
    +"PRIMER_MAX_POLY_X=3"+'\n'\
    +"PRIMER_LEFT_NUM_RETURNED=5"+'\n'\
    +"PRIMER_RIGHT_NUM_RETURNED=5"+'\n'\
    +"PRIMER_OPT_SIZE=%s" %args.primer_opt_size+'\n'\
    +"PRIMER_MIN_SIZE=%s" %args.primer_min_size+'\n'\
    +"PRIMER_MAX_SIZE=%s" %args.primer_max_size+'\n'\
    +"PRIMER_MIN_TM=%s" %args.primer_min_tm+'\n'\
    +"PRIMER_OPT_TM=%s" %args.primer_opt_tm+'\n'\
    +"PRIMER_MAX_TM=%s" %args.primer_max_tm+'\n'\
    +"PRIMER_MAX_NS_ACCEPTED=1"+'\n'\
    +"PRIMER_NUM_RETURN=5"+'\n'\
    +"P3_FILE_FLAG=1"+'\n'\
    +"PRIMER_EXPLAIN_FLAG=1"+'\n'\
    +"PRIMER_MISPRIMING_LIBRARY=%s" %args.mispriming+'\n'\
    +"PRIMER_MIN_GC=%s" %args.primer_min_gc+'\n'\
    +"PRIMER_OPT_GC_PERCENT=%s" %args.primer_opt_gc+'\n'\
    +"PRIMER_MAX_GC=%s" %args.primer_max_gc+'\n'\
    +"PRIMER_PAIR_MAX_DIFF_TM=3"+'\n'\
    +"="+'\n'
    return primer3_input_string
