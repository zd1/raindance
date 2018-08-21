'''
This is confirguration file
Specify file paths and parameters in this file
'''

param = {

    # reference fa file.
    "ref": "/file_path_to/hg19_GRCh37/hg19.fa",

    # list of sample. provided in the meta directory
    "samplelist" : "meta/samplelist",

    # fastqindex file has one file path per row. each path is a fastq of
    # one read, eigher 1 or 2, for one sample. Here's an example
    # /file_path_to/fastq/WTCHG_107218_04_1.fastq.gz
    # /file_path_to/fastq/WTCHG_107218_04_2.fastq.gz
    # /file_path_to/fastq/WTCHG_107218_06_1.fastq.gz
    # /file_path_to/fastq/WTCHG_107218_06_2.fastq.gz

    "fastqindex":"meta/fastq.index",

    "table" : "/file_path_to/meta/amplicons.csv",

    # output directory
    "fastq_outdir" : "/file_path_to/fastq.qc/",

    # parameters:
    "kmersize" : 10,
    "minQscore" : 20,  # PhredQ base quality score
    "lowqualbase" : 20, # read with more than this amount of bases below the minQscore above would be removed
    "readlength": 100,

    # fastq filters
    "minreadlength": 50,
    "probe_search_length": 25,
    "trim_read": 0,

    # Directory for the BAM files.
    "bamdir" : "/file_path_to/bams",

    # Output directories for pileups
    "pileup_outdir":"/file_path_to/pileups",
    "pileup_sum": "/file_path_to/sum/all.hd5",
    "pileup_agg": "/file_path_to/sum/",
    "vcf": "/file_path_to/sum//all.vcf.gz",
    "baseQ_cutoff": 0, # base quality cutoff
    "mapQ_cutoff": 0, # mapping quality cutoff
    "max_noncon_per_read": 10, # max amount of non concensus nucleotide in a read
    "min_coverage" : 5000, # min coverage
    "min_non_consensus_freq": 0.005, # counts of all non consensus bases together needs to reach this frequency to be considered as elevated.
    "min_het_frq" : 0.10, # considered a het if non consensus counts exceeds this threshold

    # call
    "call": "call.R",

    #alphabet
    "alphabet": ["A", "C", "G", "T"]
}
