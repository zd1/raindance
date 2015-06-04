'''
This is confirguration file
Specify file paths and parameters in this file
'''

param = {
    "ref": "/Users/zd1/volumn/wt/common/ref/old/hg19_GRCh37/hg19.fa",
    "samplelist" : "/Users/zd1/volumn/wt/raindance/meta/samplelist",
    "fastqindex":"/Users/zd1/volumn/wt/raindance/meta/fastq.index",
    "table" : "/Users/zd1/volumn/wt/raindance/meta/amplicons.csv", 
    "ampshortlist": "/users/mcvean/zd1/workspace/raindance/amp_shortlist",
    "fastq_outdir" : "/Users/zd1/volumn/wt/raindance/fastq.qc/",
    # parameters:
    "kmersize" : 10,
    "lowqualbase" : 20, # a read with more than this amount of bases below minQscore is considered as bad quality
    "minQscore" : 20,  # PhredQ base quality score
    "readlength": 100,
    
    # fastq filters
    "minreadlength": 50,
    "probe_search_length": 25,
    "trim_read": 0,

    # BAM directory, where each sample is stored bamdir/sample/sample.sorted.bam
    "bamdir" : "/Users/zd1/volumn/wimm/raindance/bams",
    
    # pileup output dirs
    "pileup_outdir":"/Users/zd1/volumn/wt/raindance/pileups",
        
    # pileup filters
    "baseQ_cutoff": 0, # base quality cutoff
    "mapQ_cutoff": 0, # mapping quality cutoff
    "max_noncon_per_read": 10, # max amount of non concensus nucleotide in a read
    
    # pileup exporting
    "min_coverage" : 5000, # min coverage
    "min_non_consensus_freq": 0.005, # counts of all non consensus bases together needs to reach this frequency to be considered as elevated. 
    "min_het_frq" : 0.10, # considered a het if non consensus counts exceeds this threshold
}
