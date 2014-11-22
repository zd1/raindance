'''This is confirguration file
'''

#basedir = "/Users/zd1/Work/projects2/raindance"
#basedir = "/users/mcvean/zd1/volumn/raindance"
basedir = "/home/crangen/zding/projects/raindance"

param = {
    "ref": "/Users/zd1/mount/wimmhts/common/ref/old/hg19_GRCh37/hg19.fa",
    "samplelist" : "/users/mcvean/zd1/volumn/raindance/meta/samplelist",
    # "ref" : "/users/mcvean/zd1/volumn/common/ref/old/hg19_GRCh37/hg19.fa",
#    "fastqindex":"%s/fastq.index"%basedir,
    "fastqindex":"%s/fastq.lib1.index"%basedir,
    "table" : "%s/amplicons.csv"%basedir, 
    "outdir" : "%s/fastq.qc"%basedir,
    "imageoutdir" : "%s/pileups_images"%basedir, 
    # parameters:
    "kmersize" : 10,
    "lowqualbase" : 20, # reads with more than this amount of bad bases are considered low quality
    "minQscore" : 20, # reads with more than this amount of bad bases are considered low quality
    "readlength": 100,
    
    # fastq filters
    "minreadlength": 50,
    "probe_search_length": 25,
    "trim_read": 0,

    # pileup filters
    "baseQ_cutoff": 0, # use same cutoff as the one used to process fastq file
    "mapQ_cutoff": 0, # 0-255, -10log10{mapping position is wrong}
    "max_noncon_per_read": 5, # 0-255, -10log10{mapping position is wrong}
    
    # pileup exporting
    "min_coverage" : 5000, # min coverage
    "min_non_consensus_freq": 0.005, # counts of all other bases together need to reach this frequency to be considered as elevated. 
    "min_het_frq" : 0.10, # considered a het if non consensus counts exceed this threshold
    
    # scripts
    "aligner":"/home/crangen/zding/workspace/rain/call.sh"
}

