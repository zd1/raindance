import os
import sys
import pdb
import argparse
import re
import subprocess


#short.qc:compC073:bioinf:zd1:rain:1560682:sge:"
HEADER = ["qname","hostname","group","owner","jobname","jobnumber","account",
          "priority","qsub_time","start_time","end_time", "failed","exit_status",
          "ru_wallclock","ru_utime","ru_stime","ru_maxrss",
          "ru_ixrss","ru_ismrss","ru_idrss","ru_isrss",
          "ru_minflt","ru_majflt","ru_nswap","ru_inblock","ru_oublock",
          "ru_nvcsw", "ru_nivcsw", "project", "department", "granted_pe","slots",
          "taskid","cpu", "mem","io", "cmd", "iow","arid", "maxvmem", 
          "ru_msgrcv","ru_nsignals"]

def load_log(jobname, jobid):
    #$SGE_ROOT/$SGE_CELL/common/accounting
    #accounting = os.environ['$SGE_ROOT'] + "/" + os.environ['$SGE_CELL'] +"/common/accounting"
    accounting = "/gpfs0/mgmt/sge/6.2u5/default/common/accounting"

    p1 = subprocess.Popen(["grep","%s:%s"%(jobname, jobid),accounting],stdout=subprocess.PIPE )
    targetjob = []
    while (p1.poll() == None):
        if p1.stdout == None:
            break
        for line in p1.stdout:
            d = line.strip().split(":")
            job = {HEADER[i]:d[i] for i in range(len(HEADER))}
            if job["exit_status"] == '137':
                targetjob.append(job['taskid'])
                print job
        if not line:
            break
    return targetjob

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Scan job status')
    parser.add_argument('--jobname', required = True,  default=None)
    parser.add_argument('--jobid', required = True,  default=None)
    parser.add_argument('--outfailedindex', required = True,  default=None)
    args = parser.parse_args()
    
    targetjob = load_log(args.jobname, args.jobid)
    ofh = open(args.outfailedindex, 'w')
    for job in targetjob:
        ofh.write(job+"\n")
    ofh.close()
    
