#!/bin/bash -x

# [[file:~/dev/org/pipe.org::*MIP%20analysis][rsync_miseq]]

rsync -avz --rsh='ssh -p5111 -l zding' --exclude=".git/*" ./ localhost:/home/crangen/zding/projects/raindance/mip/scan/

# rsync_miseq ends here
