#!/bin/bash -x

# [[file:~/dev/org/pipe.org::*MIP%20analysis][rsync_miseq]]

rsync -avz --rsh='ssh -p5111 -l zding' \
      --exclude=".git/*" \
      --exclude="sync.sh" \
      ./ localhost:/home/crangen/zding/projects/raindance/mip/code/

# rsync_miseq ends here
