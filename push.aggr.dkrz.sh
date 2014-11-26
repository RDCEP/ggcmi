#!/bin/bash

rsync --omit-dir-times --delete -av -e 'ssh -i ~/.ssh/id_rsa' processed b380079@vre1.dkrz.de:/gpfs_750/projects/ISI_MIP/data/AgMIP.output
