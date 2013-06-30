#!/bin/bash
for f in *qe ; do tail --lines=+3 $f > `basename $f .qe`; done
# Copies file.qe to file, leaving out first 2 lines
# echo $PATH to see where to put, i.e. /usr/local/bin
