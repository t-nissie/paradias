#!/usr/bin/csh
# Usage: submit EM general -jobtype parallel -n 4 -B -N -u t-nissie@imr.tohoku.ac.jp -exec coordinates_check.csh -J coordinates_check
##
setenv MEMORY_AFFINITY MCM
setenv MALLOCMULTIHEAP true

poe ./coordinates_check
