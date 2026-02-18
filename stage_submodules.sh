#!/bin/bash
#
# Script to stage the full contents of a repo with submodules.
# The staged tree can later be tarred, or rsync'ed to another computer.
#
# Usage:
#
#   STAGE_ROOT=/path/to/staging/area sh stage_submodules.sh
#
# The main repo's contents will be those of its index. Make sure you know
# what you are doing.
#
# You need to be at the root of the super-project when executing this script.
#
# Make sure that your submodules are 'updated' (pointed to the right checkout)
#
# This is the path at which the source will be put
#
STAGE_ROOT=${STAGE_ROOT:-/tmp/_stage_root}
#
echo "Staging sources (incl submodules) to: $STAGE_ROOT"
echo "You can later tar or rsync them for use on a different computer"
#
# NOTE the trailing slash !
git checkout-index -a -f --prefix=${STAGE_ROOT}/
#
# Note that $displaypath contains the relative path
# from the current working directory to the submodules root directory. This
# is essential for proper placement of submodules of a submodule, enabled
# by the --recursive option
# The escape before $displaypath avoids its evaluation by the outer shell
#
git submodule foreach --recursive \
       "git checkout-index -a -f --prefix=${STAGE_ROOT}/\$displaypath/"

#
# Generate a SIESTA.release file 
#
git describe --first-parent --tags > ${STAGE_ROOT}/SIESTA.release
