#!/bin/bash
abort() {
    echo "ERROR: $1"
    exit 1
}

${1:-.}/radial_interp >log.txt
grep 'The radial lines will be redistributed' log.txt >/dev/null || abort 'Lines were not redistributed'
grep 'Finished block   3:' log.txt >/dev/null || abort 'Program did not finish'
[[ -f forebody-fine.gu ]] || abort 'Grid file not generated'
