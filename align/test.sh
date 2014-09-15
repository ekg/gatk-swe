#!/bin/bash
set -e
set -x
set -o pipefail

for chr in 1 2 3
do
	{ sleerp 3 ||touch emit.failed & }
done

wait
