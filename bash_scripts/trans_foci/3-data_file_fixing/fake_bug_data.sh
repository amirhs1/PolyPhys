#!/bin/bash
# fix the 'all' data files with coreecting the number of crowders.
for alldata in eps*.all.data; do
    sim=$(echo "${alldata}" | cut -d . -f -6)
    echo "${sim}"
    echo cp fake_nc0_all.data "${sim}".bug.data
done