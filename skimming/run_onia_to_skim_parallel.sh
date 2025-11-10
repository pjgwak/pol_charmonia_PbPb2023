#!/usr/bin/env bash
set -e

J=4  # number of jobs in parallel

cat <<'LIST' | xargs -I{} -P "${J}" sh -c 'root -b -q "{}"'
'onia_to_skim_data.C(24)'
'onia_to_skim_data.C(25)'
'onia_to_skim_data.C(26)'
'onia_to_skim_mc.C(24, 0)'
'onia_to_skim_mc.C(25, 0)'
'onia_to_skim_mc.C(26, 0)'
'onia_to_skim_mc.C(24, 1)'
'onia_to_skim_mc.C(25, 1)'
'onia_to_skim_mc.C(26, 1)'
LIST