#!/bin/bash
set -e  # Exit on errors
set -x  # Debug mode: print commands as they are executed

cd /afs/cern.ch/work/k/ksavva/public/TauTrgSF/NanoAODTools/ 
source PhysicsTools/NanoAODTools/standalone/env_standalone.sh
cd /afs/cern.ch/work/k/ksavva/public/TauTrgSF/NanoAODTools/PhysicsTools/Tau-Trigger/

# Run the Python post-processing script
python3 scripts/nano_postproc.py "$@"

