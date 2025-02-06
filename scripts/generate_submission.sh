#!/bin/bash

# Useful Commands:
# ------------------
# If you want to create multiple sub files at the same time for different samples you can do the following:

# year=2023BPix  # Set the year variable
# 
# for data_type in $(ls Run3_"$year"/*.txt | xargs -n 1 basename -s .txt); do
#     echo "Running command for: $data_type"
#     bash scripts/generate_submission.sh --year "$year" --data_type "$data_type"
# done

# ------------------
# Generate sub files: bash generate_submission.sh --year 2023BPix --data_type Muon1_v1
# Submit all sub files: find submissions -name "*.sub" | xargs -I {} condor_submit {}

# Default values
year=""
data_type=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --year) year="$2"; shift ;;   # Assign year from argument
        --data_type) data_type="$2"; shift ;; # Assign data_type from argument
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Check if required arguments are provided
if [[ -z "$year" || -z "$data_type" ]]; then
    echo "Usage: $0 --year <YEAR> --data_type <DATA_TYPE>"
    exit 1
fi

# Determine isMC based on data_type
if [ "$data_type" = "dy" ]; then
    isMC=1
else
    isMC=0
fi

# Create the submissions directory if it doesn't exist
mkdir -p submissions/Run3_${year}
mkdir -p logs/Run3_${year}/${data_type}
mkdir -p /eos/cms/store/group/phys_tau/ksavva/TauTrgSF/Run3_${year}/${data_type}

# Generate the submission file
cat <<EOF > submissions/Run3_${year}/job_submission_${data_type}.sub
executable              = scripts/run_job.sh
arguments               = --input \$(inputFile) --output /eos/cms/store/group/phys_tau/ksavva/TauTrgSF/Run3_${year}/${data_type}/ --isMC ${isMC} --era ${year}
log                     = logs/Run3_${year}/${data_type}/\$(ProcId).log
error                   = logs/Run3_${year}/${data_type}/\$(ProcId).err
output                  = logs/Run3_${year}/${data_type}/\$(ProcId).out

# Job runtime flavor
+JobFlavour             = "longlunch"

# Queue jobs from input file list
queue inputFile from Run3_${year}/${data_type}.txt
EOF

echo "Submission file 'submissions/Run3_${year}/job_submission_${data_type}.sub' has been generated."

