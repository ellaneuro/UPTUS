#!/bin/bash

# If using K-Plan, run the kplan_prep.m script on all three images

SUBJECT_ID=$(python -c "
import yaml
with open('config.yaml') as f:
    config = yaml.safe_load(f)
print(config.get('subject_id', 'default_subject'))
")

### COREGISTRATION ####

echo "Activating coregistration environment..."
if conda env list | grep -q "coreg_env"; then
    conda activate coreg_env > /dev/null 2>&1
else
    conda create -n coreg_env python=3.10 -y > /dev/null 2>&1
    conda activate coreg_env > /dev/null 2>&1
fi

echo "Installing dependencies from coregistration_requirements.txt..."
pip install --upgrade pip > /dev/null 2>&1  
pip install -r requirements/coregistration_requirements.txt > /dev/null 2>&1

echo "Running coregistration script..."
python scripts/coregistration.py
COREG_STATUS=$?

conda deactivate > /dev/null 2>&1

if [ $COREG_STATUS -ne 0 ]; then
    echo "Coregistration script failed with status $COREG_STATUS."
    exit $COREG_STATUS
fi

# Get file paths from config
FIXED_IMAGE=$(python -c "
import yaml
with open('config.yaml') as f:
    config = yaml.safe_load(f)
print(config['paths']['T1w'])
")

T1W_FS_OUTPUT=$(python -c "
import yaml, os
with open('config.yaml') as f:
    config = yaml.safe_load(f)
input_dir = os.path.dirname(config['paths']['T1w'])
fs_basename = os.path.basename(config['paths']['T1w_fs'])
if not fs_basename.endswith('.gz'):
    fs_basename += '.gz'
print(os.path.join(input_dir, 'r' + fs_basename))
")

T2W_OUTPUT=$(python -c "
import yaml, os
with open('config.yaml') as f:
    config = yaml.safe_load(f)
input_dir = os.path.dirname(config['paths']['T1w'])
t2w_basename = os.path.basename(config['paths']['T2w'])
if not t2w_basename.endswith('.gz'):
    t2w_basename += '.gz'
print(os.path.join(input_dir, 'r' + t2w_basename))
")

echo "Opening Freeview for coregistration quality check..."
freeview -v \
    "$FIXED_IMAGE" \
    "$T1W_FS_OUTPUT":colormap=heat \
    "$T2W_OUTPUT":colormap=jet &
FREEVIEW_PID=$!

echo "Freeview is running. Please evaluate the coregistration quality and respond to the prompt below."

while true; do
    read -p "Is the coregistration quality acceptable? (yes/no): " yn
    case $yn in
        [Yy]* ) 
            echo "Proceeding with the pipeline..."
            break;;
        [Nn]* ) 
            echo "Pipeline stopped due to unacceptable coregistration quality."
            kill $FREEVIEW_PID 2>/dev/null
            exit 1;;
        * ) echo "Please answer yes or no.";;
    esac
done

kill -9 $FREEVIEW_PID 2>/dev/null

### FREESURFER RECON -ALL ###

echo "Activating FreeSurfer environment..."
if conda env list | grep -q "freesurfer_env"; then
    conda activate freesurfer_env > /dev/null 2>&1
else
    conda create -n freesurfer_env python=3.10 -y > /dev/null 2>&1
    conda activate freesurfer_env > /dev/null 2>&1
fi

FREESURFER_SETUP_SCRIPT=$(python -c "
import yaml
import sys

with open('config.yaml') as f:
    config = yaml.safe_load(f)

if 'freesurfer' not in config or 'freesurfer_setup_script' not in config['freesurfer']:
    print('Error: freesurfer_setup_script not found in config.yaml', file=sys.stderr)
    sys.exit(1)

print(config['freesurfer']['freesurfer_setup_script'])
")

echo "Sourcing FreeSurfer setup script..."
source "$FREESURFER_SETUP_SCRIPT"

echo "Running FreeSurfer script..."
python scripts/freesurfer_orig.py &
FREESURFER_PID=$!

conda deactivate > /dev/null 2>&1

### MR to PCT ###

echo "Activating MR to pCT environment..."
if conda env list | grep -q "mr_to_pct_env"; then
    conda activate mr_to_pct_env > /dev/null 2>&1
else
    conda create -n mr_to_pct_env python=3.10 -y > /dev/null 2>&1
    conda activate mr_to_pct_env > /dev/null 2>&1
fi

echo "Installing dependencies from mr_to_pct_requirements.txt..."
pip install --upgrade pip > /dev/null 2>&1  
pip install monai==1.2.0 torch antspyx==0.4.2 > /dev/null 2>&1

echo "Verifying PyTorch installation and device..."
python -c "
import torch
try:
    print(f'PyTorch version: {torch.__version__}')
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Using device: {device}')
except ImportError:
    print('Error: PyTorch is not installed correctly.')
    exit(1)
" || {
    echo "Error: PyTorch or device configuration failed."
    exit 1
}

PRETRAINED_WEIGHTS="requirements/pretrained_net_final_20220825.pth"

if [ ! -f "$PRETRAINED_WEIGHTS" ]; then
    echo "Error: Pretrained weights not found in $PRETRAINED_WEIGHTS."
    exit 1
fi

INPUT_MR_FILE=$(python -c "
import yaml
import sys

with open('config.yaml') as f:
    config = yaml.safe_load(f)

if 'paths' not in config or 'T1w' not in config['paths']:
    print('Error: T1w path not found in config.yaml', file=sys.stderr)
    sys.exit(1)

print(config['paths']['T1w'])
")

OUTPUT_DIR=$(python -c "
import yaml
import sys

with open('config.yaml') as f:
    config = yaml.safe_load(f)

if 'paths' not in config or 'output_directory' not in config['paths']:
    print('Error: output_directory not found in config.yaml', file=sys.stderr)
    sys.exit(1)

print(config['paths']['output_directory'])
")

OUTPUT_PCT_FILE="${OUTPUT_DIR}/${SUBJECT_ID}_pct.nii"

echo "Running MR to pCT script..."
python "mr-to-pct-main/mr-to-pct.py" "$INPUT_MR_FILE" "$OUTPUT_PCT_FILE" "$PRETRAINED_WEIGHTS" &
MR_TO_PCT_PID=$!

conda deactivate > /dev/null 2>&1

### SIMNIBS CHARM SEGMENTATION ###

echo "Activating SimNIBS environment..."
if conda env list | grep -q "simnibscharm_env"; then
    conda activate simnibscharm_env > /dev/null 2>&1
else
    conda create -n simnibscharm_env python=3.9 -y > /dev/null 2>&1
    conda activate simnibscharm_env > /dev/null 2>&1
fi

echo "Installing dependencies from simnibs_requirements.txt..."
pip install --upgrade pip > /dev/null 2>&1 
pip install -r requirements/simnibs_requirements.txt > /dev/null 2>&1

echo "Running SimNIBS script..."
python "scripts/segmentation.py" &
SIMNIBS_PID=$!

wait $FREESURFER_PID
FREESURFER_STATUS=$?
wait $MR_TO_PCT_PID
MR_TO_PCT_STATUS=$?
wait $SIMNIBS_PID
SIMNIBS_STATUS=$?

if [ $FREESURFER_STATUS -ne 0 ]; then
    echo "FreeSurfer script failed with status $FREESURFER_STATUS."
    exit $FREESURFER_STATUS
fi

if [ $MR_TO_PCT_STATUS -ne 0 ]; then
    echo "MR to PCT script failed with status $MR_TO_PCT_STATUS."
    exit $MR_TO_PCT_STATUS
fi

if [ $SIMNIBS_STATUS -ne 0 ]; then
    echo "SimNIBS script failed with status $SIMNIBS_STATUS."
    exit $SIMNIBS_STATUS
fi

GREEN=$(tput setaf 2)
RED=$(tput setaf 1)
RESET=$(tput sgr0)

echo "----- Process Summary -----"
if [ $FREESURFER_STATUS -eq 0 ] && [ $MR_TO_PCT_STATUS -eq 0 ] && [ $SIMNIBS_STATUS -eq 0 ]; then
    echo "${GREEN}All scripts ran successfully.${RESET}"
else
    echo "${RED}One or more scripts failed:${RESET}"
    [ $FREESURFER_STATUS -ne 0 ] && echo " - FreeSurfer script failed with status $FREESURFER_STATUS." || echo " - FreeSurfer script ran successfully."
    [ $MR_TO_PCT_STATUS -ne 0 ] && echo " - MR to PCT script failed with status $MR_TO_PCT_STATUS." || echo " - MR to PCT script ran successfully."
    [ $SIMNIBS_STATUS -ne 0 ] && echo " - SimNIBS script failed with status $SIMNIBS_STATUS." || echo " - SimNIBS script ran successfully."
fi

### TRANSDUCER PLACEMENT ### 

# echo "Activating TransducerPlacement environment..."
# if conda env list | grep -q "transducerplacement_env"; then
    #conda activate transducerplacement_env > /dev/null 2>&1
# else
    #conda create -n transducerplacement_env python=3.9 -y > /dev/null 2>&1
    #conda activate transducerplacement_env > /dev/null 2>&1

#fi

#echo "Installing dependencies from transducerplacement_requirements.txt..."
#pip install --upgrade pip > /dev/null 2>&1
#pip install -r requirements/transducerplacement_requirements.txt > /dev/null 2>&1

#echo "Running Transducerplacement.py..."
#python "scripts/transducerplacement.py" 

### RAS TO VOXEL ###

#echo "Extracting and converting RAS coordinates..."
#python "scripts/coordinates.py"

### K-WAVE ###









