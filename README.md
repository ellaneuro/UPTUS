# UPTUS
Transcranial Ultrasound Stimulation Planning 

# Requirements
- MATLAB
- Freesurfer 
- FSL
- Connectome Workbench
- SimNIBS for Python (simnibs)

# Inputs
Input variables in the config.yaml file that is within the UPTUS main folder. Make sure to name the final version with your updated paths as 'config.yaml'. 

# Running the Scripts

You can run this project in various environments:
- Linux/Mac: Open a terminal and navigate to the project directory.
- Windows: Download FSL, choose Linux. 

# Commands

1. Navigate to the UPTUS directory. 
2. make the zshell script executable by running "chmod +x main.sh"  
3. run the wrapper script with "./main.sh"


# Jupyter Notebook
To run the scripts in a Jupyter notebook, use the `!` command prefix to execute shell commands:
```python
1. "!chmod +x main.sh"
2. "!./main.sh"

# User input
The script will ask for COM coordinates from target and coordinates from the entry point

# The output found in the UPTUS/output folder will contain 
- sub-001.msh
- roi_mesh.nii.gz
- roi_coordinates.txt
- pCT.nii

# choose an entry point

# input entry point coordinates (RAS) when prompted

# prep is done - next step is K-wave simulation 
