import yaml
import os
import shutil
import subprocess

def load_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def find_coregistered_image(image_path):

    variants = [image_path, image_path + '.gz'] 
    for variant in variants:
        if os.path.exists(variant):
            return variant
    raise FileNotFoundError(f"No coregistered image found for: {image_path}")

def run_simnibs_charm(coreg_t1, coreg_t2, simnibs_dir, output_dir, subject_id):
    
    # Define the command for running SimNIBS CHARM
    charm_command = ['charm', subject_id, coreg_t1, coreg_t2]
    
    # Run the CHARM command
    subprocess.run(charm_command, cwd=simnibs_dir, check=True)
    
    # Define the path to the subject's m2m directory
    subject_m2m_dir = os.path.join(simnibs_dir, f'm2m_{subject_id}')
    
    # Check if the m2m directory was created
    if os.path.exists(subject_m2m_dir):
        # Define the path to the output directory for this subject
        subject_output_dir = os.path.join(output_dir, f'm2m_{subject_id}')
        
        # Copy the entire m2m folder to the output directory
        shutil.copytree(subject_m2m_dir, subject_output_dir, dirs_exist_ok=True)
        print(f"Copied entire m2m folder to: {subject_output_dir}")
    else:
        print(f"Warning: m2m directory for {subject_id} not found in the expected output directory.")

# Load configuration
config_path = 'config.yaml'
config = load_config(config_path)

# Define paths from config
t1w_fs = config['paths']['T1w_fs']
t2w = config['paths']['T2w']
simnibs_dir = config['paths']['simnibs_directory']
output_dir = config['paths']['output_directory']
subject_id = config['subject_id']

# Construct coregistered image paths
coreg_t1 = os.path.join(os.path.dirname(t1w_fs), 'r' + os.path.basename(t1w_fs))
coreg_t2 = os.path.join(os.path.dirname(t2w), 'r' + os.path.basename(t2w))

# Verify coregistered images
coreg_t1 = find_coregistered_image(coreg_t1)
coreg_t2 = find_coregistered_image(coreg_t2)

# Run CHARM and copy entire m2m folder
run_simnibs_charm(coreg_t1, coreg_t2, simnibs_dir, output_dir, subject_id)
