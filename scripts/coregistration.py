import os
import subprocess
import yaml
import sys

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
config_path = os.path.join(base_path, 'config.yaml')

def run_fsl_flirt(moving_image_path, fixed_image_path, output_image_path, transform_matrix_path):
   
    command = [
        'flirt',
        '-in', moving_image_path,
        '-ref', fixed_image_path,
        '-out', output_image_path,
        '-omat', transform_matrix_path
    ]
    try:
        subprocess.run(command, check=True)
        print(f"Coregistration completed successfully: {output_image_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error running FLIRT: {e}")
        sys.exit(1)

with open(config_path, 'r') as config_file:
    config = yaml.safe_load(config_file)

fixed_image_path = config['paths']['T1w'] 
moving_image_t1w_fs = config['paths']['T1w_fs'] 
moving_image_t2w = config['paths']['T2w'] 

input_dir = os.path.dirname(fixed_image_path)  

print("Fixed Image Path:", fixed_image_path)
print("Moving Image (T1w_fs) Path:", moving_image_t1w_fs)
print("Moving Image (T2w) Path:", moving_image_t2w)
print("Input Directory Path:", input_dir)

os.makedirs(input_dir, exist_ok=True)

output_image_t1w_fs = os.path.join(input_dir, 'r' + os.path.basename(moving_image_t1w_fs))
output_image_t2w = os.path.join(input_dir, 'r' + os.path.basename(moving_image_t2w))

transform_matrix_t1w_fs = os.path.join(input_dir, 'transform_matrix_T1w_fs.mat')
transform_matrix_t2w = os.path.join(input_dir, 'transform_matrix_T2w.mat')

# run FLIRT
run_fsl_flirt(moving_image_t1w_fs, fixed_image_path, output_image_t1w_fs, transform_matrix_t1w_fs)
run_fsl_flirt(moving_image_t2w, fixed_image_path, output_image_t2w, transform_matrix_t2w)

print("All coregistration tasks completed successfully.")
