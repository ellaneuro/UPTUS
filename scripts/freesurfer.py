import yaml
import os
import subprocess
import numpy as np
import nibabel as nib

def load_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def run_freesurfer_recon(subject_id, input_file, fs_dir):
    subprocess.run(['recon-all', '-s', subject_id, '-i', input_file, '-all'], cwd=fs_dir, check=True)

def save_roi_mask(subject_id, fs_dir, roi_code, output_directory, roi_name):
    roi_mask_path = os.path.join(fs_dir, subject_id, 'mri', 'aparc+aseg.mgz')
    roi_mask_img = nib.load(roi_mask_path)

    roi_mask_data = roi_mask_img.get_fdata()
    roi_binary_mask = np.zeros_like(roi_mask_data)
    roi_binary_mask[roi_mask_data == roi_code] = 1  

    roi_mask_nifti = nib.Nifti1Image(roi_binary_mask, affine=roi_mask_img.affine)
    roi_mask_output_path = os.path.join(output_directory, f'{roi_name}_mask_T1w.nii.gz')
    nib.save(roi_mask_nifti, roi_mask_output_path)
    print(f"ROI mask saved to: {roi_mask_output_path}")

def calculate_center_of_mass(roi_name, subject_id, output_file, fs_dir, roi_code):
    roi_mask_path = os.path.join(fs_dir, subject_id, 'mri', 'aparc+aseg.mgz')
    roi_mask_img = nib.load(roi_mask_path)
    roi_mask_data = roi_mask_img.get_fdata()

    roi_indices = np.argwhere(roi_mask_data == roi_code)

    if roi_indices.size > 0:
        com_voxel = np.mean(roi_indices, axis=0)

        # Apply affine transformation to get coordinates in T1w space
        com_t1w = nib.affines.apply_affine(roi_mask_img.affine, com_voxel)

        ras_affine_path = os.path.join(fs_dir, subject_id, 'mri', 'transforms', 'talairach.xfm')
        if os.path.exists(ras_affine_path):
            ras_affine = []
            with open(ras_affine_path, 'r') as file:
                for line in file:
                    if line.startswith(('MNI Transform File', '%', 'Transform_Type', 'Linear_Transform')):
                        continue
                    if line.strip(): 
                        cleaned_line = line.strip().rstrip(';') 
                        row = list(map(float, cleaned_line.split()))
                        ras_affine.append(row)

            if len(ras_affine) == 4 and len(ras_affine[0]) == 4:  # Check for a valid 4x4 matrix
                ras_affine = np.array(ras_affine)  # Convert to numpy array
                com_ras = nib.affines.apply_affine(ras_affine, com_voxel)
            else:
                print(f"Expected 4x4 matrix, got shape {np.array(ras_affine).shape}.")
                com_ras = com_t1w 
        else:
            com_ras = com_t1w  

        com_ras_rounded = [round(coord, 2) for coord in com_ras]

        with open(output_file, 'w') as f:
            f.write(f"Center of Mass for {roi_name} (RAS coordinates): {com_ras_rounded}\n")
        print(f"Center of mass coordinates saved to: {output_file}")
    else:
        print(f"No voxels found for ROI {roi_name}.")

config_path = 'config.yaml'
config = load_config(config_path)

fs_dir = config['freesurfer']['fs_directory']
subject_id = config['subject_id']
input_file = config['paths']['T1w']
output_directory = config['freesurfer']['output_directory']
selected_roi_name = config['freesurfer']['roi']['name']
selected_roi_code = int(config['freesurfer']['roi']['mni_code'])
output_file = os.path.join(output_directory, f'{selected_roi_name}_coordinates_RAS.txt')

# Run the FreeSurfer recon all
run_freesurfer_recon(subject_id, input_file, fs_dir)

# Calculate the center of mass and save coordinates in T1w and RAS space
calculate_center_of_mass(selected_roi_name, subject_id, output_file, fs_dir, selected_roi_code)

# Save the ROI mask in T1w space
save_roi_mask(subject_id, fs_dir, selected_roi_code, output_directory, selected_roi_name)
