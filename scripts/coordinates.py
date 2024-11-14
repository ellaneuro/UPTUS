import nibabel as nib
import numpy as np
import yaml
import os
import re

# Function to load config file
def load_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

# Function to convert RAS coordinates to voxel coordinates
def ras_to_voxel(ras_coords, affine_matrix):
    # Convert RAS to homogeneous coordinates
    ras_coords_homogeneous = np.append(ras_coords, 1)
    
    # Inverse affine matrix to convert from RAS to voxel coordinates
    inv_affine_matrix = np.linalg.inv(affine_matrix)
    
    # Get the voxel coordinates (in homogeneous coordinates)
    voxel_coords_homogeneous = inv_affine_matrix.dot(ras_coords_homogeneous)
    
    # Extract and round the voxel coordinates
    voxel_coords = np.round(voxel_coords_homogeneous[:3]).astype(int)
    
    return voxel_coords

# Function to load target coordinates from file
def load_target_coords(file_path):
    if not os.path.exists(file_path):
        print(f"Error: Target coordinates file {file_path} does not exist.")
        return None
    
    try:
        with open(file_path, 'r') as f:
            # Read the line and use regular expressions to extract numbers
            line = f.readline().strip()
            # Regular expression to extract all numbers (including decimal points)
            coordinates = re.findall(r"[-+]?\d*\.\d+|\d+", line)
            # Convert the extracted values to floats and return as a list of coordinates
            target_coords = list(map(float, coordinates))
        
        return target_coords
    except Exception as e:
        print(f"Error reading target coordinates: {e}")
        return None

# Function to prompt user for entry point RAS coordinates
def get_entry_point_coords():
    print("Please enter the entry point coordinates (RAS) from Connectome Workbench (x y z):")
    entry_point_str = input()
    # Convert to float instead of int
    entry_point_coords = list(map(float, entry_point_str.strip().split()))
    return entry_point_coords

# Function to read the existing config file and update it
def update_config_with_voxel_coords(config_file, entry_point_voxel, target_voxel):
    try:
        # Make voxel coordinates positive by taking the absolute value
        entry_point_voxel_abs = np.abs(entry_point_voxel).tolist()
        target_voxel_abs = np.abs(target_voxel).tolist()

        # Load the existing config.yaml
        with open(config_file, 'r') as file:
            config_data = yaml.safe_load(file)
        
        # Store the voxel coordinates as positive values [x, y, z] format
        config_data['coordinates'] = {
            'entry_point_voxel': entry_point_voxel_abs,
            'target_voxel': target_voxel_abs
        }
        
        # Write the updated config back to the YAML file
        with open(config_file, 'w') as file:
            yaml.dump(config_data, file, default_flow_style=False)
        
        print(f"Voxel coordinates successfully written to {config_file}")
    
    except Exception as e:
        print(f"Error updating the config file: {e}")

# Main script logic
def main():
    # Load the config file
    config_file = 'config.yaml'  # Path to your config.yaml
    
    # Load configuration data from the config.yaml
    config = load_config(config_file)
    
    # Extract necessary paths from the config file
    output_directory = config['freesurfer']['output_directory']
    selected_roi_name = config['freesurfer']['roi']['name']
    
    # Get T1w file path from config file
    nifti_file = config['paths']['T1w']
    
    # Construct the output file name from ROI name and output directory
    output_file = os.path.join(output_directory, f'{selected_roi_name}_coordinates_RAS.txt')
    
    # Load the NIfTI image (for affine matrix)
    img = nib.load(nifti_file)
    affine_matrix = img.affine
    
    # Load target coordinates from the text file
    target_coords_file = os.path.join(output_directory, f'{selected_roi_name}_coordinates_RAS.txt')  # Pull the target coordinates file from output folder
    target_coords = load_target_coords(target_coords_file)
    if target_coords is None:
        print("Failed to load target coordinates.")
        return
    
    # Convert target coordinates from RAS to voxel
    target_voxel = ras_to_voxel(target_coords, affine_matrix)
    print(f"Target coordinates (RAS): {target_coords}")
    print(f"Converted target voxel coordinates: {target_voxel}")
    
    # Get entry point coordinates from user input
    entry_point_coords = get_entry_point_coords()
    print(f"Entry point coordinates (RAS): {entry_point_coords}")
    
    # Convert entry point coordinates from RAS to voxel
    entry_point_voxel = ras_to_voxel(entry_point_coords, affine_matrix)
    print(f"Converted entry point voxel coordinates: {entry_point_voxel}")
    
    # Update the existing config file with voxel coordinates
    update_config_with_voxel_coords(config_file, entry_point_voxel, target_voxel)
    
    # Save the RAS coordinates to the output file
    with open(output_file, 'w') as f:
        f.write(f"Target RAS coordinates: {target_coords}\n")
        f.write(f"Entry point RAS coordinates: {entry_point_coords}\n")
    print(f"Entry + target coordinated (RAS) updated: {output_file}")

if __name__ == "__main__":
    main()
