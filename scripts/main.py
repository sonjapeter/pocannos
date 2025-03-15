"""
main.py: Main script that executes the data processing workflow.
"""

from functions import get_GPCRdb_generic_numbers, process_dataframe


def main():
    import sys
    import os

    # Check if the filename and output folder are provided
    if len(sys.argv) < 3:
        print("Usage: python script.py <csv_file> <output_folder> [threshold] [residues_ortho]")
        sys.exit(1)

    # Get the file path from command-line arguments
    csv_file_path = sys.argv[1]

    # Get the output folder from command-line arguments
    output_folder = sys.argv[2]

    # Ensure the output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Output saved in folder: {output_folder}")

    # Check if a threshold parameter is provided (optional)
    threshold = None
    if len(sys.argv) >= 4:
        try:
            threshold = float(sys.argv[3])  # Convert threshold to float if provided
            if not (0 <= threshold <= 1):
                raise ValueError("Threshold must be between 0 and 1.")
        except ValueError:
            print("Invalid threshold value. It should be a number.")
            sys.exit(1)

    # Check if orthosteric residues are provided (optional)
    residues_ortho = []
    if len(sys.argv) >= 5:
        residues_ortho = sys.argv[4] # Expect residues to be comma-separated
        print(f"Orthosteric residues provided: {residues_ortho}")
      
    # Process the CSV file with the optional threshold, residues, and output folder
    csv_file_gn = process_dataframe(csv_file_path, output_folder, threshold, residues_ortho)
    
# Ensure script runs only when executed directly
if __name__ == "__main__":
    main()



