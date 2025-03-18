"""
main.py: Main script that executes the pocannos processing workflow.
"""

from functions import process_dataframe


def main():
    import sys
    import os
    def print_help():
        help_text = """
        Usage: python pocannos.py <csv_file> <output_folder> [threshold] [residues_ortho]

        Arguments:
        csv_file
                        Path to file containing one row per binding site
        output_folder
                        Path to the output folder
        
    
        Optional Arguments:
        -h, --help      Show this help message and exit
        threshold_clustering
                        A float value between 0 and 1 can be used to cluster binding sites
                        (default: every site treated as separate cluster)
        residue_ortho
                        A list of residues that define the orthosteric site
                        (default: orthosteric major and minor pocket residues:
                        1.39, 2.53, 2.56, 2.57, 2.60, 2.64, 3.32, 3.33,
                        3.35, 3.36, 45.52, 5.43, 6.48, 6.51, 6.52, 6.55,
                        7.38, 7.39, 7.42, 7.46)
        """
        print(help_text)
        sys.exit(0)

    # Check if help flag is provided
    if "-h" in sys.argv or "--help" in sys.argv:
        print_help()

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
            print("Invalid threshold value. It should be a number between 0 and 1.")
            sys.exit(1)

    # Check if orthosteric residues are provided (optional)
    residues_ortho = []
    if len(sys.argv) >= 5:
        residues_ortho = sys.argv[4]  # Expect residues to be comma-separated
        print(f"Orthosteric residues provided: {residues_ortho}")
    
    #Define the GPCRDB API
    GPCRdb_SERVER = 'https://gpcrdb.org/'
    # Process the CSV file with the optional threshold, residues, and output folder
    csv_file_gn = process_dataframe(csv_file_path, output_folder, GPCRdb_SERVER, threshold, residues_ortho)

    
    
# Ensure script runs only when executed directly
if __name__ == "__main__":
    main()



