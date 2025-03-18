import os
import pandas as pd
import numpy as np
import subprocess
from Bio import PDB
import pymol
import ast
from Bio.PDB import PDBParser
import seaborn as sns

def process_dataframe(csv_file_path: os.PathLike, output_folder: os.PathLike, GPCRdb_SERVER: str, threshold_clustering: float = None, residues_ortho: list = None):
    """
    Reads a CSV file containing PDB IDs, finds corresponding PDB files in the same folder, 
    retrieves generic residue numbers, and writes an updated CSV.
    
    Parameters:
        csv_file (str): Path to the input CSV file.

    
    Returns: 
        Annotated csv file
    """

    # Read CSV file
    try:
        csv_file = pd.read_csv(csv_file_path)
        print(csv_file.head())  # Display first few rows
    except Exception as e:
        print(f"Error: {e}")


    # Get the folder of the CSV file
    csv_dir = os.path.dirname(os.path.abspath(csv_file_path))
   
    if "PDB" not in csv_file.columns:
        print("Error: CSV file must contain a 'PDB' column with PDB filenames.")
        return

    generic_numbers = []
    generic_numbers_ortho = []
    opm_reference_files = []
    pdb_out = []

    for index, row in csv_file.iterrows():
        pdb_id = row["PDB"]
        chain = row["chain"]
        residues = row["residues"]
        uniprot_id= row["uniprot_comb"].split("_")[0]
        state= row["state"]

        pdb_path = os.path.join(csv_dir, f"{pdb_id}.pdb") 
        out_path = os.path.join(csv_dir, "{}/{}_gn.pdb".format(output_folder, pdb_id)) # Look for PDB file in the same directory
        
        if os.path.exists(pdb_path):
            gen_residues, gen_residues_ortho, opm_file = get_GPCRdb_generic_numbers(pdb_path, out_path, chain, pdb_id, residues, uniprot_id ,state, residues_ortho  =None)
            generic_numbers.append(gen_residues)
            generic_numbers_ortho.append(gen_residues_ortho)
            opm_reference_files.append(opm_file)
            pdb_out.append(out_path)
        else:
            print(f"Warning: {pdb_path} not found. Skipping.")
            generic_numbers.append(None)
            generic_numbers_ortho.append(None)
            opm_reference_files.append(None)
            pdb_out.append(None)

    csv_file["residues_gn"] = generic_numbers
    csv_file["residues_gn_ortho"] = generic_numbers_ortho
    csv_file["opm_file"] =  opm_reference_files
    csv_file["pdb_out"] =  pdb_out
    result_opm = get_opm_location(csv_file)
    result_opm_clustered = clustering_result(result_opm, output_folder, threshold_clustering)
    
def get_GPCRdb_generic_numbers(in_file: os.PathLike, out_file: os.PathLike, chain: str, pdb_id: str, residues: list, uniprot_id: str, state:str, verbose: bool = False, residues_ortho: list = None) -> os.PathLike:
    """
    Use GPCRdb post services to assign GPCRdb numbers to a PDB file returning new file name if successful.
    :param in_file: Path to .pdb file to be sent to GPCRdb server
    :param out_file: Path to save output (e.g., <PDBCODE>_GN.pdb)
    :param verbose: Print any errors
    :return: Union[out_file, None]
    """
    
    if not os.path.exists(out_file):
        
        # Purge b-column values and initiate at 150 (randomly chosen)
        update_b_column(in_file, out_file, 150.0, verbose=True)

        # Assign generic residue numbers
        proc = subprocess.run([
            "curl", "-X", "POST", "-F",
            f"pdb_file=@{out_file}",
            os.path.join(GPCRdb_SERVER, "services/structure/assign_generic_numbers")], capture_output=True)

        stdout = proc.stdout.decode()
            
        if stdout == '':
            errmsg = proc.stderr.decode().split("\n")[3]
            if verbose:
                print(f"Error: {errmsg}")
        else:
            # Save generic residue number file
            with open(out_file, 'wt') as f:
                f.write(stdout)
    else:
        print(f"{out_file} already exists. Skipping processing.")

    

    #Find the corresponding file in OPM 
    current_directory = os.getcwd()
    opm_database_directory = os.path.join(current_directory, '..', 'opm_database')
    # Resolve the full path
    search_directory= os.path.abspath(opm_database_directory)
    opm_file = find_opm_reference(search_directory, in_file, uniprot_id, state)
        
    #Align and save the output
    out_file_aligned = align_to_opm_reference(current_directory, out_file, opm_file, chain)
    #Extract the generic residue numbers
    b_factors = extract_b_factors(out_file, residues, chain)
    b_factors_ortho = extract_orthosteric_site(out_file,chain, residues_ortho = None)
       
    return  b_factors, b_factors_ortho, opm_file

def update_b_column(in_file: os.PathLike, out_file: os.PathLike, new_b_value: float = 150.0, verbose: bool = False) -> os.PathLike:
    """
    This function updates the B-factor column in a PDB file to a new value (randomly initated to 150).

    Args:
    - in_file: Path to the input PDB file.
    - out_file: Path to the output PDB file with updated B-factor values.
    - new_b_value: The new value to set in the B-column.
    - verbose: If True, prints additional information.

    Returns:
    - Path to the output file if successful, None if an error occurs.
    """
    try:
        # Open the input PDB file for reading
        with open(in_file, 'r') as pdb_file:
            lines = pdb_file.readlines()

        # Open the output PDB file for writing
        with open(out_file, 'w') as out_pdb_file:
            for line in lines:
                # Check if the line contains atomic information (starts with "ATOM" or "HETATM")
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # Replace the B-factor (columns 61-66) with the new value
                    updated_line = line[:60] + f"{new_b_value:6.2f}" + line[66:]
                    out_pdb_file.write(updated_line)
                else:
                    # For non-atomic lines, just write them as they are
                    out_pdb_file.write(line)

        if verbose:
            print(f"B-column values purged and saved to {out_file}")
        
        return out_file

    except Exception as e:
        if verbose:
            print(f"Unexpected error: {e}")
        return None

def extract_b_factors(pdb_file: os.PathLike, residues: list, chain: str) -> dict:
    """
    Extracts the B-factor values for C-alpha atoms of specified residues in a given chain.
    
    :param pdb_file: Path to the PDB file (after GPCRdb assignment)
    :param residues: List of residue numbers to extract B-factors for
    :param chain: Chain name to filter the residues
    :return: Dictionary with residue number as the key and B-factor value as the value
    """
    b_factors = {}
    
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract chain, residue number, and atom name
                chain_id = line[21]
                residue_number = int(line[22:26].strip())
                
                atom_name = line[12:16].strip()
                
                
                # Ensure it's a C-alpha atom (CA) and matches the chain and residue number
                if atom_name == "CA" and chain_id == chain and int(residue_number) in list(map(int, residues.strip("[]").split(", "))):
                    b_factor = float(line[60:66].strip())
                    b_factors[residue_number] = b_factor
        
    return b_factors

def extract_orthosteric_site(pdb_file: os.PathLike, chain: str, residues_ortho: list) -> dict:
    """
    Extracts the residue values for C-alpha atoms of specified B-factor in a given chain.
    
    Parameters:
        pdb_file: Path to the PDB file (after GPCRdb assignment)
        residues (list): List of residue numbers to extract B-factors for
        chain: Chain name to filter the residues
    Returns: 
        Dictionary with residue number as the key and B-factor value as the value
    """
    b_factors_ortho = {}
    orthosteric_residues_major = [3.32, 3.33, 3.36, 45.52, 5.43, 6.48, 6.51, 6.52, 6.55, 7.38, 7.42]
    orthosteric_residues_minor = [1.39, 2.53, 2.56, 2.57, 2.60, 2.64, 3.32, 3.35, 7.35, 7.38, 7.39, 7.42, 7.46]
    orthosteric_residues = list(set(orthosteric_residues_major + orthosteric_residues_minor))
    
    with open(pdb_file, 'r') as file:
        for line in file:
            
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract chain, residue number, and atom name
                chain_id = line[21]
                b_factor= float(line[60:66].strip())
                atom_name = line[12:16].strip()
                
                # Ensure it's a C-alpha atom (CA) and matches the chain and residue number
                if atom_name == "CA" and chain_id == chain and float(b_factor) in orthosteric_residues:
                    residue_number = int(line[22:26].strip())
                    b_factors_ortho[residue_number] = b_factor
        
    return b_factors_ortho

def find_opm_reference(search_directory: os.PathLike, pdb_path: str, uniprot_id: str, state:str):
    """
    The matching structure is searched for in the folder opm_database. If the exact PDB is not found, the first PDB matching with the UNIPROT name and state is chosen. 
    Other wise an error is raised and the row is skipped. 
    Parameters: 
        search_directory (path): Path to folder where OPM reference files are stored. 
        pdb_path (path): Path to query file
        uniprot_id (str): Protein Family
        state (str): Protein State
    Returns:
        Representative OPM reference file either identical or matching the protein family and statess
    
    """
    # Check if the directory exists
    if not os.path.exists(search_directory):
        raise FileNotFoundError(f"The OPM database directory was not found at {search_directory}")
    
    # List all files in the directory
    pdb_substring = pdb_path.split("/")[-1].replace(".pdb", "")
    files_in_directory = os.listdir(search_directory)

    # Find files that contain the substring in their filename
    matching_files = [f for f in files_in_directory if pdb_substring in f]
    
    # If no files match, raise an error
    if not matching_files:
        matching_files = [f for f in files_in_directory if uniprot_id in f and state in f]
        if not matching_files:
            raise FileNotFoundError(f"No OPM reference file found containing the substring '{pdb_substring}' in the directory.")
    
    # Return the first matching file (you can modify this to return multiple if needed)
    opm_file = os.path.join(search_directory, matching_files[0])
    return opm_file

def align_to_opm_reference(current_directory: os.PathLike, out_file: os.PathLike, opm_file: os.PathLike, chain: str):
    """
    Aligns a given protein structure (query) to a reference structure and saves the aligned query.
    
    Args:
    - current_directory: Directory where the structure files are located.
    - out_file: The path where the aligned query structure will be saved.
    - opm_file: The path to the reference structure file (structure to align to).
    - chain: The chain ID of the query structure to retain (removes other chains).
    
    Returns:
    - The path to the aligned query structure (out_file).
    """
   
    # Step 1: Start PyMOL in quiet mode (no GUI)
    pymol.finish_launching(['pymol', '-cq'])
    pymol.cmd.reinitialize()
    pymol.cmd.load(opm_file, 'reference')
    pymol.cmd.load(out_file, 'query')
    
    # Delete chains that are not part of the specified chain in the query structure
    delete_not_necessary_chains = "query and not chain {}".format(chain)
    pymol.cmd.remove(delete_not_necessary_chains)
    
    # Step 2: Superimpose structures
    
    pymol.cmd.super('query', 'reference')
    
    # Step 3: Save the aligned query structure to the specified file path
    pymol.cmd.save(out_file, 'query')
    
    # Step 4: Delete loaded structures from PyMOL session
    pymol.cmd.delete('reference')
    pymol.cmd.delete('query')
    
    print("Completed superposition.")
    
    # Return the path to the aligned structure file
    return out_file

def custom_function(row: pd.Series):
    """
    Computes the relation of the Distance to the OPM membrane boundaries 

    Parameters:
    row (pd.Series): Pandas row containing information about the OPM boundary location

    Returns:
    pd.Series: Calculates the mean for each pocket 
    """
    if pd.isna(row['Mem_out']) or pd.isna(row['Mem_in']) or row['Mem_in'] == None:
        return [np.nan] * len(row['Z_CA'])  # Return a list of NaNs if the function was not applied

    #Location in membrane
    results = []
    for z_ca in row['Z_CA']:
        if z_ca > row['Mem_out']:
            results.append(11)
        elif z_ca < row['Mem_in']:
            results.append(-11)
        else:
            results.append(round((z_ca / row['Mem_out']) * 10, 0))
    
    if results:
        mean_results = np.mean(results)
    else:
        mean_results = np.nan
    
    #Distance to center
    results_distance = []
    for dist in row['Distance_CA']:
        results_distance.append(dist)
    if results_distance:
        mean_results_distance = np.mean(results_distance)
    else:
        mean_results_distance = np.nan
    

    #Distance to center with respect to orthosteric site
    results_distance_shift = []
    for dist in row['Distance_CA_shift']:
        results_distance_shift.append(dist)
    if results_distance_shift:
        mean_results_distance_shift = np.mean(results_distance_shift)
    else:
        mean_results_distance_shift = np.nan

    
    return results, mean_results, mean_results_distance, mean_results_distance_shift

def distance_to_center(x: float, y: float):
    #Calculate the distance from a point (x, y) to the center of the circle (0, 0).
    """
    Parameters:
    x (float): x-coordinate of the point
    y (float): y-coordinate of the point
    
    Returns:
    float: Distance from the point to the center of the circle
    """
    return np.sqrt(x**2 + y**2)

def distance_to_center_shifted(x: float, y: float, center_x: float, center_y: float):
    """
    Calculate the distance from a point (x, y) to the center of the circle (center_x, center_y).
    Parameters:
    x (float): x-coordinate of the point
    y (float): y-coordinate of the point
    center_x (float): x-coordinate of the center of the circle
    center_y (float): y-coordinate of the center of the circle
    
    Returns:
    float: Distance from the point to the center of the circle
    """
    return np.sqrt((x - center_x)**2 + (y - center_y)**2)

def str_to_list(s: str):
    """"
    Conver string values to a list
    Parameter:
    s (str): String value
    Returns: 
    List of values
    """
    return ast.literal_eval(s)

def get_opm_location(df: pd.Series):
    """
    Calculate the relation to the OPM defined membrane
    Parameters:
    df (pd.Series): Pandas Dataframe
    
    Returns:
    df (pd.Series): Pandas Dataframe modifed with location relative to OPM reference boundaries
    """
   
    current_directory = os.getcwd()
    # Create a new columns
    df["Z_CA"] = None
    df["Distance_CA"] = None
    df["Distance_CA_shift"] = None
    df["Mem_out"] = None
    df["Mem_in"] = None
    # Apply the function to the DataFrame column
    df['residues'] = df['residues'].apply(str_to_list)
    
    print("Let's start assigning the opm reference: ")
    for index, row in df.iterrows():
        #Find the aligned structure 
       
        keyword_1 = row['PDB']
        keyword_2 = row['uniprot_comb'].split("_")[0]
        keyword_3 = row['chain']
        list_residues = row['residues']
        
        print("This is the total number of residues", len(list_residues))
        list_residues_orthosteric = list(row['residues_gn_ortho'].keys())
        list_caz = []
        distance_ca = []
        distance_ca_shifted = []
        
        print("This file is processed:", keyword_1 , keyword_2, keyword_3)

        # Walk through the directory and search for files
        search_directory_aligned = row["pdb_out"]
       
        #Load the aligned PDB of interest
        file_path = os.path.join(search_directory_aligned)
      
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(keyword_1, file_path)
        model = structure[0]
        x_total = 0
        y_total = 0
        z_total = 0
        count = 0

        # Step1: Loop through the orthosteric residues to detect center
        for residue_number in list_residues_orthosteric:
            for model in structure.get_list():
                for chain in model.get_list():
                        if chain.id == keyword_3: 
                                for residue in chain:
                                        if (residue.get_id()[1] == int(residue_number)) and residue.has_id("CA"):
                                            ca = residue["CA"]
                                            #Calcualte the depth of the residues
                                            # Calculate the coordinates of the residues
                                            x_total += ca.coord[0]
                                            y_total += ca.coord[1]
                                            z_total += ca.coord[2]
                                            count += 1
        # Calculate the center if any residues were found
        if count > 0:
            center_x = x_total / count
            center_y = y_total / count
            center_z = z_total / count
            center = (center_x, center_y, center_z)
        else:
            center = None
                    

        # Step3: Loop through each row in the DataFrame
        for residue_number in list_residues:
            for model in structure.get_list():
                for chain in model.get_list():
                    if chain.id == keyword_3: 
                        for residue in chain:
                            if (residue.get_id()[1] == int(residue_number)) and residue.has_id("CA"):
                                ca = residue["CA"]
                                #Calcualte the depth of the residues
                                x_cooridnate = ca.coord[0]
                                y_coordinate = ca.coord[1]
                                #The opm reference dataframe aligns the structure so that 0, 0, 0 is the center of the membrnae
                                #Therefore distance to the middle axis can be calculated by just taking x and y axis into account
                                distance = distance_to_center(x_cooridnate, y_coordinate)
                                distance_shifted = distance_to_center_shifted(x_cooridnate, y_coordinate, center_x, center_y)
                                distance_ca.append(distance)
                                distance_ca_shifted.append(distance_shifted )
                                #Exctract the location in respect to the membrane
                                z_coordinate = ca.coord[2]
                                list_caz.append(z_coordinate)

                        
        #Add information to the dataframe                   
        df.at[index, "Distance_CA"] = distance_ca
        df.at[index, "Distance_CA_shift"] = distance_ca_shifted
        df.at[index, "Z_CA"] = list_caz
        
        #Step 4: Find OPM boundaries
        pdb_reference= row["opm_file"]
        # Open the text file (replace 'your_file.txt' with your actual file path)
        search_string = "DUM"
        values_all = []
        with open(pdb_reference, 'r') as file:
            for line in file:
                if search_string in line:
                    values = line.split()
                    if len(values) == 10:
                        values_all.append(float(values[6]))
                    elif len(values) == 11:
                        #OPM e.g. 7LD3 has a different format than 5UEN
                        values_all.append(float(values[7]))
                    elif len(values) == 12:
                        #OPM e.g. 7LD3 has a different format than 5UEN
                        values_all.append(float(values[8]))
                    else:
                        values_all.append(float(values[-1]))
                            
        unique_values = np.unique(values_all)
        df.at[index, "Mem_out"] = unique_values[1]
        df.at[index, "Mem_in"] = unique_values[0]  
                    
    # Filter out rows where 'Z_CA' is empty
    df_filtered = df[df['Z_CA'].notna()]

    # Apply custom function to the filtered DataFrame
    df_filtered['Mem_loc'], df_filtered['Mem_loc_avg'], df_filtered['Mem_distance_avg'], df_filtered['Mem_distance_avg_shift']= zip(*df_filtered.apply(custom_function, axis=1))
    
    # Merge the new columns back into the original DataFrame
    df = df.merge(df_filtered[['Mem_loc', 'Mem_loc_avg','Mem_distance_avg', 'Mem_distance_avg_shift']], left_index=True, right_index=True, how='left')
    
    print("Z coordinates calculated and stored in CSV files.")
    return df   
    
def clustering_result(df: pd.Series, output_path: os.PathLike, threshold_clustering: float =None):
    """
    Cluster pockets based on generic residue number from GPCRdb
    Parameters:
    df (pd.Series): Pandas Dataframe
    threshold_clustering (float): Number used as clustering threshold (between 0 and 1 )
    
    Returns:
    df (pd.Series): Pandas Dataframe modifed with the binding site annotations
    """
    cluster_csv = []
    # Replace NaN or empty strings with 'Other'
    df['state'] = df['state'].fillna('Other')  # Replace NaN with 'Other'
    df['state'] = df['state'].replace('', 'Other')  # Replace empty strings with 'Other'
    # Map the class values to the corresponding letters
    

    
    # Create the new column 'apple' with the mapped values joined with the value in 'state' column
    df['Class_state'] = df['class'] + '-' + df['state']
    df['Class_state'] = df['Class_state'].replace('A-Intermediate', 'A-Inactive_Intermediate')
    df['Class_state'] = df['Class_state'].replace('A-Inactive', 'A-Inactive_Intermediate')
    df['Class_state'] = df['Class_state'].replace('B1-Intermediate', 'B1-Inactive_Intermediate')
    df['Class_state'] = df['Class_state'].replace('B1-Inactive', 'B1-Inactive_Intermediate')
    df['Class_state'] = df['Class_state'].replace('B2-Intermediate', 'B2-Inactive_Intermediate')
    df['Class_state'] = df['Class_state'].replace('B2-Inactive', 'B2-Inactive_Intermediate')
    df['Class_state'] = df['Class_state'].replace('C-Intermediate', 'C-Inactive_Intermediate')
    df['Class_state'] = df['Class_state'].replace('C-Inactive', 'C-Inactive_Intermediate')
    df['Class_state'] = df['Class_state'].replace('F-Intermediate', 'F-Inactive_Intermediate')
    df['Class_state'] = df['Class_state'].replace('F-Inactive', 'F-Inactive_Intermediate')
    df['Class_state'] = df['Class_state'].replace('T-Intermediate', 'T-Inactive_Intermediate')
    df['Class_state'] = df['Class_state'].replace('T-Inactive', 'T-Inactive_Intermediate')
    df['Class_state'] = df['Class_state'].replace('D1-Intermediate', 'D1-Inactive_Intermediate')
    df['Class_state'] = df['Class_state'].replace('D1-Inactive', 'D1-Inactive_Intermediate')
    
   
    folder_names = df['Class_state'].unique().tolist()
    df['gn_residues'] = df['residues_gn'].apply(lambda x: list(x.values()) if isinstance(x, dict) else [])

    
     # loop over the folder names
    for folder in folder_names:
        print("We are processing:", folder)
        # create the folder if it doesn't exist
        current_dir = output_path+"/"+folder
        dir_name = current_dir+"/"
        folder_name = os.path.join(current_dir, folder)
        folder_name = folder_name+"/"
    
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        
        # filter the dataframe based on the folder name
        filtered_df = df[df['Class_state'] == folder]
        sites = filtered_df.reset_index(drop=True)

        #################################
        #Plot interactions_bw against PDB
        #################################
        interactions_bw_extend = []
        interactions_bw = []
        sites["gn_residues_filtered"]  = None
        for index, row in sites.iterrows():
            #Get all generic residues that are present in Class but separated by list
            sites.at[index, "gn_residues_filtered"] = [num for num in sites.at[index, "gn_residues"] if isinstance(num, (int, float)) and num <= 100]
            
            interactions_bw.append(sites.at[index, "gn_residues_filtered"] )
            #Get the all generic residues that are present in aClass
            interactions_bw_extend.extend(sites.at[index, "gn_residues_filtered"]  )
        
        
        # Order the extended labels
        label = interactions_bw_extend
        x_labels_unq = np.unique(label).tolist()
        
        x_labels_ordered = []
        lst = ['N', "TM1", '1.', 'ICL1', "TM2", "2.", "ECL1", "23.", "TM3", "3.", "ICL2", "34.", "TM4",
            "4.", "ECL2", "45.", "TM5", "5.", "ICL3", "56.", "6.", "ECL3", "78.", "TM7", "7.", "H8",
            "8.", "C-", "ICL4", "-"]

        # Only use x_labels that are indeed in the list
        for value in lst:
            items = [k for k in x_labels_unq if str(value) in str(k)]
            x_labels_ordered.extend(items)

        x_labels = list(set(x_labels_ordered)) 
    
        # Make a matrix with 0 if residue is not part of pocket or 1 if it is part of the pocket
        n_conform = len(interactions_bw)
        Figure_data = np.zeros((n_conform, len (x_labels)))
        for i in range(0,len(interactions_bw)):
            for inters in interactions_bw[i]:
                idx = x_labels.index(inters) 
                Figure_data[i,idx] = 1

        Figure_data_pd = pd.DataFrame(Figure_data, columns =x_labels)
        
        ######################
        #Clustering
        ######################
        #Jaccard
        res  = sns.clustermap(Figure_data_pd, method = "average", metric = "jaccard",    cbar_pos = None)
        
        #clustered dataframe, res2.data2d, use index 2d to order the original dataframe according 
        data = res.data2d

        #fcluster to identify clusters based on threshold
        if threshold_clustering != None: 
            cluster_array = scipy.cluster.hierarchy.fcluster(res.dendrogram_row.linkage, t=0.75, criterion='distance' )
            data["cluster"] = cluster_array
            sites["cluster"] = cluster_array

        else:#Treating each ligand as a separate site
            number_list = []
            for i in range(1, len(Figure_data_pd) + 1):
                number_list.append(i)
            data["cluster"] = number_list
            sites["cluster"] = number_list
        
        #########################
        #Annotation
        ########################
        results = sites

        unique_value = results.cluster.unique()
        out_arr = np.argsort(unique_value)
        unique_value = unique_value[out_arr ]

        #Iterate over cluster sites
        for cluster in unique_value: 
            df_sites = results.loc[results['cluster'] == cluster]
            filtered_rows = results.loc[results['cluster'] == cluster].index
            
            #Determine sublocation based on OPM 
            location_OPM_avg = df_sites['Mem_loc_avg'].dropna().mean()
            location_IH_avg_shift =df_sites['Mem_distance_avg_shift'].dropna().mean()
            location_IH_avg =df_sites['Mem_distance_avg'].dropna().mean()
            
            #####
            #logic of the thresholds: from pocket_annotation.py (line 129: round((row['Z_CA'] / row['Mem_out']) * 10, 0)) therefore if row['Z_CA'] == row['Mem_out'] = 1*10 =10
            #therefore we define mid is 0 +/- 5
            
            if location_OPM_avg == 11:
                location_suffix = ""
                location_suffix2 = "EC"
            elif location_OPM_avg > 5: 
                location_suffix = "ext"
                location_suffix2 = "IH"
            elif 5 >= location_OPM_avg >= -5:
                location_suffix = "mid"
                location_suffix2 = "IH"
            elif -5> location_OPM_avg > -11:
                location_suffix = "int"
                location_suffix2 = "IC"
            else:
                location_suffix = "int"
                location_suffix2 = "IC"
               
            num = len(df_sites)
            #Determine number which corresponds to 2/3 of the sites
            over_twothird = int(num/3*2)
    
            sites = []
                
            for row, value in df_sites.gn_residues_filtered.items():
                sites.extend(value)
                
            #Determine unique values 
            values, counts = np.unique(sites, return_counts=True)
            unique_val = dict(zip(values, counts))
            unique_val_filtered = [] 
                
            #Retain all GN numbers that occur more than 2/3 of all sites
            unique_val_filtered = {k: v for k, v in unique_val.items() if v >= over_twothird}
            list_val = []
            for k, v in  unique_val_filtered.items():

                split_string = str(k).split("x", 1)
                k = split_string[0]
                value = k.replace('.' , 'x')
                list_val.append(value)
                
            print(cluster, *list_val, sep=", ")   
                
            #Get tm values
            tm_val = []
            for entry in list_val: 
                value = entry[0]
                tm_val.append(value)
            #Filter for unique values
            values, counts = np.unique(tm_val, return_counts=True)
            
            
            #Filter for numeric values
            output = ''
            for character in values:
                if character.isnumeric():
                    output += character
                else:
                    continue
            
            print("This is cluster", cluster,":",  output, location_suffix, location_suffix2)
            #Compare to Kolb annotation
            dict_Kolb = {"EH-TM456":456 ,"EH-TM2345":2345, "EH-TM12": 12,  "EH-TM35":35, "EH-TM356":356,"EH-TM56": 56,"EH-TM67": 67, "EH-TM17": 17,"EH-TM178":178, "EH-TM18": 18, "EH-TM23": 23, "EH-TM34": 34, "EH-TM24": 24, "EH-TM45": 45, "EH-TM345": 345, "EH-TM78": 78, "EH-TM234": 234, "EH-TM124": 124,"EH-TM167": 167, "EH-TM123": 123,  "EH-TM2678": 2678, "EH-TM567": 567, "EH-TM5678": 5678, "EH-TM678": 678}
            matches = []

            for key,value in dict_Kolb.items():
                # list data '4','1','5' are in string and dictonery value 4,1,5 are in integer form
                # hence you need to compare the same data type 
                if value != '':
                    try: 
                        if value == int(float(output)):
                            matches.append(key)
                        else:
                            matches.append("no match")
                    except:
                        pass
            #print("-----------------------------------------")
            matches = np.unique(matches)

            #If there are no matches site will be declared as external site
            if len(matches) == 0:
                if folder[0] == "A" or folder[0] == "C" or folder[0] =="F" or folder[0] =="T":
                    matches_nam = "NoTMs"
                    matches = folder[0]+"-"+matches_nam
                    value_list = [matches]* len(filtered_rows)
                    results.loc[filtered_rows, 'Site'] = value_list
                else:
                    matches_nam = "NoTMS"
                    matches = folder[0:2]+"-"+matches_nam
                    value_list = [matches]* len(filtered_rows)
                    results.loc[filtered_rows, 'Site'] = value_list
                
            #If there are matches the class will be added to the name
            elif matches[0] != "no match" and location_suffix != None:
                if folder[0] == "A" or folder[0] == "C" or folder[0] =="F" or folder[0] =="T":
                    #Determined based on evaluation 8th August 2024
                    if folder[0] =="F" and location_IH_avg_shift < 8.4:
                        matches = folder[0]+"-"+location_suffix2+"-TM"+output
                        value_list = [matches]* len(filtered_rows)
                        results.loc[filtered_rows, 'Site'] = value_list
                    elif folder[0] =="A" and matches[0] =="EH-TM456" and location_IH_avg_shift < 14.0:
                        matches = folder[0]+"-"+location_suffix2+"-TM"+output
                        value_list = [matches]* len(filtered_rows)
                        results.loc[filtered_rows, 'Site'] = value_list
                    else:
                        matches = folder[0]+"-"+matches[0]+"_"+location_suffix
                        value_list = [matches]* len(filtered_rows)
                        results.loc[filtered_rows, 'Site'] = value_list
                else:
                    matches = folder[0:2]+"-"+matches[0]+"_"+location_suffix
                    value_list = [matches]* len(filtered_rows)
                    results.loc[filtered_rows, 'Site'] = value_list
                
            #If there are GN, but no matches are found the user wil have to declare if it is IH/IC or EC
            else:
                if folder[0] == "A" or folder[0] == "C" or folder[0] =="F"or folder[0] =="T":
                    if folder[0] =="A" and  location_OPM_avg > 6.45:
                        matches = folder[0]+"-"+location_suffix2+"-ECV-TM"+output
                        value_list = [matches]* len(filtered_rows)
                        results.loc[filtered_rows, 'Site'] = value_list
                    else:
                        matches = folder[0]+"-"+location_suffix2+"-TM"+output
                        value_list = [matches]* len(filtered_rows)
                        results.loc[filtered_rows, 'Site'] = value_list
                else:
                    matches = folder[0:2]+"-"+location_suffix2+"-TM"+output
                    value_list = [matches]* len(filtered_rows)
                    results.loc[filtered_rows, 'Site'] = value_list

            output_filename = folder_name + 'FINAL_Clustered_results.csv'
            results.to_csv(output_filename)
            
        cluster_csv.append(output_filename)
        
    # Loop over the file paths and concatenate the CSV files
    concatenated_df = pd.DataFrame()
    for file_path in cluster_csv:
        df = pd.read_csv(file_path)
        concatenated_df = pd.concat([concatenated_df, df])

    # Write the concatenated DataFrame to a new CSV file
    # Delete all columns with the name "Unnamed"
    concatenated_df = concatenated_df.loc[:, ~concatenated_df.columns.str.contains('^Unnamed')]
    # Drops temporary columns required for pocannos
    concatenated_df.drop(columns=[
    "opm_file", "pdb_out", "Z_CA", "Distance_CA", "Distance_CA_shift",
    "Mem_out", "Mem_in", "Mem_loc", "Mem_loc_avg", "Mem_distance_avg",
    "Mem_distance_avg_shift", "Class_state", "gn_residues"
    ], inplace=True, errors='ignore')
    concatenated_df.to_csv("{}/Annotated_file.csv".format(output_path), index=False)
    print(" Congrats, you did it! Please check in your output folder the Annotated_file.csv column 'Site'." )
    return concatenated_df
    

