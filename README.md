# GPCR POCket ANNotationS (pocannos)

---

This repository contains Python code to annotate GPCR binding sites following the procedure published in:

**"Comparative Study of Allosteric GPCR Binding Sites and Their Ligandability Potential"**  
[Publication Link](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00819)

The GPCR pocannos includes:
- Receptor class
- Transmembrane contacts
- For membrane-facing sites: membrane location

![GPCR Pocannos](pocannos.png)

---

## Installation and Use

### Clone the Repository
```bash
git clone https://github.com/sonjapeter/pocannos.git
cd pocannos
```

### Create and Activate a Conda Environment
```bash
conda create -n pocannos_env python=3.9 -y
conda activate pocannos_env
```

### Install Required Packages
```bash
conda install -y -c conda-forge os pandas numpy subprocess biopython pymol-open-source ast seaborn
```

---

## External Dependencies
This code requires access to the **GPCRdb generic residue numbering API** ([GPCRdb](https://gpcrdb.org/)) and **OPM reference files** included in the `opm_database` folder.

---

## Scripts

### `main.py`
This script allows a user to pass a CSV file containing the residues of a pocket of interest and returns the GPCR pocannos annotations.

```bash
python ./main.py --help
```

#### Add GPCR pocannos to a CSV file

```bash
Usage: python main.py [arguments]
```

#### Arguments:
```
csv_file            Path to input CSV file
output_folder       Path to output folder
```

#### Optional Arguments:
```
-h, --help          Show this help message and exit
threshold_clustering
                    A float value between 0 and 1 can be used to cluster binding sites
                    (default: every site treated as separate cluster)
residue_ortho
                    A list of residues that define the orthosteric site
                    (default: orthosteric major and minor pocket residues:
                    1.39, 2.53, 2.56, 2.57, 2.60, 2.64, 3.32, 3.33,
                    3.35, 3.36, 45.52, 5.43, 6.48, 6.51, 6.52, 6.55,
                    7.38, 7.39, 7.42, 7.46)
```

---

## License
This project is licensed under the MIT License.

## Contact
For any inquiries, please contact sonja.peter3@hispeed.ch.

