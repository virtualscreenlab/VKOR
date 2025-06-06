# Chemomodeling analysis of C.sativa compounds as VKOR inhibitors

## Overview
This repository contains code for analyzing and comparing molecular properties of various compounds, including Quercetin, Gallic Acid, Chlorogenic Acid, Catechol, Caffeic Acid, Coumaric Acid, Rosmarinic Acid, and Warfarin, with a control compound 4-Hydroxycoumarin. The analysis includes molecular weight, LogP, hydrogen bond donors, hydrogen bond acceptors, aromatic rings, and free energy (ΔG), with Tanimoto similarity calculations and Lipinski Rule-of-5 compliance checks.

## Features
- Normalization of molecular descriptors.
- Calculation of Tanimoto similarity scores relative to Warfarin and 4-Hydroxycoumarin.
- Visualization of molecular structures and similarity scores.
- Drug-likeness assessment using Lipinski's Rule-of-5.

## Requirements
- Python 3.x
- NumPy
- pandas
- RDKit
- matplotlib

## Installation
1. Clone the repository:
   ```
   git clone https://github.com/virtualscreenlab/VKOR.git
   ```
2. Install the required packages:
   ```
   pip install numpy pandas rdkit matplotlib
   ```

## Usage
Run the script `Td.py` or `Ts.py` to perform the analysis and generate visualizations:
```
python Td.py
```
or
```
python Ts.py
```

## Output
- `vkor_compounds.png`: Grid image of compound structures.
- `similarity_to_4_hydroxycoumarin.png`: Bar plot of Tanimoto similarity to 4-Hydroxycoumarin.
- `similarity_to_quercetin.png`: Bar plot of Tanimoto similarity to Quercetin.

## Contributing
Contributions are welcome! Please fork the repository and submit pull requests.

## License
[Add license information here if applicable.]

## Contact
For contact purposes, you may reach out to Prof. Sergey Shityakov via email at shityakoff@hotmail.com.
