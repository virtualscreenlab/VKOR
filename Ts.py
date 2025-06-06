import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, Lipinski
import matplotlib.pyplot as plt

# =============================================
# 1. Define Compounds & Control (4-Hydroxycoumarin)
# =============================================
compounds = {
    "Quercetin": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O",
    "Gallic Acid": "C1=C(C=C(C(=C1O)O)O)C(=O)O",
    "Chlorogenic Acid": "C1C(C(C(CC1(C(=O)O)O)OC(=O)/C=C/C2=CC(=C(C=C2)O)O)O)O",
    "Catechol": "C1=CC=C(C(=C1)O)O",
    "Caffeic Acid": "C1=CC(=C(C=C1C=CC(=O)O)O)O",
    "Coumaric Acid": "C1=CC(=CC=C1C=CC(=O)O)O",
    "Rosmarinic Acid": "C1=CC(=C(C=C1CC(C(=O)O)OC(=O)C=CC2=CC(=C(C=C2)O)O)O)O",
    "Warfarin": "CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O",
    "4-Hydroxycoumarin (Control)": "C1=CC=C2C(=C1)C(=C(C(=O)O2)O)"
}


# =============================================
# 2. Compute Molecular Properties
# =============================================
def analyze_compounds(compounds):
    data = []
    for name, smiles in compounds.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Error: Invalid SMILES for {name}: {smiles}")
            continue
        fp = Chem.RDKFingerprint(mol)
        data.append({
            "Name": name,
            "SMILES": smiles,
            "Molecular Weight": Descriptors.MolWt(mol),
            "LogP": Descriptors.MolLogP(mol),
            "H-Bond Donors": Lipinski.NumHDonors(mol),
            "H-Bond Acceptors": Lipinski.NumHAcceptors(mol),
            "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
            "Aromatic Rings": Lipinski.NumAromaticRings(mol),
        })
    return pd.DataFrame(data)


df = analyze_compounds(compounds)
if df.empty:
    raise ValueError("No valid compounds processed. Check SMILES strings.")
print("\nMolecular Properties:")
print(df[["Name", "Molecular Weight", "LogP", "H-Bond Donors","H-Bond Acceptors", "Aromatic Rings"]].to_string(index=False))

# =============================================
# 3. Visualize Structures
# =============================================
try:
    mols = [Chem.MolFromSmiles(smiles) for smiles in compounds.values() if Chem.MolFromSmiles(smiles)]
    img = Draw.MolsToGridImage(mols, legends=list(compounds.keys()), molsPerRow=4)
    img.save("vkor_compounds.png")
    print("\nCompound structures saved as 'vkor_compounds.png'")
except Exception as e:
    print(f"Error saving compound structures: {e}")


# =============================================
# 4. Tanimoto Similarity to 4-Hydroxycoumarin and Quercetin
# =============================================
def compute_similarity_to_reference(df, reference_name):
    if reference_name not in df["Name"].values:
        raise ValueError(f"Reference compound {reference_name} not found in DataFrame")

    ref_smiles = df[df["Name"] == reference_name]["SMILES"].values[0]
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    if ref_mol is None:
        raise ValueError(f"Invalid SMILES for reference {reference_name}")
    ref_fp = Chem.RDKFingerprint(ref_mol)

    similarities = []
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["SMILES"])
        if mol is None:
            similarities.append(0.0)
            continue
        fp = Chem.RDKFingerprint(mol)
        sim = Chem.DataStructs.TanimotoSimilarity(ref_fp, fp)
        similarities.append(sim)

    df[f"Similarity to {reference_name}"] = similarities
    return df.sort_values(f"Similarity to {reference_name}", ascending=False)


# Compute similarity to 4-Hydroxycoumarin (Control)
df = compute_similarity_to_reference(df, reference_name="4-Hydroxycoumarin (Control)")
print("\nTanimoto Similarity to 4-Hydroxycoumarin (Control):")
print(df[["Name", "Similarity to 4-Hydroxycoumarin (Control)"]].to_string(index=False))

# Compute similarity to Quercetin
df = compute_similarity_to_reference(df, reference_name="Quercetin")
print("\nTanimoto Similarity to Quercetin:")
print(df[["Name", "Similarity to Quercetin"]].to_string(index=False))


# =============================================
# 5. Rule-of-5 and Drug-Likeness Filters
# =============================================
def check_drug_likeness(df):
    def is_lipinski_compliant(row):
        return (row["Molecular Weight"] <= 500 and
                row["LogP"] <= 5 and
                row["H-Bond Donors"] <= 5 and
                row["H-Bond Acceptors"] <= 10)

    df["Lipinski Compliant"] = df.apply(is_lipinski_compliant, axis=1)
    return df


df = check_drug_likeness(df)
print("\nLipinski Rule-of-5 Compliance:")
print(df[["Name", "Lipinski Compliant"]].to_string(index=False))

# =============================================
# 6. Plot Results
# =============================================
try:
    plt.figure(figsize=(12, 6))
    plt.bar(df["Name"], df["Similarity to 4-Hydroxycoumarin (Control)"], color='skyblue', edgecolor='black')
    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.title("Tanimoto Similarity to 4-Hydroxycoumarin (VKOR Inhibitor)", fontsize=12)
    plt.ylabel("Similarity Score (0-1)", fontsize=10)
    plt.xlabel("Compound", fontsize=10)
    plt.tight_layout()
    plt.savefig("similarity_to_4_hydroxycoumarin.png", dpi=300)
    print("\nSimilarity plot saved as 'similarity_to_4_hydroxycoumarin.png'")
except Exception as e:
    print(f"Error saving similarity plot: {e}")

try:
    plt.figure(figsize=(12, 6))
    plt.bar(df["Name"], df["Similarity to Quercetin"], color='lightgreen', edgecolor='black')
    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.title("Tanimoto Similarity to Quercetin", fontsize=12)
    plt.ylabel("Similarity Score (0-1)", fontsize=10)
    plt.xlabel("Compound", fontsize=10)
    plt.tight_layout()
    plt.savefig("similarity_to_quercetin.png", dpi=300)
    print("\nSimilarity plot saved as 'similarity_to_quercetin.png'")
except Exception as e:
    print(f"Error saving similarity plot: {e}")