import numpy as np

# Data: List of dictionaries containing molecule names and descriptors
molecules = [
    {"name": "Quercetin", "MW": 302.24, "LogP": 1.99, "HBD": 5, "HBA": 7, "AR": 3, "ΔG": -10.77},
    {"name": "Gallic acid", "MW": 170.12, "LogP": 0.50, "HBD": 4, "HBA": 4, "AR": 1, "ΔG": -6.42},
    {"name": "Chlorogenic acid", "MW": 354.31, "LogP": -0.64, "HBD": 6, "HBA": 8, "AR": 1, "ΔG": -11.60},
    {"name": "Catechol", "MW": 110.11, "LogP": 1.09, "HBD": 2, "HBA": 2, "AR": 1, "ΔG": -5.39},
    {"name": "Caffeic acid", "MW": 180.16, "LogP": 1.19, "HBD": 3, "HBA": 3, "AR": 1, "ΔG": -7.07},
    {"name": "Coumaric acid", "MW": 164.16, "LogP": 1.49, "HBD": 2, "HBA": 2, "AR": 1, "ΔG": -6.04},
    {"name": "Rosmarinic acid", "MW": 360.32, "LogP": 1.76, "HBD": 5, "HBA": 7, "AR": 2, "ΔG": -10.91},
    {"name": "Warfarin", "MW": 308.33, "LogP": 3.61, "HBD": 1, "HBA": 4, "AR": 3, "ΔG": -10.11}
]

# Descriptors to use
descriptors = ["MW", "LogP", "HBD","HBA", "AR", "ΔG"]


# Normalize descriptors
def normalize_descriptors(molecules, descriptors):
    normalized_data = []
    for descriptor in descriptors:
        values = [mol[descriptor] for mol in molecules]
        min_val, max_val = min(values), max(values)
        range_val = max_val - min_val if max_val != min_val else 1  # Avoid division by zero
        for mol in molecules:
            mol[f"{descriptor}_norm"] = (mol[descriptor] - min_val) / range_val
    return molecules


# Calculate Tanimoto similarity
def tanimoto_similarity(mol1, mol2, descriptors):
    # Extract normalized descriptor values
    vec1 = np.array([mol1[f"{d}_norm"] for d in descriptors])
    vec2 = np.array([mol2[f"{d}_norm"] for d in descriptors])

    # Calculate dot product and sums of squares
    dot_product = np.sum(vec1 * vec2)
    sum_sq1 = np.sum(vec1 ** 2)
    sum_sq2 = np.sum(vec2 ** 2)

    # Tanimoto coefficient
    if sum_sq1 + sum_sq2 - dot_product == 0:  # Avoid division by zero
        return 0.0
    return dot_product / (sum_sq1 + sum_sq2 - dot_product)


# Main execution
def main():
    # Normalize descriptors
    normalized_molecules = normalize_descriptors(molecules, descriptors)

    # Find warfarin
    warfarin = next(mol for mol in normalized_molecules if mol["name"] == "Warfarin")

    # Calculate similarity scores
    similarity_scores = []
    for mol in normalized_molecules:
        if mol["name"] != "Warfarin":
            score = tanimoto_similarity(warfarin, mol, descriptors)
            similarity_scores.append({"name": mol["name"], "similarity": score})

    # Sort by similarity score (descending)
    similarity_scores.sort(key=lambda x: x["similarity"], reverse=True)

    # Print results
    print("Tanimoto Similarity Scores relative to Warfarin:")
    print("-" * 45)
    print(f"{'Molecule':<20} {'Similarity Score':<15}")
    print("-" * 45)
    for item in similarity_scores:
        print(f"{item['name']:<20} {item['similarity']:.3f}")


if __name__ == "__main__":
    main()