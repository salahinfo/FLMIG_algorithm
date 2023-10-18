from sklearn.metrics import normalized_mutual_info_score

# Example labels (replace with your actual labels)
labels_true = [0, 1, 2, 0, 1, 2]
labels_pred = [4, 0, 2, 4, 0, 2]

# Calculate NMI
nmi = normalized_mutual_info_score(labels_true, labels_pred)
print("Normalized Mutual Information (NMI):", nmi)
