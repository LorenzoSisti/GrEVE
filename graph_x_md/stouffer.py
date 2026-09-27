# Statistical assessment of the replicas BC Z-score

import os
import numpy as np
import pandas as pd

# ==========================================
# 1. PERCORSI DEI FILE
# ==========================================
data_dir = ""

file_r1 = os.path.join(data_dir, "replica_1.csv")
file_r2 = os.path.join(data_dir, "replica_2.csv")
file_r3 = os.path.join(data_dir, "replica_3.csv")

# ==========================================
# 2. CARICAMENTO E MERGE DEI DATI
# ==========================================
df1 = pd.read_csv(file_r1)
df2 = pd.read_csv(file_r2)
df3 = pd.read_csv(file_r3)

merged = df1[["residue", "z_score"]].rename(columns={"z_score": "z_score_R1"})
merged = merged.merge(
    df2[["residue", "z_score"]].rename(columns={"z_score": "z_score_R2"}),
    on="residue",
)
merged = merged.merge(
    df3[["residue", "z_score"]].rename(columns={"z_score": "z_score_R3"}),
    on="residue",
)

# ==========================================
# 3. CALCOLO METRICHE (STOUFFER, ABS E STD)
# ==========================================
# Stouffer Meta-Z = (Z1 + Z2 + Z3) / sqrt(3)
merged["stouffer_z"] = (
    merged["z_score_R1"] + merged["z_score_R2"] + merged["z_score_R3"]
) / np.sqrt(3)

# Valore assoluto per catturare deviazioni estreme sia positive che negative
merged["abs_stouffer_z"] = merged["stouffer_z"].abs()

# Deviazione standard tra le tre repliche
merged["inter_replica_std"] = merged[
    ["z_score_R1", "z_score_R2", "z_score_R3"]
].std(axis=1)

# Ordiniamo per ampiezza assoluta (|Z| decrescente)
merged = merged.sort_values(by="abs_stouffer_z", ascending=False).reset_index(
    drop=True
)

# ==========================================
# 4. CALCOLO PERCENTILI ASSOLUTI E FILTRAGGIO
# ==========================================
p95 = np.percentile(merged["abs_stouffer_z"], 95.0)
p97_5 = np.percentile(merged["abs_stouffer_z"], 97.5)
p99 = np.percentile(merged["abs_stouffer_z"], 99.0)

print("=== SOGLIE PERCENTILI (AMPIEZZA ASSOLUTA |Z|) ===")
print(f"Top 5.0% (95° percentile)  : |Z| >= {p95:.4f}")
print(f"Top 2.5% (97.5° percentile): |Z| >= {p97_5:.4f}")
print(f"Top 1.0% (99° percentile)  : |Z| >= {p99:.4f}\n")

# Estrazione del Top 2.5% in valore assoluto (include sia valori + che -)
top_2_5 = merged[merged["abs_stouffer_z"] >= p97_5].copy()

print(f"Residui totali analizzati: {len(merged)}")
print(f"Residui selezionati (Top 2.5% assoluto): {len(top_2_5)}\n")
print("I 10 residui con variazione più estrema (|Z| più alto):")
print(top_2_5.head(10).to_string(index=False))

# ==========================================
# 5. SALVATAGGIO DEI DATI IN CSV
# ==========================================
output_all = os.path.join(data_dir, "stouffer_all_residues.csv")
output_top = os.path.join(data_dir, "stouffer_top_2_5_percent.csv")

merged.to_csv(output_all, index=False)
top_2_5.to_csv(output_top, index=False)

print(f"\n[OK] Risultati salvati in:\n - {output_all}\n - {output_top}")