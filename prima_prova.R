library(bio3d)
library(igraph)
library(dplyr) 
library(ggplot2)


#demo("md")

#trj_file <- read.dcd("/Users/lorenzosisti/lysozyme_tutorial_gromacs/prova_traiettoria_lisozima.dcd")
trj_file <- read.dcd("/Users/lorenzosisti/lysozyme_tutorial_gromacs/lysozyme_Protein.dcd")
pdb_file <- read.pdb("/Users/lorenzosisti/lysozyme_tutorial_gromacs/1AKI_clean.pdb")

ca_indexes <- atom.select(pdb_file, elety = "CA")

xyz <- fit.xyz(fixed = pdb_file$xyz, mobile = trj_file, fixed.inds = ca_indexes$xyz, mobile.inds = ca_indexes$xyz)

rmsd <- rmsd(xyz[1, ca_indexes$xyz], xyz[, ca_indexes$xyz])
plot(rmsd, typ = "l", ylab = "RMSD", xlab = "Frame No.")
#points(lowess(rmsd), typ = "l", col = "red", lty = 2, lwd = 2)
summary(rmsd)

rmsf <- rmsf(xyz[, ca_indexes$xyz])
plot(rmsf, ylab = "RMSF", xlab = "Residue Position", typ = "l")

pc <- pca.xyz(xyz[, ca_indexes$xyz])
plot(pc, col = bwr.colors(nrow(xyz)))

# Ora faccio la parte con i grafi, in particolare implemento una funzione che calcola i nodi sui centroidi delle side chains

# -----------------------------------------------------------------
# Inizio della nuova sezione: Analisi di Rete sulla Traiettoria
# -----------------------------------------------------------------

# Mantengo la tua funzione per i descrittori di rete ESATTAMENTE come l'hai fornita.
# Verrà usata per ogni frame.
Local_Network_Des_Bin <- function(g){
  BetweennessCentrality <- betweenness(g, v= V(g), directed = FALSE, normalized = TRUE)
  ClosenessCentrality <- closeness(g, v= V(g), normalized = TRUE)
  ShortestPath <- mean_distance(g)
  Degree <- degree(g)
  ClusteringCoefficient <- transitivity(g, type = "global")
  Density <- edge_density(g, loops = FALSE)
  Modularity <- modularity(g, membership(cluster_walktrap(g)))
  
  df_des <- as.data.frame(cbind(BetweennessCentrality, ClosenessCentrality, ShortestPath, Degree, ClusteringCoefficient, Density, Modularity))
  return((df_des))
}

# Non useremo le altre due funzioni ('compute_graph_energy' e 
# 'geometrical_topological_descriptors_old') perché la loro logica 
# sarà integrata direttamente nel ciclo sui frame.

## 1. Definizione dei Nodi (Centroidi)

# Imposta un cutoff di distanza per definire un "contatto" (edge nel grafo)
# Un valore comune per i contatti tra side-chain è ~7 Ångström.
DistCutoff <- 8.5 

#cat("Definizione dei nodi (centroidi side chain / CA per GLY)...\n")

pdb_atoms <- pdb_file$atom
ca_atoms <- atom.select(pdb_file, elety = "CA")

# Ottieni l'elenco univoco di residui (numero e nome)
unique_residues <- unique(pdb_atoms[ca_atoms$atom, c("resno", "resid")])
n_nodes <- nrow(unique_residues)

# Questa lista conterrà gli *indici xyz* degli atomi che compongono ciascun nodo
node_xyz_indices <- vector("list", n_nodes)
node_names <- rep(NA, n_nodes) # Per salvare i nomi dei residui (es. "LYS1")

for (i in 1:n_nodes) {
  r_num <- unique_residues$resno[i]
  r_id <- unique_residues$resid[i]
  
  node_names[i] <- paste0(r_id, r_num)
  
  if (r_id == "GLY") {
    # Per la Glicina, usiamo il Carbonio Alpha (CA)
    inds <- atom.select(pdb_file, resno = r_num, elety = "CA")
  } else {
    # Per tutti gli altri, usiamo la side chain
    inds <- atom.select(pdb_file, resno = r_num, string = "sidechain")
    
    # Fallback: se una side chain non ha atomi (improbabile, ma sicuro), usa il CA
    #if (length(inds$atom) == 0) {
    #  inds <- atom.select(pdb_file, resno = r_num, elety = "CA")
    #}
  }
  
  # Salviamo gli indici XYZ che useremo per estrarre le coordinate dalla matrice 'xyz'
  node_xyz_indices[[i]] <- inds$xyz
}

cat("Nodi definiti:", n_nodes, "\n")


## 2. Loop sui Frame della Traiettoria

n_frames <- nrow(xyz)
# Creiamo una lista per contenere i data.frame dei risultati di ciascun frame
results_list <- vector("list", n_frames)

cat("Inizio elaborazione dei frame (Totali:", n_frames, ")...\n")

for (f in 1:n_frames) {
  
  # Aggiornamento sullo stato di avanzamento
  if (f %% 100 == 0 || f == 1 || f == n_frames) {
    cat("Elaborazione frame:", f, "/", n_frames, "\n")
  }
  
  # Estrai le coordinate XYZ per il frame corrente
  current_frame_xyz <- xyz[f, ]
  
  # Calcola i centroidi per questo frame
  centroids_f <- matrix(NA, nrow = n_nodes, ncol = 3)
  
  for (i in 1:n_nodes) {
    xyz_inds_for_node_i <- node_xyz_indices[[i]]
    
    # Estrai le coordinate degli atomi per questo nodo
    # Dobbiamo gestirli come matrice (un atomo per riga)
    node_coords_matrix <- matrix(current_frame_xyz[xyz_inds_for_node_i], ncol = 3, byrow = TRUE)
    
    # Calcola il centroide (media delle coordinate x, y, z)
    centroids_f[i, ] <- colMeans(node_coords_matrix)
  }
  
  colnames(centroids_f) <- c("x", "y", "z")
  rownames(centroids_f) <- node_names
  
  # --- Costruzione del Grafo per questo Frame ---
  
  # 1. Matrice delle distanze tra i centroidi
  dist_mat <- as.matrix(dist(centroids_f))
  
  # 2. Matrice di adiacenza binaria (logica presa dalla tua 'compute_graph_energy')
  adj_matrix <- ifelse(dist_mat <= DistCutoff, 1, 0)
  diag(adj_matrix) <- 0 # Un nodo non ha un edge con se stesso
  
  # 3. Calcolo Graph Energy (logica presa dalla tua 'compute_graph_energy')
  eig <- eigen(adj_matrix, only.values = TRUE)
  graph_energy <- sum(abs(eig$values))
  
  # 4. Creazione dell'oggetto grafo 'igraph'
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  
  # --- Calcolo Descrittori di Rete ---
  
  # 5. Usa la tua funzione per calcolare i descrittori
  df_des <- Local_Network_Des_Bin(g)
  
  # 6. Aggiungi le informazioni supplementari
  df_des$GraphEnergy <- graph_energy # Aggiungiamo l'energia calcolata
  df_des$frame <- f                   # Il numero del frame corrente
  df_des$residue <- node_names        # I nomi dei residui (nodi)
  
  # 7. Salva il data.frame di questo frame nella lista
  results_list[[f]] <- df_des
}

cat("Elaborazione di tutti i frame completata.\n")


## 3. Combinazione dei Risultati

cat("Combinazione dei risultati in un unico data.frame...\n")

# Carica 'dplyr' per 'bind_rows', che è più efficiente di do.call(rbind, ...)
final_results_df <- bind_rows(results_list)

cat("Fatto.\n\n")

# -----------------------------------------------------------------
# Fine della nuova sezione
# -----------------------------------------------------------------

## 4. Analisi dei Risultati

# Ora hai un singolo, grande data.frame 'final_results_df'
# con (numero_di_frame * numero_di_residui) righe.
# Ogni riga contiene i descrittori di rete per un residuo specifico 
# in un frame specifico.

print(head(final_results_df))
summary(final_results_df)

# Puoi usarlo per plottare come un descrittore cambia nel tempo per un residuo.
# Ad esempio, come cambia la Betweenness del primo residuo (es. LYS1)

# Carica ggplot2 per grafici migliori

# Estrai i dati per il primo residuo (puoi cambiare 'node_names[1]' con "LYS1" o quello che ti serve)
residue_to_plot <- node_names[1]
data_residue_1 <- subset(final_results_df, residue == residue_to_plot)

ggplot(data_residue_1, aes(x = frame, y = BetweennessCentrality)) +
  geom_line(color = "blue") +
  labs(title = paste("Betweenness Centrality per", residue_to_plot),
       x = "Numero Frame",
       y = "Betweenness (normalizzata)") +
  theme_minimal()


# Oppure, puoi plottare un descrittore MEDIO sull'intera proteina nel tempo
# Ad esempio, il Degree medio per frame
avg_degree_per_frame <- final_results_df %>%
  group_by(frame) %>%
  summarize(MeanDegree = mean(Degree),
            GlobalDensity = first(Density)) # 'Density' è già globale

ggplot(avg_degree_per_frame, aes(x = frame, y = MeanDegree)) +
  geom_line(color = "red") +
  labs(title = "Degree Medio per Frame",
       x = "Numero Frame",
       y = "Degree Medio") +
  theme_minimal()
