# Required packages: dynetan, networkx, mdanalysis

from __future__ import annotations

import pickle
from pathlib import Path
import numpy as np
import networkx as nx
import dynetan
from dynetan.proctraj import DNAproc

# --- Load the two states trajectory (necessary for delta analyses) -------------------------------------------------
APO_TOPOLOGY = "/Users/lorenzosisti/GrEVE/invisible_data/apo.pdb" 
APO_TRAJECTORY = "/Users/lorenzosisti/GrEVE/invisible_data/apo_nowt.xtc"
APO_SEG_IDS = ["A", "B"]  # Segment ID(s) (chain A and chain B)
HOLO_TOPOLOGY = "/Users/lorenzosisti/GrEVE/invisible_data/holo.pdb"
HOLO_TRAJECTORY = "/Users/lorenzosisti/GrEVE/invisible_data/holo_nowt.xtc"
HOLO_SEG_IDS = ["A", "B"]

# --- Parameters ---------------------------------------------
CORRELATION_STRIDE = 1      
CUTOFF_DIST = 4.5           # Cutoff distance to define a contact (in Angstroms)
CONTACT_PERSISTENCE = 0.75  # A contact is defined if it is present in at least 75% of the frames
N_CORES = 7                 # Core number for parallel GC computing
WEIGHT_EPSILON = 1e-6       # Floor to avoid -log(0)

OUTPUT_DIR = "/Users/lorenzosisti/GrEVE/invisible_data/all_atom_weighted_graphs"


# ----------------------------------------------------------------------
# Define a function to generate the All-Atom node groups (1 node per heavy atom)
# ----------------------------------------------------------------------

def build_all_atom_groups(dnap: DNAproc, output_txt: str = "node_groups_inspection.txt") -> dict[str, dict[str, set[str]]]:
    """
    Genera il dizionario All-Atom normalizzando i nomi degli atomi C-terminali 
    (OT1/OXT -> O) per non perdere l'atomo O del backbone nei calcoli di rete.
    Salva un report leggibile in formato .txt per l'ispezione manuale.
    """
    u = dnap.getU()

    # 1. Normalizza i nomi degli atomi C-terminali nell'universo MDAnalysis
    for atom in u.atoms:
        if atom.name in ("OT1", "O1", "OXT"):
            atom.name = "O"

    # 2. Mappa i residui per tipologia (resname)
    res_type_map = {}
    for res in u.residues:
        resname = res.resname
        if resname not in res_type_map:
            res_type_map[resname] = []
        res_type_map[resname].append(res)

    # 3. Seleziona gli atomi pesanti standard (escludendo l'atomo terminale extra OT2/O2)
    node_groups = {}
    for resname, residues in res_type_map.items():
        heavy_sets = []
        for r in residues:
            heavy_atoms = r.atoms.select_atoms("not element H and not name OT2 O2")
            heavy_sets.append({a.name for a in heavy_atoms})

        if not heavy_sets:
            continue

        common_atoms = set.intersection(*heavy_sets)
        node_groups[resname] = {atom_name: {atom_name} for atom_name in common_atoms}

    # 4. Salva il dizionario formattato in un file .txt per ispezione
    with open(output_txt, "w") as f:
        for resname, nodes in sorted(node_groups.items()):
            f.write(f"=== AMMINOACIDO: {resname} ({len(nodes)} nodi) ===\n")
            for node_name, atoms in sorted(nodes.items()):
                atoms_str = ", ".join(sorted(atoms))
                f.write(f"  • Nodo '{node_name}': atomi = [{atoms_str}]\n")
            f.write("\n")

    print(f" -> Dizionario dei nodi salvato per ispezione in: {output_txt}")

    return node_groups


# ----------------------------------------------------------------------
# Define a function to run the dynetan pipeline for a single trajectory
# ----------------------------------------------------------------------
def run_dynetan_pipeline(
    topology: str,
    trajectory: str,
    seg_ids: list[str],
    stride: int = 1,
    cutoff: float = 4.5, # For contact definition between CA to 4.5 Å. For all-atoms simulations it would be 4.5
    persistence: float = 0.75,
    ncores: int = 4,
) -> DNAproc:
    """
    DNAproc:Dynamic Network Analysis processing class. It contains the
    infrastructure to carry out the data analysis.
    1. Parameters initialisation
    2. Loading topology and trajectory
    3. Alignment
    4. Contact calculation and filtering
    5. Generalized Correlation calculation
    """
    dnap = DNAproc()

    dnap.setNumWinds(1)  # 1 single window for the entire trajectory
    dnap.setCutoffDist(cutoff) 
    dnap.setContactPersistence(persistence)
    dnap.setSegIDs(seg_ids)

    # System loading
    print(f"System loading: {topology}, {trajectory}")
    dnap.loadSystem(topology, trajectory)

    # Decomment the following line if you want to use all-atom representation instead of CA representation
    dnap.setNodeGroups(build_all_atom_groups(dnap))

    # Preparation of the network nodes (1 node per protein residue / C-alpha in this case)
    dnap.prepareNetwork()

    # Alignment of the trajectory (in memory if possible)
    print("Alignment of the trajectory...")
    dnap.alignTraj(inMemory=True)

    # Contact calculation 
    print(f"Contact calculation (cutoff={cutoff}Å, stride={stride})...")
    dnap.findContacts(stride=stride)

    # Filtering of intra-residue and isolated contacts
    dnap.filterContacts(
        notSameRes=True, notConsecutiveRes=False, removeIsolatedNodes=False
    )

    # Calculation of the Generalized Correlation (GC) in parallel
    print(f"Generalized Correlation calculation on {ncores} cores...")
    dnap.calcCor(ncores=ncores, verbose=1)

    return dnap


# ----------------------------------------------------------------------
# Processing and saving the NetworkX graph with GC and weight attributes
# ----------------------------------------------------------------------
def process_and_save_networkx(
    dnap: DNAproc, output_dir: str, label: str, epsilon: float = 1e-6
) -> None:
    """
    Extracts the NetworkX graph generated by dynetan, assigns the 'GC' and 'weight' attributes to each edge, and saves it in both .graphml and .pkl formats.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Save matrix in numpy (.npy)
    gc_matrix = dnap.corrMatAll[0]  # N x N matrix of the window (N = number of nodes) 
    npy_path = out_dir / f"{label}_GC_matrix.npy"
    np.save(npy_path, gc_matrix)
    print(f" -> Matrix saved in: {npy_path}")

    # 2. Create a graph from the correlation matrix
    dnap.calcGraphInfo()
    G = dnap.nxGraphs[0]

    # Assign resid/resname/chain as node attributes, so the saved graph
    # is self-contained (no need for the dnap object to know "who is who").
    resid_map = {i: int(atom.resid) for i, atom in enumerate(dnap.nodesAtmSel)}
    resname_map = {i: atom.resname for i, atom in enumerate(dnap.nodesAtmSel)}

    # Extract the chain from the segid or chainID of the MDAnalysis atom
    chain_map = {
        i: getattr(atom, "segid", getattr(atom, "chainID", "A"))
        for i, atom in enumerate(dnap.nodesAtmSel)
    }

    nx.set_node_attributes(G, resid_map, "resid")
    nx.set_node_attributes(G, resname_map, "resname")
    nx.set_node_attributes(G, chain_map, "chain")

    # The following line saves a CSV file with node information (node_id, chain, resid, resname, atom_name).
    # This was necessary as in the original documentation it is not well explained how to select CA over all-atoms.
    # Thus, I used this .csv file to check whether I was running a CG or all-atom simulation.
    open(out_dir / f"{label}_nodes.csv", "w").write("node_id,chain,resid,resname,atom_name\n" + "\n".join(f"{i},{getattr(a, 'segid', getattr(a, 'chainID', 'A'))},{a.resid},{a.resname},{a.name}" for i, a in enumerate(dnap.nodesAtmSel)))

    # 3. Compute 'GC' and 'weight' = -log(GC) directly on the edges
    for u, v in G.edges():
        gc_val = float(gc_matrix[u, v])
        gc_clamped = max(gc_val, epsilon)
        weight_val = float(-np.log(gc_clamped))

        G[u][v]["GC"] = gc_val
        G[u][v]["weight"] = weight_val

    # 4. Save in Pickle (.pkl)
    pickle_path = out_dir / f"{label}_weighted.pkl"
    with open(pickle_path, "wb") as f:
        pickle.dump(G, f)
    print(f" -> Graph NetworkX saved (.pkl): {pickle_path}")

    # 5. Save in GraphML (.graphml)
    graphml_path = out_dir / f"{label}_weighted.graphml"
    G_graphml = G.copy()

    for _, data in G_graphml.nodes(data=True):
        for k, val in list(data.items()):
            if isinstance(val, (set, list, tuple, dict)):
                data[k] = str(val)

    for _, _, data in G_graphml.edges(data=True):
        for k, val in list(data.items()):
            if isinstance(val, (set, list, tuple, dict)):
                data[k] = str(val)

    nx.write_graphml(G_graphml, str(graphml_path))
    print(f" -> Graph NetworkX saved (.graphml): {graphml_path}")


# ----------------------------------------------------------------------
# Run the pipeline
# ----------------------------------------------------------------------
def process_system(
    topology: str,
    trajectory: str,
    seg_ids: list[str],
    label: str,
) -> None:
    print(f"\n=================== System: {label} ===================")

    # Run dynetan pipeline
    dnap = run_dynetan_pipeline(
        topology=topology,
        trajectory=trajectory,
        seg_ids=seg_ids,
        stride=CORRELATION_STRIDE,
        cutoff=CUTOFF_DIST,
        persistence=CONTACT_PERSISTENCE,
        ncores=N_CORES,
    )

    # Process and save the graph in NetworkX
    process_and_save_networkx(
        dnap=dnap,
        output_dir=OUTPUT_DIR,
        label=label,
        epsilon=WEIGHT_EPSILON,
    )


def main() -> None:
    systems = [
        dict(
            topology=APO_TOPOLOGY,
            trajectory=APO_TRAJECTORY,
            seg_ids=APO_SEG_IDS,
            label="apo_open",
        ),
        dict(
            topology=HOLO_TOPOLOGY,
            trajectory=HOLO_TRAJECTORY,
            seg_ids=HOLO_SEG_IDS,
            label="holo_closed",
        ),
    ]

    for sys_info in systems:
        process_system(**sys_info)

    print("\nGood job! You completed the single-trajectory pipeline.")


if __name__ == "__main__":
    main()