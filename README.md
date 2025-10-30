# GrEVE: GRaph Evolution Visualization and Evaluation

A pipeline developed during my internship at the [Structural Bioinformatics Group at Sapienza University of Rome](https://schubert.bio.uniroma1.it/index.html).

Per strippare una traiettoria dall'acqua e da eventuali idrogeni:

```
gmx trjconv -f md_0_1_noPBC.xtc -s md_0_1.tpr -o lysozyme_Protein.xtc
```

Per convertire un xtc in dcd:

```
pip install mdtraj
mdconvert -o lysozyme_Protein.dcd  lysozyme_Protein.xtc
```
Work in progress...
