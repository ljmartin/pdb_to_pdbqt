# pdb_to_pdbqt
 download and prepare a pdb file for docking wth Vina-flavoured algorithms

# run notebook
user must input the PDB code (5WIU in the test case), and the ligand resname (AQD in the test case). If performing the optional step of parameterizing the ligand with GAFF (via openmmforcefields), user must also input the ligand isomeric SMILES code.

outputs are:
- 5wiu.pdb.gz - not really an output, it's just the compressed pdb file from the PDB
- protein.pdb - protein atoms only, aligned to the ligand's principal moment of inertia
- ligand.pdb - ligand only, aligned to the ligand's principal moment of inertia
- proteinH.pdb - protein only, missing atoms added by PDBFixer, and hydrogen's minimized by OpenMM
- ligandH.pdb - ligand only, bond topology fixed by RDKit
- combined.pdb - protein, ligand, both with hydrogens, before minimizing.
- gaff_minimized.pdb - protein, ligand, both with hydrogens and minimized in vacuum without any constraints using GAFF-parameterized ligand

After this, one can generate the docking input file by:
```
obabel proteinH.pdb -xr -O proteinH.pdbqt
```

or use ODDT [MolToMolToPDBQTBlock](https://oddt.readthedocs.io/en/latest/rst/oddt.toolkits.extras.rdkit.html#oddt.toolkits.extras.rdkit.MolToPDBQTBlock)
