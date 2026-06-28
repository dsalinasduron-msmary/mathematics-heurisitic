tar -zxvf charmm36-feb2026_cgenff-5.0.ff.tgz

gmx pdb2gmx -f 1iep_no_ligands.pdb -o 1iep_processed.gro -ter
