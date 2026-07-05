# add hydrogens
scrub.py $1.smi -o $1.sdf --skip_tautomers --ph_low 7.35 --ph_high 7.45

# make pdbqt with -w option so that waters are added to ligand
mk_prepare_ligand.py -w -i $1.sdf 
