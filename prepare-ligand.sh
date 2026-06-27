scrub.py $1.smi -o $1.sdf --skip_tautomers --ph_low 7.35 --ph_high 7.45
mk_prepare_ligand.py -i $1.sdf
