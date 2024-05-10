for file in *.pdbqt; do
    grep "ATOM" "$file" > "${file%.pdbqt}_atom_only.pdbqt"
done
