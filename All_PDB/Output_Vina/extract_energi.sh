for file in *.pdbqt; do
    filename=$(basename -- "$file")
    filename="${filename%.*}"
    energy=$(sed '2q;d' "$file" | awk '{print $4}')
    echo "$filename,$energy" >> fichier.csv
done
