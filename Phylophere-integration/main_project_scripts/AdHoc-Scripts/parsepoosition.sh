while IFS=$' ' read -r species amino_acids; do
    # Extract species name
    species_name=$(echo "$species")

    # Extract the first 5 amino acids
    first_five_amino_acids=$(echo "$amino_acids" | cut -c 1-5)

    # Print or use the extracted information
    echo "Species Name: $species_name"
    echo "First Five Amino Acids: $first_five_amino_acids"
done < ./2.Alignments/Primate_alignments/HMBOX1.Homo_sapiens.filter2.phy

