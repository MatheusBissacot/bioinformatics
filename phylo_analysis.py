import argparse
from datetime import datetime
from utils import (
    preprocess_sequence,
    run_blast,
    parse_blast_results,
    perform_msa,
    build_phylogenetic_tree,
    build_tree_image,
    generate_report,
    create_report, 
    parse_fasta, 
    parse_aln_file
)

def main():
    parser = argparse.ArgumentParser(description='Phylogenetic Analysis Tool')
    parser.add_argument('input_file', help='Path to the input FASTA sequence file')
    parser.add_argument('num_species', type=int, help='Number of species to include in the analysis')
    args = parser.parse_args()

    print("[*] Preprocessing sequence...")
    protein_sequence, id, frame = preprocess_sequence(args.input_file)

    print("[*] Running BLAST search...")
    blast_xml = run_blast(protein_sequence, args.num_species)

    print("[*] Parsing BLAST results...")
    homologs_fasta = parse_blast_results(blast_xml, args.num_species, protein_sequence)

    print("[*] Performing multiple sequence alignment...")
    msa_file = perform_msa(homologs_fasta)

    print("[*] Building phylogenetic tree...")
    tree_file = build_phylogenetic_tree(msa_file)
    tree = build_tree_image(tree_file)

    print("[*] Generating report...")
    generate_report(args.input_file, homologs_fasta, msa_file, tree_file)

    homologs = parse_fasta(homologs_fasta, protein_sequence)
    msa_file = parse_aln_file("results/alignment.aln")
    
    data = {
        'gene_name': id,  # Nome do gene
        'best_frame': frame,  # Frame de leitura selecionado
        'translated_sequence': protein_sequence,  # Sequência traduzida simulada
        'blast_hits': homologs,  # Hits do BLAST
        'blast_hits_number': len(homologs) - 1,
        'tree_image_path': tree,  # Caminho para a árvore filogenética gerada
        'clustal_alignment': msa_file,
    }

    create_report(data)

    print("[+] Analysis complete! Check the 'results/' directory.")

if __name__ == '__main__':
    main()
