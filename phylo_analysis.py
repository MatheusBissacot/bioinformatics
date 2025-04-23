import argparse
from datetime import datetime
from utils import (
    preprocess_sequence,
    run_blast,
    parse_blast_results,
    perform_msa,
    build_phylogenetic_tree,
    generate_report
)

def main():
    parser = argparse.ArgumentParser(description='Phylogenetic Analysis Tool')
    parser.add_argument('input_file', help='Path to the input FASTA sequence file')
    parser.add_argument('num_species', type=int, help='Number of species to include in the analysis')
    args = parser.parse_args()

    print("[*] Preprocessing sequence...")
    protein_sequence = preprocess_sequence(args.input_file)

    print("[*] Running BLAST search...")
    blast_xml = run_blast(protein_sequence)

    print("[*] Parsing BLAST results...")
    homologs_fasta = parse_blast_results(blast_xml, args.num_species, protein_sequence)

    print("[*] Performing multiple sequence alignment...")
    msa_file = perform_msa(homologs_fasta)

    print("[*] Building phylogenetic tree...")
    tree_file = build_phylogenetic_tree(msa_file)

    print("[*] Generating report...")
    generate_report(args.input_file, homologs_fasta, msa_file, tree_file)

    data = {
        'gene_name': protein_sequence.id,  # Nome do gene
        'best_frame': '+1',  # Frame de leitura selecionado
        'translated_sequence': protein_sequence,  # Sequência traduzida simulada
        'blast_hits': homologs_fasta,  # Hits do BLAST
        'tree_image_path': tree_file,  # Caminho para a árvore filogenética gerada
        'current_year': datetime.now().year  # Ano atual
    }

    print("[+] Analysis complete! Check the 'results/' directory.")

if __name__ == '__main__':
    main()
