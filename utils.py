# utils.py

from Bio import SeqIO, AlignIO, Seq, Phylo
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
import requests
import os
import uuid
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib
import matplotlib.pyplot as plt
from jinja2 import Environment, FileSystemLoader

RESULTS_DIR = "results"

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

ensure_dir(RESULTS_DIR)

def preprocess_sequence(input_file):
    """Reads the input FASTA, finds the longest ORF if DNA, returns protein sequence."""
    record = SeqIO.read(input_file, "fasta")

    #Some sequences from NCBI GeneBank contain letter 'N', which illustrates that these nucleotide bases are not 
    #deciphered correctly, leaving an unidentified nucleotide. 
    #For this we decid to remove the N's basead in this https://www.researchgate.net/post/How_to_handle_N_in_Nucleotide_Genes_Sequences_retrieved_from_NCBI_GeneBank 
    seq = Seq(str(record.seq).replace("N", "")).upper()
    seq = seq[:len(seq) - len(seq) % 3]
    
    # If DNA, translate to protein
    if set(seq.upper()).issubset({"A", "T", "C", "G"}):
        # Generate 6 frames
        frames = [
            seq.translate(to_stop=False),
            seq[1:].translate(to_stop=False),
            seq[2:].translate(to_stop=False),
            seq.reverse_complement().translate(to_stop=False),
            seq.reverse_complement()[1:].translate(to_stop=False),
            seq.reverse_complement()[2:].translate(to_stop=False),
        ]

        longest_peptide = ""
        id = 0
        for idx, frame in enumerate(frames):
            peptides = str(frame).split("*")  # split at stop codon
            longest_in_frame = max(peptides, key=len)
            if len(longest_in_frame) > len(longest_peptide):
                longest_peptide = longest_in_frame
                id = idx

        print(f"[+] Preprocessing complete. Longest ORF length: {len(longest_peptide)} aa. Frame number: {id+1}")
        print("Please choose the protein sequence:")
        for i, prot in enumerate(frames):
            print(f"Frame {i+1}: {prot}")
        try:
            choice = int(input("Enter frame number: "))
            if not 1 <= choice <= 6:
                raise ValueError
        except ValueError:
            print("Invalid frame. Please enter a number between 1 and 6.")
            exit(1)
        return frames[choice-1], record.id, choice
    print("[*] Sequence appears to already be protein or contains invalid bases.")
    return record.seq, record.id, None


def run_blast(protein_sequence, num_species):
    """Runs BLASTP for the given protein sequence and saves XML results."""
    print("[*] Running BLASTP query...")
    hitlist_size = num_species * 3
    result_handle = NCBIWWW.qblast(
        program="blastp",
        database="nr",
        sequence=protein_sequence,
        hitlist_size=hitlist_size,
        expect=1e-5,
        format_type="XML"
    )

    output_file = os.path.join(RESULTS_DIR, f"blast_result_{uuid.uuid4().hex}.xml")
    with open(output_file, "w") as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()

    print(f"[+] BLAST results saved to {output_file}")
    return output_file

def parse_blast_results(blast_xml, num_species, human_protein_seq):
    """Parses BLAST XML and saves top unique species sequences to FASTA."""
    print("[*] Parsing BLAST results...")
    species_sequences = {}
    with open(blast_xml) as result_handle:
        blast_record = NCBIXML.read(result_handle)

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                desc = alignment.hit_def
                if "[" in desc and "]" in desc:
                    species_name = desc.split("[")[-1].strip("]")
                    if species_name not in species_sequences:
                        species_sequences[species_name] = alignment.hsps[0].sbjct

                if len(species_sequences) >= num_species:
                    break
            if len(species_sequences) >= num_species:
                break

    if not species_sequences:
        raise ValueError("No species found in BLAST results. Try increasing hitlist size or check input.")

    output_fasta = os.path.join(RESULTS_DIR, "homologs.fasta")
    with open(output_fasta, "w") as f:
        f.write(">Human_input\n")
        f.write(f"{human_protein_seq}\n")

        for idx, (species, seq) in enumerate(species_sequences.items(), 1):
            f.write(f">seq{idx}_{species.replace(' ', '_')}\n")
            f.write(f"{seq}\n")
    return output_fasta


def perform_msa(input_fasta):
    """Performs Multiple Sequence Alignment (MSA) using the EBI Clustal Omega API."""
    print("[*] Performing MSA via EBI Clustal Omega API...")

    # Read input FASTA
    with open(input_fasta, "r") as file:
        fasta_content = file.read()

    # Submit job
    submit_url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
    params = {
        "email": "your_email@example.com",  
        "sequence": fasta_content,
        "stype": "protein",
        "outfmt": "clustal"
    }
    response = requests.post(submit_url, data=params)

    if response.status_code != 200:
        raise Exception(f"Failed to submit MSA job. Response: {response.text}")

    job_id = response.text.strip()
    print(f"[+] MSA job submitted. Job ID: {job_id}")

    # Check job status
    status_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}"
    while True:
        status_response = requests.get(status_url)
        status = status_response.text.strip()
        if status == "FINISHED":
            print("[+] MSA job finished!")
            break
        elif status in ["RUNNING", "PENDING"]:
            print("[*] Job still running... waiting")
            import time
            time.sleep(3)
        else:
            raise Exception(f"Job failed with status: {status}")

    # Get result
    result_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-clustal"
    result_response = requests.get(result_url)

    msa_output = os.path.join(RESULTS_DIR, "alignment.aln")
    with open(msa_output, "w") as msa_file:
        msa_file.write(result_response.text)

    print(f"[+] MSA saved to {msa_output}")
    return msa_output

def build_phylogenetic_tree(msa_file):
    print("[*] Building phylogenetic tree from MSA...")

    # Read the alignment
    alignment = AlignIO.read(msa_file, "clustal")

    # Calculate distance matrix
    calculator = DistanceCalculator("blosum62")
    distance_matrix = calculator.get_distance(alignment)

    # Construct the tree
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)  # You can also try .upgma()

    # Save the tree
    tree_file = os.path.join(RESULTS_DIR, "phylogenetic_tree.xml")
    Phylo.write(tree, tree_file, "phyloxml")

    print(f"[+] Phylogenetic tree saved to {tree_file}")
    return tree_file

def build_tree_image(tree_file):

    matplotlib.use('Agg')
    output_file = "results/phylogenetic_tree.png"

    tree = Phylo.read(tree_file, "phyloxml")

    for clade in tree.get_terminals():
        if clade.name:
            clade.name = "_".join(clade.name.split("_")[1:])

    fig = plt.figure(figsize=(16, 12))
    axes = fig.add_subplot(1, 1, 1)

    Phylo.draw(tree, do_show=False, axes=axes)

    plt.savefig(output_file, dpi=300)
    print(f"[+] Phylogenetic tree stored at: {output_file}")
    return output_file

def generate_report(input_file, homologs_fasta, msa_file, tree_file):
    """Generates a simple text report of the analysis."""
    report_file = os.path.join(RESULTS_DIR, "report.txt")
    with open(report_file, "w") as f:
        f.write("Phylogenetic Analysis Report\n")
        f.write("===========================\n\n")
        f.write(f"Input FASTA: {input_file}\n")
        f.write(f"Homologs FASTA: {homologs_fasta}\n")
        f.write(f"MSA File: {msa_file}\n")
        f.write(f"Tree File: {tree_file}\n\n")
        f.write("Analysis complete.\n")

    print(f"[+] Report generated at {report_file}")

def parse_fasta(fasta_file, human_protein):
    hits = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        hit = {
            'species':  "_".join(record.id.split("_")[1:]),
            'description': record.seq,
        }
        hits.append(hit)
    return hits

def parse_aln_file(filepath):
    sequences = dict()

    with open(filepath, 'r') as file:
        for line in file:
            line = line.rstrip()
            if not line or line.startswith('CLUSTAL') or line.startswith(' '):
                continue  
            parts = line.split()
            if len(parts) < 2:
                continue  
            seq_id, seq_fragment = parts[0], parts[1]
            if seq_id not in sequences:
                sequences[seq_id] = ''
            sequences[seq_id] += seq_fragment

    max_length = max(len(seq) for seq in sequences.values())
    for seq_id in sequences:
        sequences[seq_id] = sequences[seq_id].ljust(max_length)

    return sequences

def create_report(data):
    env = Environment(loader=FileSystemLoader(''))
    template = env.get_template('report.html')
    html_content = template.render(data)

    with open('results/report_generated.html', 'w', encoding='utf-8') as html_file:
        html_file.write(html_content)