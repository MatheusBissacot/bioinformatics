# Phylogenetics of a Taste 2 Receptor Member Gene (T2R38)

2024/2025 2ºSemestre 1ºAno MIA

This repository contains the scripts and report for the phylogenetic analysis of the **T2R38** gene, developed as part of the Bioinformatics course at the University of Porto, Master in Artifical Inteligence.

## 📄 Description

The main goal of this project is to analyze the evolution of the T2R38 gene, which encodes a bitter taste receptor in humans. The workflow includes:

- Retrieving the reference DNA sequence (FASTA) from NCBI
- Translating the DNA sequence into a protein
- Searching for homologous sequences using BLASTp
- Performing multiple sequence alignment (MSA)
- Building a phylogenetic tree
- Analyzing and interpreting the evolutionary relationships

## 🧪 Dependencies

This project was developed using Python 3 and requires the following libraries:

- `biopython`
- `matplotlib`
- `ete3`

You can install the dependencies via:

```bash
pip install -r requirements.txt
```

## How to run 

```bash
python phyloanalysis.py sequence.fasta number_of_differente_species
```

## Report
A complete LaTeX report with figures, phylogenetic trees, and conclusions is available in the report/ folder.
Also inclue a HTML file "report_generated" which contains a web-style report for the basic analysis carried out in this project.

## Authors

Gonçalo Brochado (Faculty of Engineering, University of Porto)
Matheus Bissacot (Faculty of Engineering, University of Porto)
