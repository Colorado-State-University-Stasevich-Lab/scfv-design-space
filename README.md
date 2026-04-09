---
title: Scfv Design Space
emoji: 🦀
colorFrom: blue
colorTo: indigo
sdk: docker
pinned: false
license: mit
short_description: Two-step scFv antibody redesign with ProteinMPNN
---

# scFv Design Space

A two-step pipeline for scFv antibody redesign using ProteinMPNN and structure-based CDR analysis.

## Overview

**Step 1 — scFv Builder:** Provide an antibody heavy chain (VH) and light chain (VL) sequence. The app runs [ColabFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) (via AlphaFold2) to fold the scFv, identifies CDR loops using ANARCI, and produces a summary HTML with numbered sequences and a FASTA ready for Step 2.

**Step 2 — scFv-MPNN-Light:** Provide a folded scFv structure (PDB or CIF — fold using any method, e.g. the [AlphaFold server](https://alphafoldserver.com)). ProteinMPNN redesigns CDR and Vernier positions while holding the framework fixed. Designs are scored with scfvtools and ranked by BLOSUM substitution score. Top, bottom, and random selections are returned as a downloadable HTML summary.

If you already have a folded scFv structure, skip Step 1 and go directly to Step 2.

## Features

- CDR determination via ANARCI (Martin, Chothia, Kabat, or IMGT numbering)
- ProteinMPNN soluble or vanilla model
- Configurable sampling temperature, backbone noise, and number of designs
- scfvtools BLOSUM-based scoring relative to a reference antibody distribution
- Summary HTML with ANARCI-numbered sequences, score table, and FASTA subsets

## Citation

If you use this tool, please cite:

> Galindo et al., *AI-assisted protein design to rapidly convert antibody sequences to intrabodies targeting diverse peptides and histone modifications.* Science Advances (2026). [doi:10.1126/sciadv.adx8352](https://www.science.org/doi/10.1126/sciadv.adx8352)

This app also uses:

> Dauparas et al., *Robust deep learning–based protein sequence design using ProteinMPNN.* Science 378, 49–56 (2022).

> Dunbar & Deane, *ANARCI: antigen receptor numbering and receptor classification.* Bioinformatics 32, 298–300 (2016).

## Source code

- HuggingFace Space: [timstasevich/scfv-design-space](https://huggingface.co/spaces/timstasevich/scfv-design-space)
- Local pipeline (WSL/Linux): [Colorado-State-University-Stasevich-Lab/scfv-mpnn-light](https://github.com/Colorado-State-University-Stasevich-Lab/scfv-mpnn-light)
- Scoring library: [Colorado-State-University-Stasevich-Lab/scfvtools](https://github.com/Colorado-State-University-Stasevich-Lab/scfvtools)
