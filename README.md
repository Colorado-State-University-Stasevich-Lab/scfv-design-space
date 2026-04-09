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

# Antibody to Intrabody Pipeline

A two-step pipeline for converting antibody sequences into scFv intrabodies using ProteinMPNN.

## Overview

**Step 1 — scFv Builder:** Upload an antibody structure (PDB/CIF) or sequence (FASTA). The app identifies VH and VL domains, assembles them into a single-chain variable fragment (scFv) with a linker, and outputs a FASTA ready for structure prediction.

**→ Fold the scFv:** Use any tool to obtain a 3D structure from the FASTA, e.g. the [AlphaFold server](https://alphafoldserver.com). If you have an epitope, we recommend folding the scFv together with the epitope as a multimer.

**Step 2 — scFv-MPNN-Light:** Upload the folded scFv structure (PDB or CIF). ProteinMPNN redesigns the framework while keeping CDRs and nearby residues fixed. Designs are ranked by scfvtools, which scores each sequence by comparison to a consensus derived from intrabodies proven to work.

**If you already have a folded scFv structure, skip Step 1 and go directly to Step 2.**

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
