"""
scFv Design Space — Hugging Face Gradio app
Two-tab interface:
  Tab 1: scFv Builder  — PDB/CIF or FASTA → scFv FASTA for AlphaFold Server
  Tab 2: scFv-MPNN-Light — CIF from AlphaFold → ranked redesign + scoring
"""

import gradio as gr
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Generator, Optional, Tuple

SCRIPTS_DIR = Path(__file__).resolve().parent / "scripts"
PMPNN_DIR   = Path(__file__).resolve().parent / "ProteinMPNN"


# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────

def _normalize_epitope_chain(value) -> str:
    if value is None:
        return "None"
    if isinstance(value, (list, tuple)):
        items = [str(v).strip() for v in value if str(v).strip()]
        items = [v for v in items if v.lower() != "none"]
        return ",".join(items) if items else "None"
    s = str(value).strip()
    if not s or s.lower() == "none":
        return "None"
    parts = re.split(r"[\s,;]+", s)
    parts = [p.strip() for p in parts if p.strip() and p.lower() != "none"]
    return ",".join(parts) if parts else "None"


def _stream(log: list, line: str = "", zip_path=None, html=None):
    if line:
        log.append(line)
    return "".join(log), str(zip_path) if zip_path else None, html


# ─────────────────────────────────────────────────────────────────────────────
# Tab 1 — scFv Builder
# ─────────────────────────────────────────────────────────────────────────────

def run_builder(
    input_file,
    determine_cdrs: str,
    linker_seq: str,
    epitope_chain: str,
) -> Generator[Tuple[str, Optional[str], Optional[str]], None, None]:
    """
    Run make_scfv_vernierV3.py on the uploaded file.
    Outputs scfv_output.fasta (and optional multimer FASTA if epitope supplied).
    """
    workdir = Path(tempfile.mkdtemp())
    input_dir = workdir / "input"
    input_dir.mkdir()
    script_dst = workdir / "scripts"
    shutil.copytree(SCRIPTS_DIR, script_dst)

    infile = input_dir / Path(input_file.name).name
    shutil.copy(input_file.name, infile)

    ext = infile.suffix.lower()
    log = []

    yield _stream(log, f"Working directory: {workdir}\n")
    yield _stream(log, f"Input: {infile.name}\n")

    # CIF → PDB if needed
    if ext in (".cif", ".mmcif"):
        pdb_path = infile.with_suffix(".pdb")
        yield _stream(log, "Converting CIF → PDB...\n")
        result = subprocess.run(
            ["python", "-c",
             f"import gemmi, pathlib; "
             f"doc = gemmi.cif.read_file('{infile}'); "
             f"st = gemmi.make_structure_from_block(doc.sole_block()); "
             f"st.write_pdb('{pdb_path}')"],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            yield _stream(log, f"ERROR converting CIF:\n{result.stderr}\n")
            return
        infile = pdb_path
        yield _stream(log, f"Converted to: {infile.name}\n")

    # Run make_scfv_vernierV3.py
    cmd = [
        "python", str(script_dst / "make_scfv_vernierV3.py"),
        "--scheme", determine_cdrs,
        "--combine",
        "--linker", linker_seq,
        "--save-anarci",
        "--simple-output", "design",
        str(infile),
    ]
    yield _stream(log, f"\nRunning scFv builder...\n")

    proc = subprocess.Popen(
        cmd, cwd=workdir, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, text=True, bufsize=1
    )
    if proc.stdout:
        for line in proc.stdout:
            log.append(line)
            yield _stream(log)
    retcode = proc.wait()

    if retcode != 0:
        yield _stream(log, "\nERROR: scFv builder failed.\n")
        return

    scfv_fasta = workdir / "scfv_output.fasta"
    if not scfv_fasta.exists():
        yield _stream(log, "\nERROR: scfv_output.fasta not found.\n")
        return

    # If epitope chain provided, append epitope sequence(s) to the FASTA
    output_fasta = scfv_fasta
    epi = _normalize_epitope_chain(epitope_chain)
    if epi != "None" and infile.suffix.lower() in (".pdb",):
        yield _stream(log, f"\nExtracting epitope chain(s): {epi}\n")
        epi_seqs = []
        for chain_id in epi.split(","):
            result = subprocess.run(
                ["python", str(script_dst / "extract_chain_seq.py"), str(infile), chain_id],
                capture_output=True, text=True
            )
            if result.returncode != 0 or not result.stdout.strip():
                yield _stream(log, f"WARNING: Could not extract chain {chain_id}\n")
            else:
                epi_seqs.append((chain_id, result.stdout.strip()))

        if epi_seqs:
            multimer_fasta = workdir / "scfv_alphafold_input.fasta"
            with open(scfv_fasta) as f:
                scfv_text = f.read().strip()
            with open(multimer_fasta, "w") as f:
                f.write(scfv_text + "\n")
                for chain_id, seq in epi_seqs:
                    f.write(f"\n>epitope_chain_{chain_id}\n{seq}\n")
            output_fasta = multimer_fasta
            yield _stream(log, f"Added {len(epi_seqs)} epitope chain(s) to FASTA.\n")

    fasta_text = output_fasta.read_text()
    yield _stream(log, f"\nDone. Download your scFv FASTA below.\n",
                  zip_path=None, html=None)
    yield "".join(log), str(output_fasta), fasta_text


# ─────────────────────────────────────────────────────────────────────────────
# Tab 2 — scFv-MPNN-Light
# ─────────────────────────────────────────────────────────────────────────────

def _log_run(params: dict):
    """Append an anonymous run record to runs_log.csv in the HF dataset."""
    token = os.environ.get("HF_TOKEN")
    if not token:
        return
    try:
        import datetime, csv, io
        from huggingface_hub import HfApi
        api = HfApi(token=token)
        dataset_repo = "timstasevich/scfv-summaries"
        api.create_repo(repo_id=dataset_repo, repo_type="dataset",
                        private=True, exist_ok=True)

        # Try to download existing log
        existing = ""
        try:
            from huggingface_hub import hf_hub_download
            path = hf_hub_download(repo_id=dataset_repo, filename="runs_log.csv",
                                   repo_type="dataset", token=token)
            with open(path) as f:
                existing = f.read()
        except Exception:
            pass

        # Append new row
        fields = ["timestamp", "tab", "scheme", "seqs_per_run", "epitope",
                  "model", "sampling_temp", "backbone_noise", "sort_column"]
        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=fields, extrasaction="ignore")
        if not existing:
            writer.writeheader()
        else:
            buf.write(existing.rstrip("\n") + "\n")
            writer = csv.DictWriter(buf, fieldnames=fields, extrasaction="ignore")
        params["timestamp"] = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
        writer.writerow(params)

        api.upload_file(
            path_or_fileobj=buf.getvalue().encode(),
            path_in_repo="runs_log.csv",
            repo_id=dataset_repo,
            repo_type="dataset",
        )
    except Exception:
        pass  # never block the pipeline over logging


def _upload_to_hf_dataset(html_path: Path, log: list):
    """Upload summary HTML to HF dataset if HF_TOKEN is set."""
    token = os.environ.get("HF_TOKEN")
    if not token:
        log.append("[share] HF_TOKEN not set — skipping upload.\n")
        return
    try:
        from huggingface_hub import HfApi
        import datetime
        api = HfApi(token=token)
        dataset_repo = "timstasevich/scfv-summaries"
        # Create dataset repo if it doesn't exist
        api.create_repo(repo_id=dataset_repo, repo_type="dataset",
                        private=True, exist_ok=True)
        timestamp = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
        dest = f"summaries/{timestamp}_{html_path.name}"
        api.upload_file(
            path_or_fileobj=str(html_path),
            path_in_repo=dest,
            repo_id=dataset_repo,
            repo_type="dataset",
        )
        log.append(f"[share] Uploaded summary → {dest}\n")
    except Exception as e:
        log.append(f"[share] Upload failed: {e}\n")


def run_mpnn_light(
    input_file,
    determine_cdrs: str,
    ss_near_cdrs: int,
    linker_seq: str,
    seqs_per_run: int,
    epitope_chain: str,
    pmpnn_seed: int,
    pmpnn_model: str,
    sampling_temp: float,
    backbone_noise: float,
    sort_column: str,
    Nt: int,
    Nb: int,
    Nr: int,
    share_results: bool,
) -> Generator[Tuple[str, Optional[str], Optional[str]], None, None]:

    workdir   = Path(tempfile.mkdtemp())
    input_dir = workdir / "input"
    output_dir = workdir / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    script_dst = workdir / "scripts"
    shutil.copytree(SCRIPTS_DIR, script_dst)

    infile = input_dir / Path(input_file.name).name
    shutil.copy(input_file.name, infile)

    log = []
    yield _stream(log, f"Working directory: {workdir}\n")

    env = os.environ.copy()
    env["input_dir"]       = str(input_dir)
    env["output_dir"]      = str(output_dir)
    env["determine_CDRs"]  = determine_cdrs
    env["ss_near_CDRs"]    = str(ss_near_cdrs)
    env["linker_seq"]      = linker_seq
    env["seqs_per_run"]    = str(seqs_per_run)
    env["epitope_chain"]   = _normalize_epitope_chain(epitope_chain)
    env["PMPNN"]           = str(PMPNN_DIR)
    env["pmpnn_seed"]        = str(pmpnn_seed)
    env["sampling_temp"]     = str(sampling_temp)
    env["backbone_noise"]    = str(backbone_noise)
    env["use_soluble_model"] = "true" if pmpnn_model == "soluble" else "false"
    env["use_gpu"]           = "false"
    env["sort_column"]     = sort_column
    env["Nt"]              = str(Nt)
    env["Nb"]              = str(Nb)
    env["Nr"]              = str(Nr)
    env["CUDA_VISIBLE_DEVICES"] = ""

    pipeline_sh = Path(__file__).resolve().parent / "run_light.sh"

    proc = subprocess.Popen(
        ["bash", str(pipeline_sh)],
        cwd=workdir, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1
    )
    if proc.stdout:
        for line in proc.stdout:
            log.append(line)
            yield _stream(log)
    retcode = proc.wait()

    log.append(f"\nPipeline exit code: {retcode}\n")

    if retcode != 0:
        yield _stream(log, "\nERROR: Pipeline failed.\n")
        return

    # Zip results
    zip_base = workdir / "results"
    shutil.make_archive(str(zip_base), "zip", output_dir)
    zip_path = zip_base.with_suffix(".zip")

    # Load HTML preview
    html_preview = None
    summary_dir = output_dir / "Summary"
    html_file = None
    if summary_dir.exists():
        html_files = sorted(summary_dir.glob("*_summary.html"))
        if html_files:
            html_file = html_files[0]
            try:
                html_text = html_file.read_text()
                html_preview = (
                    "<div style='max-height:600px;overflow:auto;"
                    "border:1px solid #ccc;padding:8px;background:white;'>"
                    + html_text + "</div>"
                )
            except Exception as e:
                log.append(f"\nWARNING: Could not load HTML preview: {e}\n")

    # Anonymous run log + opt-in upload — run in background so they don't block the response
    import threading
    def _background():
        _log_run({
            "tab": "step2",
            "scheme": determine_cdrs,
            "seqs_per_run": seqs_per_run,
            "epitope": _normalize_epitope_chain(epitope_chain),
            "model": pmpnn_model,
            "sampling_temp": sampling_temp,
            "backbone_noise": backbone_noise,
            "sort_column": sort_column,
        })
        if share_results and html_file is not None:
            _upload_to_hf_dataset(html_file, [])
    threading.Thread(target=_background, daemon=True).start()

    yield _stream(log, zip_path=zip_path, html=html_preview)


# ─────────────────────────────────────────────────────────────────────────────
# Gradio UI
# ─────────────────────────────────────────────────────────────────────────────

WORKFLOW_MD = """
## Two-step antibody to scFv intrabody design pipeline

### How it works

**Step 1 — scFv Builder:** Upload an antibody structure (PDB/CIF) or sequence (FASTA). The app identifies VH and VL domains, assembles them into a single-chain variable fragment (scFv) with your chosen linker, and outputs a FASTA ready for structure prediction.

**→ Fold the scFv:** Use any tool to obtain a 3D structure from the FASTA, e.g. [AlphaFold Server](https://alphafoldserver.com). Download the resulting CIF or PDB file. If you have an epitope sequence, we recommend folding the scFv together with the epitope as a multimer — this gives ProteinMPNN a more realistic binding context.

**Step 2 — scFv-MPNN-Light:** Upload the folded CIF (or any existing scFv PDB/CIF). ProteinMPNN redesigns the framework regions while keeping CDRs fixed. Designs are scored with scfvtools and ranked for download.

**If you already have a folded scFv structure, you can skip Step 1 and go directly to Step 2.**

---

### Citation

If you use this tool, please cite:

> Galindo et al., *AI-assisted protein design to rapidly convert antibody sequences to intrabodies targeting diverse peptides and histone modifications.* Science Advances (2026). [doi:10.1126/sciadv.adx8352](https://www.science.org/doi/10.1126/sciadv.adx8352)

Source code: [github.com/Colorado-State-University-Stasevich-Lab/scfv-design-space](https://github.com/Colorado-State-University-Stasevich-Lab/scfv-design-space)

This app also uses:
> Dauparas et al., *Robust deep learning–based protein sequence design using ProteinMPNN.* Science 378, 49–56 (2022).

> Dunbar & Deane, *ANARCI: antigen receptor numbering and receptor classification.* Bioinformatics 32, 298–300 (2016).
"""

with gr.Blocks(title="Antibody to Intrabody Pipeline") as demo:
    gr.Markdown("# Antibody to Intrabody Pipeline")
    gr.Markdown(WORKFLOW_MD)

    with gr.Tabs():

        # ── Tab 1 ──────────────────────────────────────────────────────────
        with gr.Tab("Step 1 — scFv Builder"):
            gr.Markdown(
                "Upload an antibody **PDB/CIF** (or **FASTA** with VH+VL) "
                "to generate a scFv sequence ready for AlphaFold Server."
            )
            with gr.Row():
                with gr.Column(scale=1):
                    t1_file = gr.File(
                        label="Upload structure (PDB or CIF) or sequence (FASTA)",
                        file_types=[".pdb", ".cif", ".mmcif", ".fasta", ".fa"],
                    )
                    t1_scheme = gr.Dropdown(
                        choices=["martin", "chothia", "kabat", "imgt"],
                        value="martin", label="CDR numbering scheme",
                        visible=False,
                    )
                    t1_linker = gr.Textbox(
                        value="GGGGSGGGGSGGGGS", label="Linker sequence",
                    )
                    t1_epitope = gr.Textbox(
                        value="None",
                        label="Epitope chain ID(s) — optional",
                        info="e.g. 'B' or 'B,C'. If provided, a multimer FASTA is output.",
                    )
                    t1_btn = gr.Button("Build scFv FASTA", variant="primary")

                with gr.Column(scale=1):
                    t1_log = gr.Textbox(label="Log", lines=20, interactive=False)
                    t1_file_out = gr.File(label="Download scFv FASTA")
                    t1_fasta_preview = gr.Textbox(
                        label="FASTA preview", lines=6, interactive=False,
                    )

            t1_btn.click(
                fn=run_builder,
                inputs=[t1_file, t1_scheme, t1_linker, t1_epitope],
                outputs=[t1_log, t1_file_out, t1_fasta_preview],
            )

        # ── Tab 2 ──────────────────────────────────────────────────────────
        with gr.Tab("Step 2 — scFv-MPNN-Light"):
            gr.Markdown(
                "Upload a **CIF or PDB** of a folded scFv to redesign and score "
                "framework sequences using ProteinMPNN + scfvtools. "
                "If you already have a scFv structure, you can skip Step 1 and use it directly here."
            )
            with gr.Row():
                with gr.Column(scale=1):
                    t2_file = gr.File(
                        label="Upload structure (CIF or PDB)",
                        file_types=[".cif", ".mmcif", ".pdb"],
                    )
                    t2_scheme = gr.Dropdown(
                        choices=["martin", "chothia", "kabat", "imgt"],
                        value="martin", label="CDR numbering scheme",
                    )
                    t2_dist = gr.Slider(
                        0, 10, step=1, value=3,
                        label="Fixed residue distance near CDRs (Å)",
                    )
                    t2_linker = gr.Textbox(
                        value="GGGGSGGGGSGGGGS", label="Linker sequence",
                    )
                    t2_epitope = gr.Textbox(
                        value="None",
                        label="Epitope chain ID(s) — optional",
                        info="The pipeline assumes the scFv is chain A. If your structure contains a bound antigen, you can enter its chain ID(s) here (e.g. 'B' or 'B,C') to include it in the multimer FASTA output. Epitope chains are not required for redesign or scoring.",
                    )
                    t2_seqs = gr.Slider(
                        1, 500, step=1, value=30,
                        label="Sequences per run (ProteinMPNN)",
                        info="Default 30 runs in ~60–90s on CPU. Values above 50 may timeout on the free tier — run locally for large jobs.",
                    )
                    t2_seed = gr.Number(value=37, precision=0, label="ProteinMPNN seed")
                    t2_model = gr.Dropdown(
                        choices=["soluble", "vanilla"],
                        value="soluble",
                        label="ProteinMPNN model",
                        info="Soluble model recommended for intrabodies; vanilla for general use.",
                    )
                    t2_temp = gr.Slider(
                        0.0, 1.0, step=0.05, value=0.1,
                        label="Sampling temperature",
                        info="Higher = more diversity, lower = closer to input sequence.",
                    )
                    t2_noise = gr.Slider(
                        0.0, 0.2, step=0.01, value=0.02,
                        label="Backbone noise (Å)",
                        info="Gaussian noise added to backbone coordinates during design.",
                    )
                    t2_sort = gr.Dropdown(
                        choices=[
                            "scfvtools_blosum_diff_score",
                            "scfvtools_blosum_score",
                            "scfvtools_score",
                        ],
                        value="scfvtools_blosum_diff_score",
                        label="Sort / rank by",
                    )
                    with gr.Row():
                        t2_Nt = gr.Slider(1, 50, step=1, value=10, label="Top N")
                        t2_Nb = gr.Slider(1, 50, step=1, value=10, label="Bottom N")
                        t2_Nr = gr.Slider(1, 50, step=1, value=10, label="Random N")
                    t2_share = gr.Checkbox(
                        value=False,
                        label="Share summary with developers",
                        info="If checked, your HTML summary will be anonymously uploaded to help us improve the tool.",
                    )
                    t2_btn = gr.Button("Run scFv-MPNN-Light", variant="primary")

                with gr.Column(scale=1):
                    t2_log  = gr.Textbox(label="Pipeline log (live)", lines=25, interactive=False)
                    t2_zip  = gr.File(label="Download results.zip — the best way to explore results is to open Summary/*_summary.html from this zip in your browser")
                    t2_html = gr.HTML(label="Summary preview")

            t2_btn.click(
                fn=run_mpnn_light,
                inputs=[
                    t2_file, t2_scheme, t2_dist, t2_linker,
                    t2_seqs, t2_epitope, t2_seed,
                    t2_model, t2_temp, t2_noise,
                    t2_sort, t2_Nt, t2_Nb, t2_Nr, t2_share,
                ],
                outputs=[t2_log, t2_zip, t2_html],
            )

    gr.Markdown("""
---

### CDR definitions

CDR positions are determined by [ANARCI](https://github.com/oxpig/ANARCI) (Dunbar & Deane, *Bioinformatics* 2016).
The following numbering schemes are supported:

| Scheme | VH CDR1 | VH CDR2 | VH CDR3 | VL CDR1 | VL CDR2 | VL CDR3 |
|--------|---------|---------|---------|---------|---------|---------|
| **Martin** (default) | 30–35 | 47–58 | 95–101 | 30–35 | 46–55 | 89–97 |
| **Chothia** | 26–32 | 52–56 | 95–102 | 26–32 | 50–52 | 91–96 |
| **Kabat** | 30–35 | 50–65 | 95–102 | 24–34 | 50–56 | 89–97 |
| **IMGT** | 27–38 | 56–65 | 105–117 | 27–38 | 56–65 | 105–117 |

Residues within the specified distance (Å) of any CDR atom are also held fixed during ProteinMPNN redesign.
""")

if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0", server_port=7860,
                allowed_paths=[str(Path(__file__).resolve().parent)])
