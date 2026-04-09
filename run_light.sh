#!/bin/bash
# run_light.sh — scFv-MPNN-Light pipeline
# Requires a PDB or CIF structure as input (FASTA not supported).
# No ColabFold. No struct-evo. Single conda environment: scfv_design.
#
# Steps:
#   1. CDR determination + scFv FASTA  (ANARCI / make_scfv_vernierV3.py)
#   2. ProteinMPNN design on input structure
#   3. Build combined FASTA (WT + designs)
#   4. Score: scfvtools + SWI
#   5. Merge scores, rank, select Top / Bottom / Random
#   6. Build FASTA subsets + HTML summary

###############################################################
### USER CONFIGURATION (OVERRIDABLE VIA ENV VARS)
###############################################################

input_dir="${input_dir:-input_dir}"
output_dir="${output_dir:-output_dir}"
determine_CDRs="${determine_CDRs:-martin}"
ss_near_CDRs="${ss_near_CDRs:-3}"
linker_seq="${linker_seq:-GGGGSGGGGSGGGGS}"
seqs_per_run="${seqs_per_run:-5}"
epitope_chain="${epitope_chain:-None}"
PMPNN="${PMPNN:-/path/to/ProteinMPNN}"   # clone from https://github.com/dauparas/ProteinMPNN
pmpnn_seed="${pmpnn_seed:-37}"
sampling_temp="${sampling_temp:-0.1}"
backbone_noise="${backbone_noise:-0.02}"
use_soluble_model="${use_soluble_model:-true}"  # false = vanilla model
use_gpu="${use_gpu:-false}"              # set true only if GPU environment installed

# scoring / selection controls
Nt="${Nt:-10}"
Nb="${Nb:-10}"
Nr="${Nr:-10}"
sort_column="${sort_column:-scfvtools_blosum_diff_score}"

###############################################################
### Helpers
###############################################################

parse_epitope_chains () {
    local raw="$1"
    raw="${raw//;/,}"
    raw="${raw// /,}"
    raw=$(echo "$raw" | tr -s ',')
    IFS=',' read -r -a _chains <<< "$raw"
    epitope_chains=()
    for c in "${_chains[@]}"; do
        c="$(echo "$c" | tr -d '[:space:]')"
        if [ -z "$c" ]; then continue; fi
        if [ "$c" = "None" ] || [ "$c" = "none" ]; then continue; fi
        epitope_chains+=("$c")
    done
}

join_by_colon () {
    local IFS=":"
    echo "$*"
}

###############################################################

# Activate single environment
if [ "$CONDA_DEFAULT_ENV" != "scfv_design" ]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate scfv_design
fi

# GPU / CPU control for ProteinMPNN
if [ "$use_gpu" != "true" ]; then
    export CUDA_VISIBLE_DEVICES=""
fi

mkdir -p "$output_dir"
summary_dir="$output_dir/Summary"
mkdir -p "$summary_dir"

###############################################################
### STEP 0: Load input file
###############################################################

input_file=$(find "$input_dir" -maxdepth 1 -type f | head -n 1)
if [ -z "$input_file" ]; then
    echo "ERROR: No input file found in $input_dir"
    exit 1
fi

base_name=$(basename "$input_file")
base_noext="${base_name%.*}"
ext="${input_file##*.}"
ext="${ext,,}"

echo "Input file = $input_file"

# Light mode requires a structure — FASTA is not supported
if [[ "$ext" == "fa" || "$ext" == "fasta" ]]; then
    echo "ERROR: FASTA input is not supported in light mode."
    echo "       Provide a PDB or CIF file."
    exit 1
fi

# CIF → PDB conversion
if [[ "$ext" == "cif" || "$ext" == "mmcif" ]]; then
    echo "Detected CIF input. Converting to PDB..."
    python - << EOF
import gemmi, pathlib
cif_path = pathlib.Path("$input_file")
pdb_path = cif_path.with_suffix(".pdb")
doc = gemmi.cif.read_file(str(cif_path))
structure = gemmi.make_structure_from_block(doc.sole_block())
structure.write_pdb(str(pdb_path))
EOF
    input_file="${input_file%.*}.pdb"
    base_name=$(basename "$input_file")
    base_noext="${base_name%.*}"
    ext="pdb"
    echo "Converted PDB written to: $input_file"
fi

###############################################################
### STEP 1: Determine CDRs + scFv FASTA
###############################################################

parse_epitope_chains "$epitope_chain"

design_positions="$(python scripts/make_scfv_vernierV3.py \
    --scheme "$determine_CDRs" \
    --combine \
    --linker "$linker_seq" \
    --save-anarci \
    --dist "$ss_near_CDRs" \
    --simple-output design \
    "$input_file")"

if [ $? -ne 0 ]; then
    echo "ERROR: make_scfv_vernierV3.py failed."
    exit 1
fi

if [ ! -f "scfv_output.fasta" ]; then
    echo "ERROR: scfv_output.fasta was not produced."
    exit 1
fi

echo "Design positions: $design_positions"
echo "$design_positions" > "$summary_dir/${base_noext}_design_positions.txt"

for f in \
    "${base_noext}_summary.html" \
    "${base_noext}_fixed_positions.txt" \
    "${base_noext}_anarci.txt" \
    VH_mapping.tsv \
    VL_mapping.tsv \
    pdb_linear_map.tsv; do
    [ -f "$f" ] && mv "$f" "$summary_dir/"
done

mv scfv_output.fasta "$output_dir/"
combined_fasta="$output_dir/scfv_output.fasta"
scfv_seq=$(awk '/^>scFv/{getline; print; exit}' "$combined_fasta")

###############################################################
### STEP 2: Extract epitope sequence (structure input only)
###############################################################

epitope_seq=""
if [ "${#epitope_chains[@]}" -gt 0 ]; then
    epitope_seqs=()
    for c in "${epitope_chains[@]}"; do
        seq=$(python scripts/extract_chain_seq.py "$input_file" "$c")
        if [ -z "$seq" ]; then
            echo "ERROR: Could not extract sequence for epitope chain '$c' from $input_file"
            exit 1
        fi
        epitope_seqs+=("$seq")
    done
    epitope_seq=$(join_by_colon "${epitope_seqs[@]}")
fi

rank_file="$input_file"
echo "Using input structure for MPNN: $rank_file"

###############################################################
### STEP 3: Run ProteinMPNN
###############################################################

echo "Running ProteinMPNN..."

dir="$output_dir/MPNN_input"
mkdir -p "$dir"
cp "$rank_file" "$dir/"

parsed_jsonl="$dir/parsed_pdbs.jsonl"
assigned_jsonl="$dir/assigned_pdbs.jsonl"
fixed_jsonl="$dir/fixed_pdbs.jsonl"
chains_to_design="A"

python "$PMPNN/helper_scripts/parse_multiple_chains.py" \
    --input_path "$dir" \
    --output_path "$parsed_jsonl"

python "$PMPNN/helper_scripts/assign_fixed_chains.py" \
    --input_path "$parsed_jsonl" \
    --output_path "$assigned_jsonl" \
    --chain_list "$chains_to_design"

python "$PMPNN/helper_scripts/make_fixed_positions_dict.py" \
    --input_path "$parsed_jsonl" \
    --output_path "$fixed_jsonl" \
    --chain_list "$chains_to_design" \
    --position_list "$design_positions" \
    --specify_non_fixed

pmpnn_cmd=(
    python "$PMPNN/protein_mpnn_run.py"
    --jsonl_path "$parsed_jsonl"
    --chain_id_jsonl "$assigned_jsonl"
    --fixed_positions_jsonl "$fixed_jsonl"
    --out_folder "$output_dir"
    --num_seq_per_target "$seqs_per_run"
    --sampling_temp "$sampling_temp"
    --seed "$pmpnn_seed"
    --batch_size 1
    --backbone_noise "$backbone_noise"
)
if [ "$use_soluble_model" = "true" ]; then
    pmpnn_cmd+=(--use_soluble_model)
fi
"${pmpnn_cmd[@]}"

if [ $? -ne 0 ]; then
    echo "ERROR: ProteinMPNN failed."
    exit 1
fi

if [ ! -d "$output_dir/seqs" ] || ! ls "$output_dir/seqs"/*.fa >/dev/null 2>&1; then
    echo "ERROR: ProteinMPNN produced no FASTA outputs."
    exit 1
fi

echo "ProteinMPNN completed successfully."

###############################################################
### STEP 4: Build combined FASTA (WT + designs)
###############################################################

python scripts/make_combined_multimer_fasta.py \
    "$output_dir/seqs" \
    "$epitope_seq" \
    "$output_dir/combined_multimer.fa"

summary_html="$summary_dir/${base_noext}_summary.html"
combined_fa="$output_dir/combined_multimer.fa"

python scripts/append_designs_to_summary.py \
    "$summary_html" \
    "$combined_fa"

###############################################################
### STEP 5: Score designs
###############################################################

clean_fa="$output_dir/combined_multimer_chainA_only.fa"
scfvtools_csv="$output_dir/scfvtools_scores.csv"
merged_csv="$output_dir/merged_scores.csv"

echo "Cleaning FASTA (chain A only)..."
awk '
    /^>/ { print; next }
    { split($0, arr, ":"); print arr[1] }
' "$combined_fa" > "$clean_fa"

########################################
# 5A. scfvtools scoring
########################################

reference_csv=$(python -c \
    "import scfvtools, os; print(os.path.join(os.path.dirname(scfvtools.__file__), 'example_data', 'diff_H.csv'))")

echo "Running scfvtools scoring..."
run_scfvtools_scoring.py \
    --fasta "$clean_fa" \
    --scheme martin \
    --reference_csv "$reference_csv" \
    --reference_name "diff_H" \
    --out_csv "$scfvtools_csv" \
    --summary_html "$summary_html"

echo "scfvtools scores written to: $scfvtools_csv"

###############################################################
### STEP 5B: Merge + select Top / Bottom / Random
###############################################################

sort_column="${sort_column:-scfvtools_blosum_diff_score}"

python scripts/merge_preAF2_scores.py \
    --scfv "$scfvtools_csv" \
    --sort_by "$sort_column" \
    --out "$merged_csv"

echo "Selecting Top ($Nt), Bottom ($Nb), Random ($Nr)..."
python scripts/select_top_bottom_random.py \
    --input "$merged_csv" \
    --sortcol "$sort_column" \
    --Nt "$Nt" \
    --Nb "$Nb" \
    --Nr "$Nr" \
    --outdir "$summary_dir"

###############################################################
### STEP 5D: Build FASTA subsets + update HTML summary
###############################################################

top_fa="$summary_dir/top_designs.fa"
bottom_fa="$summary_dir/bottom_designs.fa"
random_fa="$summary_dir/random_designs.fa"

extract_chainA () {
    local acc_file="$1"
    local out_fa="$2"
    : > "$out_fa"

    # Always include WT (first record)
    awk 'BEGIN{RS=">"; FS="\n"} NR==2 { print ">"$1; print $2 }' \
        "$combined_fa" >> "$out_fa"

    while IFS=, read -r accession rest; do
        [[ "$accession" =~ ^Accession ]] && continue
        design_num=$(echo "$accession" | grep -oE 'design_([0-9]+)' | grep -oE '[0-9]+')
        [ -z "$design_num" ] && continue
        awk -v num="$design_num" \
            'BEGIN{RS=">"; FS="\n"} $1 ~ ("_design_" num "$") { print ">"$1; print $2 }' \
            "$combined_fa" >> "$out_fa"
    done < "$acc_file"
}

extract_chainA "$summary_dir/selected_top.csv"    "$top_fa"
extract_chainA "$summary_dir/selected_bottom.csv" "$bottom_fa"
extract_chainA "$summary_dir/selected_random.csv" "$random_fa"

python scripts/append_designs_to_summary.py "$summary_html" "$top_fa"
python scripts/append_designs_to_summary.py "$summary_html" "$bottom_fa"
python scripts/append_designs_to_summary.py "$summary_html" "$random_fa"

###############################################################
### STEP 6: Append score table to HTML summary
###############################################################

python scripts/append_scores_to_summary.py \
    --csv "$merged_csv" \
    --summary "$summary_html" \
    --top "$summary_dir/selected_top.csv" \
    --bottom "$summary_dir/selected_bottom.csv" \
    --random "$summary_dir/selected_random.csv"

echo ""
echo "Done. Summary: $summary_html"
