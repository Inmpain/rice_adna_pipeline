#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

usage() {
  cat <<'EOF'
Compare Bowtie2 and BWA for extracting Asian rice/Oryza-related reads.

Modes:
  1. Local mode (default): run all samples in the current shell
  2. Sbatch mode: submit per-sample Bowtie2 and BWA jobs separately
  3. Collect mode: merge finished per-sample stats into summary tables

Usage:
  bash compare_rice_read_extractors.sh [options]

Common options:
  --read-glob GLOB         Repeatable. FASTQ/FQ glob(s) to scan.
  --task-list FILE         Read FASTQ tasks from a file, one FASTQ path per line.
  --ref-fasta PATH         Reference FASTA used by BWA.
  --bt2-prefix PATH        Bowtie2 index prefix. Default: same as --ref-fasta.
  --out-dir PATH           Output directory.
  --threads N              Threads per mapper. Default: 20.
  --count-input-reads      Also count input FASTQ reads. Slower on very large files.
  --force                  Re-run mapping even if output files already exist.
  --no-modules             Do not try to load environment modules.

Sbatch options:
  --sbatch                 Submit jobs to Slurm.
  --sbatch-partition NAME  Slurm partition. Default: use cluster default partition.
  --sbatch-nodelist LIST   Slurm nodelist, e.g. node01,node02,node03. This is a hard constraint.
  --sbatch-time TIME       Slurm time for Bowtie2/BWA jobs. Empty means do not pass --time.
  --sbatch-mem MEM         Slurm memory for Bowtie2/BWA jobs. Empty means do not pass --mem.
  --no-collector-job       Do not submit the final summary collector job.

Summary options:
  --collect-only           Merge finished per-sample stats into summary tables.
  --write-pending-list F   Write unfinished FASTQ tasks to file F, then exit.
  --split-pending N        Split unfinished FASTQ tasks into N task-list files, then exit.
  --split-prefix PREFIX    Prefix used with --split-pending. Output files: PREFIX.1.txt, etc.

Other:
  -h, --help               Show this help.

Default input globs:
  /home/scratch/yinmt202607/3.angkor_capture_panel1/data/reads/*.bbduk.lowcomp_filtered.fq
  /home/scratch/yinmt202607/4.mcp_reshotgun/data/reads/*.bbduk.lowcomp_filtered.fq
  /home/scratch/yinmt202607/7_angor_capture_panel2/data/reads/*.bbduk.lowcomp_filtered.fq

Examples:
  Run locally:
    bash compare_rice_read_extractors.sh --threads 32

  Submit to Slurm:
    bash compare_rice_read_extractors.sh --sbatch --threads 16

  Submit to a specific partition:
    bash compare_rice_read_extractors.sh --sbatch --threads 16 --sbatch-partition your_partition

  Restrict jobs to specific nodes:
    bash compare_rice_read_extractors.sh --sbatch --threads 16 --sbatch-nodelist node01,node02,node03

  Collect results after jobs finish:
    bash compare_rice_read_extractors.sh --collect-only

  Write unfinished tasks to a list:
    bash compare_rice_read_extractors.sh --write-pending-list pending.txt

  Split unfinished tasks into 4 lists:
    bash compare_rice_read_extractors.sh --split-pending 4 --split-prefix pending

  Run one list in a single shell:
    bash compare_rice_read_extractors.sh --task-list pending.1.txt
EOF
}

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2
}

die() {
  log "ERROR: $*"
  exit 1
}

DEFAULT_READ_GLOBS=(
  "/home/scratch/yinmt202607/3.angkor_capture_panel1/data/reads/*.bbduk.lowcomp_filtered.fq"
  "/home/scratch/yinmt202607/4.mcp_reshotgun/data/reads/*.bbduk.lowcomp_filtered.fq"
  "/home/scratch/yinmt202607/7_angor_capture_panel2/data/reads/*.bbduk.lowcomp_filtered.fq"
)

READ_GLOBS=("${DEFAULT_READ_GLOBS[@]}")
TASK_LIST_FILE=""
REF_FASTA="/home/scratch/yinmt202607/db/asian_rice_panel_index/asian_rice_panel.fa"
OUT_DIR="/home/scratch/yinmt202607/results/asian_rice_compare"
THREADS=20
COUNT_INPUT_READS=0
FORCE=0
USE_MODULES=1
BT2_PREFIX_SET=0
MODE="local"
RUN_FASTQ=""
WRITE_PENDING_LIST=""
SPLIT_PENDING=0
SPLIT_PREFIX=""
SBATCH_PARTITION=""
SBATCH_NODELIST=""
SBATCH_TIME="08:00:00"
SBATCH_MEM="120G"
SUBMIT_COLLECTOR=1

BOWTIE2_EXTRA=(
  -k 3
  -L 22
  -i S,1,1.15
  --mp 1,1
  --rdg 0,1
  --rfg 0,1
  --score-min L,0,-0.1
  --no-unal
)

BWA_ALN_EXTRA=(
  -l 1024
  -n 0.01
  -o 2
)

PRIMARY_MAPPED_FILTER="0x904"
SUMMARY_HEADER=$'sample\tinput_fastq\tinput_reads\tbowtie2_mapped_reads\tbwa_mapped_reads\tdiff_bowtie2_minus_bwa\twinner\tbowtie2_bam\tbowtie2_fastq\tbwa_bam\tbwa_fastq'

while (( $# > 0 )); do
  case "$1" in
    --read-glob)
      if [[ "${READ_GLOBS[*]}" == "${DEFAULT_READ_GLOBS[*]}" ]]; then
        READ_GLOBS=()
      fi
      [[ $# -ge 2 ]] || die "--read-glob requires a value"
      READ_GLOBS+=("$2")
      shift 2
      ;;
    --task-list)
      [[ $# -ge 2 ]] || die "--task-list requires a value"
      TASK_LIST_FILE="$2"
      shift 2
      ;;
    --ref-fasta)
      [[ $# -ge 2 ]] || die "--ref-fasta requires a value"
      REF_FASTA="$2"
      shift 2
      ;;
    --bt2-prefix)
      [[ $# -ge 2 ]] || die "--bt2-prefix requires a value"
      BOWTIE2_PREFIX="$2"
      BT2_PREFIX_SET=1
      shift 2
      ;;
    --out-dir)
      [[ $# -ge 2 ]] || die "--out-dir requires a value"
      OUT_DIR="$2"
      shift 2
      ;;
    --threads)
      [[ $# -ge 2 ]] || die "--threads requires a value"
      THREADS="$2"
      shift 2
      ;;
    --count-input-reads)
      COUNT_INPUT_READS=1
      shift
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --no-modules)
      USE_MODULES=0
      shift
      ;;
    --sbatch)
      MODE="sbatch"
      shift
      ;;
    --sbatch-partition)
      [[ $# -ge 2 ]] || die "--sbatch-partition requires a value"
      SBATCH_PARTITION="$2"
      shift 2
      ;;
    --sbatch-nodelist)
      [[ $# -ge 2 ]] || die "--sbatch-nodelist requires a value"
      SBATCH_NODELIST="$2"
      shift 2
      ;;
    --sbatch-time)
      [[ $# -ge 2 ]] || die "--sbatch-time requires a value"
      SBATCH_TIME="$2"
      shift 2
      ;;
    --sbatch-mem)
      [[ $# -ge 2 ]] || die "--sbatch-mem requires a value"
      SBATCH_MEM="$2"
      shift 2
      ;;
    --no-collector-job)
      SUBMIT_COLLECTOR=0
      shift
      ;;
    --collect-only)
      MODE="collect"
      shift
      ;;
    --write-pending-list)
      [[ $# -ge 2 ]] || die "--write-pending-list requires a value"
      MODE="write-pending-list"
      WRITE_PENDING_LIST="$2"
      shift 2
      ;;
    --split-pending)
      [[ $# -ge 2 ]] || die "--split-pending requires a value"
      MODE="split-pending"
      SPLIT_PENDING="$2"
      shift 2
      ;;
    --split-prefix)
      [[ $# -ge 2 ]] || die "--split-prefix requires a value"
      SPLIT_PREFIX="$2"
      shift 2
      ;;
    --run-bowtie2)
      [[ $# -ge 2 ]] || die "--run-bowtie2 requires a value"
      MODE="run-bowtie2"
      RUN_FASTQ="$2"
      shift 2
      ;;
    --run-bwa)
      [[ $# -ge 2 ]] || die "--run-bwa requires a value"
      MODE="run-bwa"
      RUN_FASTQ="$2"
      shift 2
      ;;
    --finalize-sample)
      [[ $# -ge 2 ]] || die "--finalize-sample requires a value"
      MODE="finalize-sample"
      RUN_FASTQ="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "Unknown option: $1"
      ;;
  esac
done

if (( BT2_PREFIX_SET == 0 )); then
  BOWTIE2_PREFIX="$REF_FASTA"
fi

SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
SUMMARY_TSV="$OUT_DIR/summary.tsv"
TOTALS_TSV="$OUT_DIR/summary_totals.tsv"
STATS_DIR="$OUT_DIR/per_sample_stats"
SBATCH_LOG_DIR="$OUT_DIR/slurm_logs"
SUBMITTED_JOBS_TSV="$OUT_DIR/sbatch_jobs.tsv"

maybe_enable_module_cmd() {
  if type module >/dev/null 2>&1; then
    return 0
  fi

  local init_script
  for init_script in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
    if [[ -r "$init_script" ]]; then
      # shellcheck disable=SC1090
      source "$init_script"
      break
    fi
  done
}

load_modules() {
  (( USE_MODULES == 1 )) || return 0

  maybe_enable_module_cmd
  if ! type module >/dev/null 2>&1; then
    log "Environment modules not available; using tools already on PATH"
    return 0
  fi

  module load bowtie2/ 2>/dev/null || module load bowtie2 2>/dev/null || true
  module load samtools/ 2>/dev/null || module load samtools 2>/dev/null || true
}

require_cmd() {
  local cmd=$1
  command -v "$cmd" >/dev/null 2>&1 || die "Required command not found: $cmd"
}

check_ref_and_indexes() {
  [[ -f "$REF_FASTA" ]] || die "Reference FASTA not found: $REF_FASTA"

  local missing_bwa=0
  local suffix
  for suffix in amb ann bwt pac sa; do
    if [[ ! -f "${REF_FASTA}.${suffix}" ]]; then
      missing_bwa=1
      break
    fi
  done
  (( missing_bwa == 0 )) || die "Missing BWA index files for: $REF_FASTA"

  if [[ ! -f "${BOWTIE2_PREFIX}.1.bt2" && ! -f "${BOWTIE2_PREFIX}.1.bt2l" ]]; then
    die "Missing Bowtie2 index files for prefix: $BOWTIE2_PREFIX"
  fi
}

sample_name_from_path() {
  local name
  name=$(basename "$1")
  name=${name%.fastq.gz}
  name=${name%.fq.gz}
  name=${name%.fastq}
  name=${name%.fq}
  name=${name%.cleaned}
  name=${name%.bbduk.lowcomp_filtered}
  printf '%s\n' "$name"
}

count_fastq_reads() {
  local fq=$1

  if [[ "$fq" == *.gz ]]; then
    gzip -cd "$fq" | awk 'END {printf "%.0f\n", NR / 4}'
  else
    awk 'END {printf "%.0f\n", NR / 4}' "$fq"
  fi
}

count_bam_reads() {
  local bam=$1
  samtools view -c "$bam"
}

set_sample_paths() {
  local sample=$1

  SAMPLE_NAME="$sample"
  BOWTIE2_BAM="$OUT_DIR/bowtie2/bam/${sample}.bowtie2.primary_mapped.bam"
  BOWTIE2_FASTQ="$OUT_DIR/bowtie2/fastq/${sample}.bowtie2.primary_mapped.fastq.gz"
  BOWTIE2_LOG="$OUT_DIR/logs/bowtie2/${sample}.log"

  BWA_BAM="$OUT_DIR/bwa/bam/${sample}.bwa.primary_mapped.bam"
  BWA_FASTQ="$OUT_DIR/bwa/fastq/${sample}.bwa.primary_mapped.fastq.gz"
  BWA_LOG="$OUT_DIR/logs/bwa/${sample}.log"

  SAMPLE_STATS_TSV="$STATS_DIR/${sample}.tsv"
}

init_output_dirs() {
  mkdir -p "$OUT_DIR" "$STATS_DIR" "$SBATCH_LOG_DIR"
}

mapper_outputs_exist() {
  local mapper=$1

  case "$mapper" in
    bowtie2)
      [[ -s "$BOWTIE2_BAM" && -s "$BOWTIE2_FASTQ" ]]
      ;;
    bwa)
      [[ -s "$BWA_BAM" && -s "$BWA_FASTQ" ]]
      ;;
    *)
      die "Unknown mapper: $mapper"
      ;;
  esac
}

sample_outputs_exist() {
  mapper_outputs_exist bowtie2 && mapper_outputs_exist bwa && [[ -s "$SAMPLE_STATS_TSV" ]]
}

run_bowtie2() {
  local fq=$1
  local bam=$2
  local mapped_fastq=$3
  local log_file=$4
  local tmp_bam="${bam}.tmp"
  local tmp_fastq="${mapped_fastq}.tmp"

  if (( FORCE == 0 )) && [[ -s "$bam" && -s "$mapped_fastq" ]]; then
    log "Reuse existing Bowtie2 outputs for $(basename "$fq")"
    return 0
  fi

  rm -f "$tmp_bam" "$tmp_fastq" "$log_file"
  mkdir -p "$(dirname "$bam")" "$(dirname "$mapped_fastq")" "$(dirname "$log_file")"

  (
    bowtie2 \
      -p "$THREADS" \
      --time \
      "${BOWTIE2_EXTRA[@]}" \
      -x "$BOWTIE2_PREFIX" \
      -U "$fq" \
    | samtools view \
      -@ "$THREADS" \
      -bh \
      -F "$PRIMARY_MAPPED_FILTER" \
      -o "$tmp_bam" \
      -
  ) 2> "$log_file"

  samtools fastq -@ "$THREADS" "$tmp_bam" | gzip -c > "$tmp_fastq"

  mv "$tmp_bam" "$bam"
  mv "$tmp_fastq" "$mapped_fastq"
}

run_bwa() {
  local fq=$1
  local bam=$2
  local mapped_fastq=$3
  local log_file=$4
  local tmp_bam="${bam}.tmp"
  local tmp_fastq="${mapped_fastq}.tmp"

  if (( FORCE == 0 )) && [[ -s "$bam" && -s "$mapped_fastq" ]]; then
    log "Reuse existing BWA outputs for $(basename "$fq")"
    return 0
  fi

  rm -f "$tmp_bam" "$tmp_fastq" "$log_file"
  mkdir -p "$(dirname "$bam")" "$(dirname "$mapped_fastq")" "$(dirname "$log_file")"

  (
    bwa aln \
      "${BWA_ALN_EXTRA[@]}" \
      -t "$THREADS" \
      "$REF_FASTA" \
      "$fq" 2> "$log_file" \
    | bwa samse "$REF_FASTA" - "$fq" 2>> "$log_file" \
    | samtools view \
      -@ "$THREADS" \
      -bh \
      -F "$PRIMARY_MAPPED_FILTER" \
      -o "$tmp_bam" \
      -
  )

  samtools fastq -@ "$THREADS" "$tmp_bam" | gzip -c > "$tmp_fastq"

  mv "$tmp_bam" "$bam"
  mv "$tmp_fastq" "$mapped_fastq"
}

write_sample_stats() {
  local sample=$1
  local fq=$2
  local input_reads=$3
  local bowtie2_count=$4
  local bwa_count=$5
  local winner=$6
  local diff=$7

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" \
    "$fq" \
    "$input_reads" \
    "$bowtie2_count" \
    "$bwa_count" \
    "$diff" \
    "$winner" \
    "$BOWTIE2_BAM" \
    "$BOWTIE2_FASTQ" \
    "$BWA_BAM" \
    "$BWA_FASTQ" \
    > "$SAMPLE_STATS_TSV"
}

finalize_one_fastq() {
  local fq=$1
  local sample bowtie2_count bwa_count input_reads winner diff

  [[ -f "$fq" ]] || die "FASTQ not found: $fq"

  sample=$(sample_name_from_path "$fq")
  set_sample_paths "$sample"

  [[ -s "$BOWTIE2_BAM" ]] || die "Missing Bowtie2 BAM for $sample: $BOWTIE2_BAM"
  [[ -s "$BOWTIE2_FASTQ" ]] || die "Missing Bowtie2 FASTQ for $sample: $BOWTIE2_FASTQ"
  [[ -s "$BWA_BAM" ]] || die "Missing BWA BAM for $sample: $BWA_BAM"
  [[ -s "$BWA_FASTQ" ]] || die "Missing BWA FASTQ for $sample: $BWA_FASTQ"

  bowtie2_count=$(count_bam_reads "$BOWTIE2_BAM")
  bwa_count=$(count_bam_reads "$BWA_BAM")

  if (( COUNT_INPUT_READS == 1 )); then
    input_reads=$(count_fastq_reads "$fq")
  else
    input_reads="NA"
  fi

  if (( bowtie2_count > bwa_count )); then
    winner="bowtie2"
  elif (( bwa_count > bowtie2_count )); then
    winner="bwa"
  else
    winner="tie"
  fi

  diff=$((bowtie2_count - bwa_count))
  write_sample_stats "$sample" "$fq" "$input_reads" "$bowtie2_count" "$bwa_count" "$winner" "$diff"
}

process_one_fastq() {
  local fq=$1
  local sample

  [[ -f "$fq" ]] || die "FASTQ not found: $fq"

  sample=$(sample_name_from_path "$fq")
  set_sample_paths "$sample"

  log "Processing $sample"
  run_bowtie2 "$fq" "$BOWTIE2_BAM" "$BOWTIE2_FASTQ" "$BOWTIE2_LOG"
  run_bwa "$fq" "$BWA_BAM" "$BWA_FASTQ" "$BWA_LOG"
  finalize_one_fastq "$fq"
}

gather_fastq_files() {
  FASTQ_FILES=()
  local -a sample_names=()
  local -a sample_paths=()
  local pattern fq sample seen_fq seen_idx already_seen line

  if [[ -n "$TASK_LIST_FILE" ]]; then
    [[ -f "$TASK_LIST_FILE" ]] || die "Task list not found: $TASK_LIST_FILE"

    while IFS= read -r line || [[ -n "$line" ]]; do
      line="${line#"${line%%[![:space:]]*}"}"
      line="${line%"${line##*[![:space:]]}"}"

      if [[ -z "$line" || "${line:0:1}" == "#" ]]; then
        continue
      fi

      fq="$line"
      already_seen=0
      for seen_fq in "${FASTQ_FILES[@]}"; do
        if [[ "$seen_fq" == "$fq" ]]; then
          already_seen=1
          break
        fi
      done

      if (( already_seen == 0 )); then
        FASTQ_FILES+=("$fq")
      fi
    done < "$TASK_LIST_FILE"
  else
    for pattern in "${READ_GLOBS[@]}"; do
      matches=( $pattern )
      if (( ${#matches[@]} == 0 )); then
        log "No FASTQ matched pattern: $pattern"
        continue
      fi

      for fq in "${matches[@]}"; do
        already_seen=0
        for seen_fq in "${FASTQ_FILES[@]}"; do
          if [[ "$seen_fq" == "$fq" ]]; then
            already_seen=1
            break
          fi
        done

        if (( already_seen == 0 )); then
          FASTQ_FILES+=("$fq")
        fi
      done
    done
  fi

  (( ${#FASTQ_FILES[@]} > 0 )) || die "No FASTQ files found. Check your --read-glob or default paths."

  for fq in "${FASTQ_FILES[@]}"; do
    [[ -f "$fq" ]] || die "FASTQ not found: $fq"
    sample=$(sample_name_from_path "$fq")
    for seen_idx in "${!sample_names[@]}"; do
      if [[ "${sample_names[$seen_idx]}" == "$sample" ]]; then
        die "Duplicate sample name after suffix trimming: $sample
  first: ${sample_paths[$seen_idx]}
  second: $fq
Use more specific filenames or split the runs."
      fi
    done

    sample_names+=("$sample")
    sample_paths+=("$fq")
  done
}

write_pending_list() {
  local output_file=$1
  local fq sample pending_count

  pending_count=0
  : > "$output_file"

  for fq in "${FASTQ_FILES[@]}"; do
    sample=$(sample_name_from_path "$fq")
    set_sample_paths "$sample"

    if ! sample_outputs_exist; then
      printf '%s\n' "$fq" >> "$output_file"
      pending_count=$((pending_count + 1))
    fi
  done

  log "Wrote $pending_count pending tasks to $output_file"
}

split_pending_lists() {
  local shard_count=$1
  local prefix=$2
  local fq sample shard_index pending_count i

  [[ "$shard_count" =~ ^[1-9][0-9]*$ ]] || die "--split-pending requires a positive integer"
  [[ -n "$prefix" ]] || die "--split-prefix is required with --split-pending"

  for (( i = 1; i <= shard_count; i++ )); do
    : > "${prefix}.${i}.txt"
  done

  pending_count=0
  for fq in "${FASTQ_FILES[@]}"; do
    sample=$(sample_name_from_path "$fq")
    set_sample_paths "$sample"

    if ! sample_outputs_exist; then
      shard_index=$(( (pending_count % shard_count) + 1 ))
      printf '%s\n' "$fq" >> "${prefix}.${shard_index}.txt"
      pending_count=$((pending_count + 1))
    fi
  done

  log "Split $pending_count pending tasks into $shard_count files with prefix $prefix"
}

collect_summaries() {
  local stats_files sample_count total_bowtie2 total_bwa total_winner

  init_output_dirs
  stats_files=( "$STATS_DIR"/*.tsv )
  (( ${#stats_files[@]} > 0 )) || die "No per-sample stats found in $STATS_DIR"

  printf '%s\n' "$SUMMARY_HEADER" > "$SUMMARY_TSV"
  printf '%s\n' "${stats_files[@]}" | sort | while IFS= read -r stats_file; do
    cat "$stats_file"
  done >> "$SUMMARY_TSV"

  read -r sample_count total_bowtie2 total_bwa < <(
    awk -F '\t' '{
      n += 1
      bt += $4
      bwa += $5
    }
    END {
      printf "%d %d %d\n", n, bt, bwa
    }' "$STATS_DIR"/*.tsv
  )

  if (( total_bowtie2 > total_bwa )); then
    total_winner="bowtie2"
  elif (( total_bwa > total_bowtie2 )); then
    total_winner="bwa"
  else
    total_winner="tie"
  fi

  printf 'samples\tbowtie2_total_mapped_reads\tbwa_total_mapped_reads\tdiff_bowtie2_minus_bwa\twinner\n' > "$TOTALS_TSV"
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$sample_count" \
    "$total_bowtie2" \
    "$total_bwa" \
    "$((total_bowtie2 - total_bwa))" \
    "$total_winner" \
    >> "$TOTALS_TSV"

  log "Finished. Sample summary: $SUMMARY_TSV"
  log "Finished. Total summary:  $TOTALS_TSV"
}

submit_mapper_job() {
  local mapper=$1
  local fq=$2
  local sample=$3
  local mode_arg log_path job_name job_id
  local -a submit_cmd=()
  local -a job_args=()

  case "$mapper" in
    bowtie2)
      mode_arg="--run-bowtie2"
      log_path="$SBATCH_LOG_DIR/${sample}.bowtie2.%j.out"
      job_name="ricecmp_bt2_${sample}"
      ;;
    bwa)
      mode_arg="--run-bwa"
      log_path="$SBATCH_LOG_DIR/${sample}.bwa.%j.out"
      job_name="ricecmp_bwa_${sample}"
      ;;
    *)
      die "Unknown mapper: $mapper"
      ;;
  esac

  job_args=(
    "$mode_arg" "$fq"
    --ref-fasta "$REF_FASTA"
    --bt2-prefix "$BOWTIE2_PREFIX"
    --out-dir "$OUT_DIR"
    --threads "$THREADS"
  )

  if (( FORCE == 1 )); then
    job_args+=( --force )
  fi

  if (( USE_MODULES == 0 )); then
    job_args+=( --no-modules )
  fi

  submit_cmd=(
    sbatch
    --parsable
    --job-name "$job_name"
    --cpus-per-task "$THREADS"
    --output "$log_path"
  )

  if [[ -n "$SBATCH_MEM" ]]; then
    submit_cmd+=( --mem "$SBATCH_MEM" )
  fi

  if [[ -n "$SBATCH_TIME" ]]; then
    submit_cmd+=( --time "$SBATCH_TIME" )
  fi

  if [[ -n "$SBATCH_PARTITION" ]]; then
    submit_cmd+=( --partition "$SBATCH_PARTITION" )
  fi

  if [[ -n "$SBATCH_NODELIST" ]]; then
    submit_cmd+=( --nodelist "$SBATCH_NODELIST" )
  fi

  submit_cmd+=( "$SCRIPT_PATH" "${job_args[@]}" )

  job_id=$("${submit_cmd[@]}")
  printf '%s\n' "${job_id%%;*}"
}

submit_finalize_job() {
  local fq=$1
  local sample=$2
  local dependency_ids=$3
  local job_id
  local -a submit_cmd=()
  local -a job_args=()

  job_args=(
    --finalize-sample "$fq"
    --out-dir "$OUT_DIR"
  )

  if (( COUNT_INPUT_READS == 1 )); then
    job_args+=( --count-input-reads )
  fi

  if (( USE_MODULES == 0 )); then
    job_args+=( --no-modules )
  fi

  submit_cmd=(
    sbatch
    --parsable
    --job-name "ricecmp_finalize_${sample}"
    --cpus-per-task 1
    --mem 4G
    --time 00:30:00
    --output "$SBATCH_LOG_DIR/${sample}.finalize.%j.out"
  )

  if [[ -n "$SBATCH_PARTITION" ]]; then
    submit_cmd+=( --partition "$SBATCH_PARTITION" )
  fi

  if [[ -n "$SBATCH_NODELIST" ]]; then
    submit_cmd+=( --nodelist "$SBATCH_NODELIST" )
  fi

  if [[ -n "$dependency_ids" ]]; then
    submit_cmd+=( --dependency "afterok:${dependency_ids}" )
  fi

  submit_cmd+=( "$SCRIPT_PATH" "${job_args[@]}" )

  job_id=$("${submit_cmd[@]}")
  printf '%s\n' "${job_id%%;*}"
}

submit_collector_job() {
  local dependency_ids=$1
  local job_id
  local -a submit_cmd=()

  submit_cmd=(
    sbatch
    --parsable
    --job-name "ricecmp_collect"
    --cpus-per-task 1
    --mem 2G
    --time 00:20:00
    --output "$SBATCH_LOG_DIR/collect.%j.out"
  )

  if [[ -n "$SBATCH_PARTITION" ]]; then
    submit_cmd+=( --partition "$SBATCH_PARTITION" )
  fi

  if [[ -n "$SBATCH_NODELIST" ]]; then
    submit_cmd+=( --nodelist "$SBATCH_NODELIST" )
  fi

  if [[ -n "$dependency_ids" ]]; then
    submit_cmd+=( --dependency "afterok:${dependency_ids}" )
  fi

  submit_cmd+=(
    "$SCRIPT_PATH"
    --collect-only
    --out-dir "$OUT_DIR"
  )

  job_id=$("${submit_cmd[@]}")
  printf '%s\n' "${job_id%%;*}"
}

submit_sbatch_jobs() {
  local fq sample
  local bowtie2_job_id bwa_job_id finalize_job_id finalize_dependencies
  local submit_count finalize_count
  local -a finalize_job_ids=()

  init_output_dirs
  printf 'sample\tinput_fastq\tstep\tjob_id\n' > "$SUBMITTED_JOBS_TSV"
  submit_count=0
  finalize_count=0

  for fq in "${FASTQ_FILES[@]}"; do
    sample=$(sample_name_from_path "$fq")
    set_sample_paths "$sample"

    bowtie2_job_id=""
    bwa_job_id=""
    finalize_job_id=""
    finalize_dependencies=""

    if (( FORCE == 0 )) && sample_outputs_exist; then
      log "Skip already finished sample: $sample"
      continue
    fi

    if (( FORCE == 1 )) || ! mapper_outputs_exist bowtie2; then
      bowtie2_job_id=$(submit_mapper_job bowtie2 "$fq" "$sample")
      printf '%s\t%s\t%s\t%s\n' "$sample" "$fq" "bowtie2" "$bowtie2_job_id" >> "$SUBMITTED_JOBS_TSV"
      log "Submitted Bowtie2 job for $sample: $bowtie2_job_id"
      submit_count=$((submit_count + 1))
    else
      log "Reuse existing Bowtie2 outputs for $sample"
    fi

    if (( FORCE == 1 )) || ! mapper_outputs_exist bwa; then
      bwa_job_id=$(submit_mapper_job bwa "$fq" "$sample")
      printf '%s\t%s\t%s\t%s\n' "$sample" "$fq" "bwa" "$bwa_job_id" >> "$SUBMITTED_JOBS_TSV"
      log "Submitted BWA job for $sample: $bwa_job_id"
      submit_count=$((submit_count + 1))
    else
      log "Reuse existing BWA outputs for $sample"
    fi

    if [[ -n "$bowtie2_job_id" && -n "$bwa_job_id" ]]; then
      finalize_dependencies="${bowtie2_job_id}:${bwa_job_id}"
    elif [[ -n "$bowtie2_job_id" ]]; then
      finalize_dependencies="$bowtie2_job_id"
    elif [[ -n "$bwa_job_id" ]]; then
      finalize_dependencies="$bwa_job_id"
    fi

    if (( FORCE == 1 )) || [[ -n "$finalize_dependencies" ]] || [[ ! -s "$SAMPLE_STATS_TSV" ]]; then
      finalize_job_id=$(submit_finalize_job "$fq" "$sample" "$finalize_dependencies")
      finalize_job_ids+=("$finalize_job_id")
      printf '%s\t%s\t%s\t%s\n' "$sample" "$fq" "finalize" "$finalize_job_id" >> "$SUBMITTED_JOBS_TSV"
      log "Submitted finalize job for $sample: $finalize_job_id"
      finalize_count=$((finalize_count + 1))
    else
      log "Reuse existing sample stats for $sample"
    fi
  done

  if (( submit_count == 0 && finalize_count == 0 )); then
    log "No new jobs submitted."
    if compgen -G "$STATS_DIR/*.tsv" >/dev/null 2>&1; then
      log "Existing per-sample stats found, collecting summaries now."
      collect_summaries
    fi
    return 0
  fi

  if (( SUBMIT_COLLECTOR == 1 && finalize_count > 0 )); then
    collector_job_id=$(submit_collector_job "$(IFS=:; printf '%s' "${finalize_job_ids[*]}")")
    printf '%s\t%s\t%s\t%s\n' "ALL" "ALL" "collect" "$collector_job_id" >> "$SUBMITTED_JOBS_TSV"
    log "Submitted collector job: $collector_job_id"
  elif (( SUBMIT_COLLECTOR == 1 )); then
    log "No finalize jobs were submitted, so no collector job was needed."
  else
    log "Collector job disabled. Run 'bash compare_rice_read_extractors.sh --collect-only --out-dir $OUT_DIR' after sample jobs finish."
  fi

  log "Submitted jobs table: $SUBMITTED_JOBS_TSV"
}

run_local_all_samples() {
  local fq total sample_idx

  init_output_dirs
  total=${#FASTQ_FILES[@]}
  sample_idx=0

  for fq in "${FASTQ_FILES[@]}"; do
    sample_idx=$((sample_idx + 1))
    log "[$sample_idx/$total] $(sample_name_from_path "$fq")"
    process_one_fastq "$fq"
  done

  collect_summaries
}

validate_processing_prereqs() {
  load_modules
  require_cmd bowtie2
  require_cmd bwa
  require_cmd samtools
  require_cmd gzip
  require_cmd awk
}

validate_finalize_prereqs() {
  load_modules
  require_cmd samtools
  require_cmd gzip
  require_cmd awk
}

case "$MODE" in
  write-pending-list)
    gather_fastq_files
    write_pending_list "$WRITE_PENDING_LIST"
    ;;
  split-pending)
    gather_fastq_files
    split_pending_lists "$SPLIT_PENDING" "$SPLIT_PREFIX"
    ;;
  collect)
    collect_summaries
    ;;
  run-bowtie2)
    check_ref_and_indexes
    validate_processing_prereqs
    init_output_dirs
    sample=$(sample_name_from_path "$RUN_FASTQ")
    set_sample_paths "$sample"
    run_bowtie2 "$RUN_FASTQ" "$BOWTIE2_BAM" "$BOWTIE2_FASTQ" "$BOWTIE2_LOG"
    ;;
  run-bwa)
    check_ref_and_indexes
    validate_processing_prereqs
    init_output_dirs
    sample=$(sample_name_from_path "$RUN_FASTQ")
    set_sample_paths "$sample"
    run_bwa "$RUN_FASTQ" "$BWA_BAM" "$BWA_FASTQ" "$BWA_LOG"
    ;;
  finalize-sample)
    validate_finalize_prereqs
    init_output_dirs
    finalize_one_fastq "$RUN_FASTQ"
    ;;
  sbatch)
    require_cmd sbatch
    check_ref_and_indexes
    gather_fastq_files
    submit_sbatch_jobs
    ;;
  local)
    check_ref_and_indexes
    validate_processing_prereqs
    gather_fastq_files
    run_local_all_samples
    ;;
  *)
    die "Unsupported mode: $MODE"
    ;;
esac

