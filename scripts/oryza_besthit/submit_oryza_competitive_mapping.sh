#!/usr/bin/env bash
#
# Submit and run the Oryza competitive-mapping stage on SLURM.
#
# check mode validates every input/index and calculates the memory range without
# submitting jobs. test mode submits a smoke test using a small FASTQ subset.
# submit mode:
#   - discovers *.oryza_candidates.combined.fastq.gz samples
#   - calculates --mem independently for every Bowtie2 database from the six
#     index files: mem_GiB = (0.7 * index_GiB + 40) * MEMORY_MULTIPLIER
#   - submits one mapping job per sample x database
#   - submits one merge + query-name-sort job per sample with afterok dependencies
#
# Mapping databases:
#   - wgs_eukaryota.1 .. wgs_eukaryota.129
#   - asian_rice_panel
#   - IRGSP
#
# The script is also its own SLURM worker (internal map/merge modes).

set -Eeuo pipefail

report_error() {
    local exit_code="$?"
    echo "ERROR: line ${BASH_LINENO[0]} failed with exit ${exit_code}: ${BASH_COMMAND}" >&2
    exit "$exit_code"
}
trap report_error ERR

# -----------------------------------------------------------------------------
# Paths. Environment variables with the same names override these defaults.
# -----------------------------------------------------------------------------

READ_DIR="${READ_DIR:-/home/scratch/yinmt202607/gene/results/oryza_candidates_combined}"
READ_SUFFIX="${READ_SUFFIX:-.oryza_candidates.combined.fastq.gz}"

WGS_DB_DIR="${WGS_DB_DIR:-/home/database/ref20250728/cph_euk}"
RICE_DB_DIR="${RICE_DB_DIR:-/home/scratch/yinmt202607/db/asian_rice_panel_index}"

OUT_DIR="${OUT_DIR:-/home/scratch/yinmt202607/gene/results/oryza_competitive_mapping}"
BAM_DIR="${OUT_DIR}/bam_by_database"
FINAL_DIR="${OUT_DIR}/by_sample"
LOG_DIR="${OUT_DIR}/logs"
SUBMIT_DIR="${OUT_DIR}/submissions"

# -----------------------------------------------------------------------------
# SLURM resources. Defaults mirror new_single_multi/step4.euk.mapping.smk.
# Override at submission time if needed, for example:
#   SLURM_PARTITION=compregular SLURM_ACCOUNT= bash script.sh submit
# -----------------------------------------------------------------------------

# These defaults match the current cluster shown by `squeue`. Leave account
# empty to use the user's default SLURM association.
SLURM_PARTITION="${SLURM_PARTITION:-comp}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"

MAP_CPUS="${MAP_CPUS:-20}"
BOWTIE2_THREADS="${BOWTIE2_THREADS:-18}"
SAMTOOLS_VIEW_THREADS="${SAMTOOLS_VIEW_THREADS:-1}"
MAP_TIME="${MAP_TIME:-04:00:00}"

MERGE_CPUS="${MERGE_CPUS:-16}"
MERGE_THREADS="${MERGE_THREADS:-2}"
SORT_THREADS="${SORT_THREADS:-11}"
SORT_MEM_PER_THREAD="${SORT_MEM_PER_THREAD:-2G}"
MERGE_MEM_MB="${MERGE_MEM_MB:-102400}"
MERGE_TIME="${MERGE_TIME:-24:00:00}"

# First submission: 1. Failed/OOM jobs can be resubmitted with 2 or 3.
MEMORY_MULTIPLIER="${MEMORY_MULTIPLIER:-1}"

# Smoke-test settings. The test is isolated under OUT_DIR/smoke_test and uses
# the first discovered sample unless TEST_FASTQ is supplied.
TEST_READS="${TEST_READS:-1000}"
TEST_DB="${TEST_DB:-irgsp}"
TEST_FASTQ="${TEST_FASTQ:-}"
TEST_OUT_DIR="${TEST_OUT_DIR:-${OUT_DIR}/smoke_test}"
TEST_MAP_TIME="${TEST_MAP_TIME:-01:00:00}"

# Bowtie2 parameters from new_single_multi/step4.euk.mapping.smk.
BOWTIE2_EXTRA=(
    -k 100
    -L 22
    -i S,1,1.15
    --mp 1,1
    --rdg 0,1
    --rfg 0,1
    --score-min L,0,-0.1
    --no-unal
)

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"

DBS=()
for shard in $(seq 1 129); do
    DBS+=("wgs_eukaryota.${shard}")
done
DBS+=("asian_rice_panel" "irgsp")

usage() {
    cat <<'EOF'
Usage:
  # 1. Preflight only: validates all 131 indexes; submits nothing.
  bash submit_oryza_competitive_mapping.sh check

  # 2. Smoke test: submit only the first 1,000 reads against IRGSP.
  bash submit_oryza_competitive_mapping.sh test

  # 3. Full submission is locked until it is explicitly confirmed.
  CONFIRM_FULL_SUBMIT=YES bash submit_oryza_competitive_mapping.sh submit

Optional smoke-test overrides:
  TEST_READS=100 TEST_DB=asian_rice_panel \
    bash submit_oryza_competitive_mapping.sh test
  TEST_DB=wgs_eukaryota.1 bash submit_oryza_competitive_mapping.sh test
  TEST_FASTQ=/absolute/path/sample.fastq.gz \
    bash submit_oryza_competitive_mapping.sh test

Internal SLURM worker modes (do not normally run by hand):
  bash submit_oryza_competitive_mapping.sh map SAMPLE DB FASTQ
  bash submit_oryza_competitive_mapping.sh merge SAMPLE

Do not rerun submit while earlier mapping jobs are still queued/running.
After an interruption, first inspect squeue and the submission TSV.

Retry failed/OOM mapping jobs with twice the calculated memory (only after the
previous jobs have finished/failed):
  MEMORY_MULTIPLIER=2 CONFIRM_FULL_SUBMIT=YES \
    bash submit_oryza_competitive_mapping.sh submit
EOF
}

timestamp() {
    date --iso-8601=seconds 2>/dev/null || date '+%Y-%m-%dT%H:%M:%S%z'
}

load_mapping_modules() {
    if command -v module >/dev/null 2>&1; then
        module load bowtie2/
        module load samtools/
    fi
    command -v bowtie2 >/dev/null 2>&1 || {
        echo "ERROR: bowtie2 not found on PATH" >&2
        exit 1
    }
    command -v samtools >/dev/null 2>&1 || {
        echo "ERROR: samtools not found on PATH" >&2
        exit 1
    }
}

load_samtools_module() {
    if command -v module >/dev/null 2>&1; then
        module load samtools/
    fi
    command -v samtools >/dev/null 2>&1 || {
        echo "ERROR: samtools not found on PATH" >&2
        exit 1
    }
}

db_prefix() {
    local db="$1"
    case "$db" in
        wgs_eukaryota.*)
            printf '%s/%s.fas.gz\n' "$WGS_DB_DIR" "$db"
            ;;
        asian_rice_panel)
            printf '%s/asian_rice_panel.fa\n' "$RICE_DB_DIR"
            ;;
        irgsp)
            printf '%s/irgsp_bt2idx\n' "$RICE_DB_DIR"
            ;;
        *)
            echo "ERROR: unknown database: $db" >&2
            return 1
            ;;
    esac
}

# Populates the global INDEX_FILES array for one Bowtie2 prefix.
resolve_index_files() {
    local prefix="$1"
    local suffix
    local all_present=1

    INDEX_FILES=()
    for suffix in .1.bt2l .2.bt2l .3.bt2l .4.bt2l .rev.1.bt2l .rev.2.bt2l; do
        INDEX_FILES+=("${prefix}${suffix}")
        [[ -f "${prefix}${suffix}" ]] || all_present=0
    done
    if [[ "$all_present" -eq 1 ]]; then
        return 0
    fi

    INDEX_FILES=()
    all_present=1
    for suffix in .1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2; do
        INDEX_FILES+=("${prefix}${suffix}")
        [[ -f "${prefix}${suffix}" ]] || all_present=0
    done
    if [[ "$all_present" -eq 1 ]]; then
        return 0
    fi

    echo "ERROR: no complete .bt2l or .bt2 index for prefix: $prefix" >&2
    return 1
}

index_bytes() {
    local prefix="$1"
    local total=0
    local file size
    resolve_index_files "$prefix"
    for file in "${INDEX_FILES[@]}"; do
        if size="$(stat -c '%s' -- "$file" 2>/dev/null)"; then
            :
        elif size="$(stat -f '%z' -- "$file" 2>/dev/null)"; then
            :
        else
            echo "ERROR: cannot determine file size: $file" >&2
            return 1
        fi
        total=$((total + size))
    done
    printf '%s\n' "$total"
}

calculated_mem_mb() {
    local bytes="$1"
    awk -v bytes="$bytes" -v multiplier="$MEMORY_MULTIPLIER" '
        BEGIN {
            index_gib = bytes / 1024 / 1024 / 1024
            mem_mib = (0.7 * index_gib + 40) * multiplier * 1024
            rounded = int(mem_mib)
            if (mem_mib > rounded) rounded++
            print rounded
        }
    '
}

sbatch_common_args() {
    SBATCH_COMMON=(--partition="$SLURM_PARTITION")
    if [[ -n "$SLURM_ACCOUNT" ]]; then
        SBATCH_COMMON+=(--account="$SLURM_ACCOUNT")
    fi
}

discover_fastqs() {
    [[ -d "$READ_DIR" ]] || {
        echo "ERROR: READ_DIR not found: $READ_DIR" >&2
        return 1
    }

    shopt -s nullglob
    FASTQS=("$READ_DIR"/*"$READ_SUFFIX")
    shopt -u nullglob
    [[ "${#FASTQS[@]}" -gt 0 ]] || {
        echo "ERROR: no FASTQ matched ${READ_DIR}/*${READ_SUFFIX}" >&2
        return 1
    }

    local fastq
    for fastq in "${FASTQS[@]}"; do
        [[ -s "$fastq" ]] || {
            echo "ERROR: FASTQ missing/empty: $fastq" >&2
            return 1
        }
    done
}

known_database() {
    local wanted="$1" db
    for db in "${DBS[@]}"; do
        [[ "$wanted" == "$db" ]] && return 0
    done
    return 1
}

validate_slurm_request() {
    local requested_mem_mb="${1:-100}"
    local requested_cpus="${2:-1}"
    local requested_time="${3:-00:01:00}"

    command -v sbatch >/dev/null 2>&1 || {
        echo "ERROR: sbatch not found on PATH" >&2
        return 1
    }
    sbatch_common_args

    if command -v sinfo >/dev/null 2>&1; then
        local available_partitions
        available_partitions="$(sinfo -h -o '%P' | tr -d '*' | sort -u)"
        if ! grep -Fxq "$SLURM_PARTITION" <<< "$available_partitions"; then
            echo "ERROR: SLURM partition does not exist: $SLURM_PARTITION" >&2
            echo "Available partitions:" >&2
            sed 's/^/  /' <<< "$available_partitions" >&2
            return 1
        fi
    fi

    # --test-only asks SLURM to validate partition/account/resources without
    # creating a queued job.
    sbatch \
        --test-only \
        "${SBATCH_COMMON[@]}" \
        --job-name=oryza.preflight \
        --cpus-per-task="$requested_cpus" \
        --mem="${requested_mem_mb}M" \
        --time="$requested_time" \
        --wrap=true >/dev/null
}

run_check() {
    [[ "$MEMORY_MULTIPLIER" =~ ^[0-9]+([.][0-9]+)?$ ]] || {
        echo "ERROR: MEMORY_MULTIPLIER must be numeric: $MEMORY_MULTIPLIER" >&2
        return 1
    }

    echo "[check] no jobs will be submitted"
    echo "[check] READ_DIR=$READ_DIR"
    echo "[check] WGS_DB_DIR=$WGS_DB_DIR"
    echo "[check] RICE_DB_DIR=$RICE_DB_DIR"
    echo "[check] partition=$SLURM_PARTITION account=${SLURM_ACCOUNT:-<default>}"

    discover_fastqs
    validate_slurm_request
    load_mapping_modules

    local db prefix bytes mem_mb
    local min_mem=0 max_mem=0 min_db='' max_db=''
    for db in "${DBS[@]}"; do
        prefix="$(db_prefix "$db")"
        bytes="$(index_bytes "$prefix")"
        mem_mb="$(calculated_mem_mb "$bytes")"
        if [[ "$min_mem" -eq 0 || "$mem_mb" -lt "$min_mem" ]]; then
            min_mem="$mem_mb"
            min_db="$db"
        fi
        if [[ "$mem_mb" -gt "$max_mem" ]]; then
            max_mem="$mem_mb"
            max_db="$db"
        fi
    done
    echo "[check] samples=${#FASTQS[@]} databases=${#DBS[@]}"
    echo "[check] first_sample=$(basename "${FASTQS[0]}")"
    echo "[check] minimum request=${min_mem}MiB (${min_db})"
    echo "[check] maximum request=${max_mem}MiB (${max_db})"
    validate_slurm_request "$max_mem" "$MAP_CPUS" "$MAP_TIME"
    echo "[check] largest mapping request accepted by SLURM test-only"
    echo "[check] PASS: inputs, tools, indexes, partition, and account validated"
}

create_test_fastq() {
    local source_fastq="$1" target_fastq="$2" reads="$3"
    local tmp_fastq="${target_fastq}.tmp.$$"
    local wanted_lines=$((reads * 4))
    local observed_lines

    mkdir -p "$(dirname "$target_fastq")"
    rm -f "$tmp_fastq"

    # The compressor is the final pipeline command. Disabling pipefail here
    # avoids treating gzip's expected SIGPIPE after `head` as a test failure.
    set +o pipefail
    gzip -cd -- "$source_fastq" | head -n "$wanted_lines" | gzip -c > "$tmp_fastq"
    local subset_status="$?"
    set -o pipefail
    [[ "$subset_status" -eq 0 ]] || {
        rm -f "$tmp_fastq"
        echo "ERROR: failed to create smoke-test FASTQ from $source_fastq" >&2
        return 1
    }

    observed_lines="$(gzip -cd -- "$tmp_fastq" | wc -l | awk '{print $1}')"
    [[ "$observed_lines" -ge 4 && $((observed_lines % 4)) -eq 0 ]] || {
        rm -f "$tmp_fastq"
        echo "ERROR: smoke-test FASTQ is incomplete (${observed_lines} lines)" >&2
        return 1
    }

    mv "$tmp_fastq" "$target_fastq"
    echo "[test] prepared_reads=$((observed_lines / 4)) test_fastq=$target_fastq"
}

submit_test() {
    [[ "$TEST_READS" =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: TEST_READS must be a positive integer: $TEST_READS" >&2
        exit 1
    }
    known_database "$TEST_DB" || {
        echo "ERROR: unknown TEST_DB: $TEST_DB" >&2
        exit 1
    }

    # Full preflight must pass before the single smoke-test job is submitted.
    run_check

    local source_fastq base source_sample test_sample test_fastq
    if [[ -n "$TEST_FASTQ" ]]; then
        source_fastq="$TEST_FASTQ"
        [[ -s "$source_fastq" ]] || {
            echo "ERROR: TEST_FASTQ missing/empty: $source_fastq" >&2
            exit 1
        }
        base="$(basename "$source_fastq")"
        source_sample="${base%%.*}"
    else
        source_fastq="${FASTQS[0]}"
        base="$(basename "$source_fastq")"
        source_sample="${base%$READ_SUFFIX}"
    fi
    test_sample="${source_sample}.smoke_test"

    # Redirect every test artifact away from production output.
    OUT_DIR="$TEST_OUT_DIR"
    BAM_DIR="${OUT_DIR}/bam_by_database"
    FINAL_DIR="${OUT_DIR}/by_sample"
    LOG_DIR="${OUT_DIR}/logs"
    SUBMIT_DIR="${OUT_DIR}/submissions"
    mkdir -p "$BAM_DIR" "$LOG_DIR/$test_sample" "$SUBMIT_DIR" "$OUT_DIR/input"

    test_fastq="${OUT_DIR}/input/${test_sample}.first_${TEST_READS}_reads.fastq.gz"
    create_test_fastq "$source_fastq" "$test_fastq" "$TEST_READS"

    local prefix bytes mem_mb index_gib map_log job_id clean_job_id
    prefix="$(db_prefix "$TEST_DB")"
    bytes="$(index_bytes "$prefix")"
    mem_mb="$(calculated_mem_mb "$bytes")"
    index_gib="$(awk -v b="$bytes" 'BEGIN{printf "%.2f", b/1024/1024/1024}')"
    map_log="${LOG_DIR}/${test_sample}/${test_sample}.${TEST_DB}.map.%j.log"
    sbatch_common_args

    job_id="$(sbatch \
        --parsable \
        "${SBATCH_COMMON[@]}" \
        --export="ALL,OUT_DIR=${OUT_DIR},WGS_DB_DIR=${WGS_DB_DIR},RICE_DB_DIR=${RICE_DB_DIR}" \
        --job-name="orytest.${TEST_DB}" \
        --cpus-per-task="$MAP_CPUS" \
        --mem="${mem_mb}M" \
        --time="$TEST_MAP_TIME" \
        --output="$map_log" \
        --error="$map_log" \
        "$SCRIPT_PATH" map "$test_sample" "$TEST_DB" "$test_fastq")"
    clean_job_id="${job_id%%;*}"

    printf 'job_id\tsample\tdatabase\tindex_GiB\tmem_MiB\tfastq\n' \
        > "${SUBMIT_DIR}/smoke_test.${clean_job_id}.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$clean_job_id" "$test_sample" "$TEST_DB" "$index_gib" "$mem_mb" "$test_fastq" \
        >> "${SUBMIT_DIR}/smoke_test.${clean_job_id}.tsv"

    echo "[test] submitted exactly one mapping job"
    echo "[test] job_id=$clean_job_id database=$TEST_DB mem=${mem_mb}MiB"
    echo "[test] production output was not touched"
    echo "[test] watch: squeue -j $clean_job_id"
    echo "[test] log:   tail -f ${map_log//%j/$clean_job_id}"
    echo "[test] BAM:   $(map_bam_path "$test_sample" "$TEST_DB")"
}

map_bam_path() {
    local sample="$1" db="$2"
    printf '%s/%s/%s.%s.bam\n' "$BAM_DIR" "$sample" "$sample" "$db"
}

map_done_path() {
    local sample="$1" db="$2"
    printf '%s/%s/%s.%s.bam.finished\n' "$BAM_DIR" "$sample" "$sample" "$db"
}

run_map() {
    [[ "$#" -eq 3 ]] || { usage >&2; exit 2; }
    local sample="$1" db="$2" fastq="$3"
    local prefix bam done tmp_bam

    prefix="$(db_prefix "$db")"
    resolve_index_files "$prefix"
    [[ -s "$fastq" ]] || { echo "ERROR: FASTQ missing/empty: $fastq" >&2; exit 1; }

    bam="$(map_bam_path "$sample" "$db")"
    done="$(map_done_path "$sample" "$db")"
    tmp_bam="${bam}.tmp.${SLURM_JOB_ID:-$$}"

    mkdir -p "$(dirname "$bam")"
    rm -f "$tmp_bam"
    trap 'rm -f "$tmp_bam"' EXIT

    load_mapping_modules

    echo "[map] sample=$sample db=$db"
    echo "[map] fastq=$fastq"
    echo "[map] index=$prefix"
    echo "[map] started=$(timestamp)"

    bowtie2 \
        --threads "$BOWTIE2_THREADS" \
        --time \
        -x "$prefix" \
        -U "$fastq" \
        "${BOWTIE2_EXTRA[@]}" \
    | samtools view \
        -@ "$SAMTOOLS_VIEW_THREADS" \
        -b \
        -o "$tmp_bam" \
        -

    samtools quickcheck -v "$tmp_bam"
    mv "$tmp_bam" "$bam"
    touch "$done"
    trap - EXIT

    echo "[map] finished=$(timestamp)"
    echo "[map] output=$bam"
}

run_merge() {
    [[ "$#" -eq 1 ]] || { usage >&2; exit 2; }
    local sample="$1"
    local final_bam="${FINAL_DIR}/${sample}.competitive.name_sorted.bam"
    local final_done="${FINAL_DIR}/${sample}.competitive.name_sorted.bam.finished"
    local tmp_bam="${final_bam}.tmp.${SLURM_JOB_ID:-$$}"
    local bam_list="${FINAL_DIR}/.${sample}.database_bams.list.${SLURM_JOB_ID:-$$}"
    local sort_tmp="${FINAL_DIR}/.${sample}.name_sort_tmp.${SLURM_JOB_ID:-$$}"
    local db bam done
    local input_bams=()

    mkdir -p "$FINAL_DIR"
    rm -f "$tmp_bam" "$bam_list"
    trap 'rm -f "$tmp_bam" "$bam_list"' EXIT

    load_samtools_module

    for db in "${DBS[@]}"; do
        bam="$(map_bam_path "$sample" "$db")"
        done="$(map_done_path "$sample" "$db")"
        [[ -s "$bam" && -f "$done" ]] || {
            echo "ERROR: mapping output incomplete: $bam" >&2
            exit 1
        }
        input_bams+=("$bam")
    done

    [[ "${#input_bams[@]}" -eq 131 ]] || {
        echo "ERROR: expected 131 BAMs, found ${#input_bams[@]}" >&2
        exit 1
    }

    samtools quickcheck -v "${input_bams[@]}"
    printf '%s\n' "${input_bams[@]}" > "$bam_list"

    echo "[merge] sample=$sample inputs=${#input_bams[@]}"
    echo "[merge] started=$(timestamp)"

    samtools merge \
        -@ "$MERGE_THREADS" \
        -u \
        -c \
        -p \
        --no-PG \
        -b "$bam_list" \
        -o - \
    | samtools sort \
        -n \
        -@ "$SORT_THREADS" \
        -m "$SORT_MEM_PER_THREAD" \
        -T "$sort_tmp" \
        -o "$tmp_bam" \
        -

    samtools quickcheck -v "$tmp_bam"
    if ! samtools view -H "$tmp_bam" | awk -F'\t' '
        $1 == "@HD" {
            for (i=2; i<=NF; i++) if ($i == "SO:queryname") ok=1
        }
        END { exit(ok ? 0 : 1) }
    '; then
        echo "ERROR: merged BAM header is not SO:queryname" >&2
        exit 1
    fi

    mv "$tmp_bam" "$final_bam"
    touch "$final_done"
    rm -f "$bam_list"
    trap - EXIT

    echo "[merge] finished=$(timestamp)"
    echo "[merge] output=$final_bam"
}

submit_all() {
    [[ "${CONFIRM_FULL_SUBMIT:-}" == "YES" ]] || {
        echo "ERROR: full submission is locked." >&2
        echo "Run 'bash $SCRIPT_PATH check', then 'bash $SCRIPT_PATH test' first." >&2
        echo "After the smoke test passes, use:" >&2
        echo "  CONFIRM_FULL_SUBMIT=YES bash $SCRIPT_PATH submit" >&2
        exit 2
    }

    # Repeat the complete no-submit preflight immediately before production.
    run_check

    mkdir -p "$BAM_DIR" "$FINAL_DIR" "$LOG_DIR" "$SUBMIT_DIR"
    sbatch_common_args

    local fastqs=("${FASTQS[@]}")

    # Calculate index sizes/memory once per database, not once per sample.
    declare -A DB_MEM_MB=()
    local plan_file="${SUBMIT_DIR}/index_memory_plan.multiplier_${MEMORY_MULTIPLIER}.tsv"
    printf 'database\tindex_GiB\trequested_mem_MiB\tmemory_multiplier\n' > "$plan_file"

    local db prefix bytes mem_mb index_gib
    for db in "${DBS[@]}"; do
        prefix="$(db_prefix "$db")"
        bytes="$(index_bytes "$prefix")"
        mem_mb="$(calculated_mem_mb "$bytes")"
        index_gib="$(awk -v b="$bytes" 'BEGIN{printf "%.2f", b/1024/1024/1024}')"
        DB_MEM_MB["$db"]="$mem_mb"
        printf '%s\t%s\t%s\t%s\n' \
            "$db" "$index_gib" "$mem_mb" "$MEMORY_MULTIPLIER" >> "$plan_file"
    done

    local stamp submission_log
    stamp="$(date '+%Y%m%d_%H%M%S')"
    submission_log="${SUBMIT_DIR}/submitted_jobs.${stamp}.tsv"
    printf 'job_type\tsample\tdatabase\tjob_id\tmem_MiB\tdependency\n' > "$submission_log"

    echo "[submit] samples=${#fastqs[@]} databases=${#DBS[@]}"
    echo "[submit] memory plan=$plan_file"
    echo "[submit] multiplier=$MEMORY_MULTIPLIER"

    local fq base sample bam done final_bam final_done
    local job_id clean_job_id dependency map_log merge_log short_sample
    local map_job_ids=()

    for fq in "${fastqs[@]}"; do
        base="$(basename "$fq")"
        sample="${base%$READ_SUFFIX}"
        [[ -n "$sample" ]] || { echo "ERROR: empty sample name from $fq" >&2; exit 1; }

        final_bam="${FINAL_DIR}/${sample}.competitive.name_sorted.bam"
        final_done="${FINAL_DIR}/${sample}.competitive.name_sorted.bam.finished"
        if [[ -s "$final_bam" && -f "$final_done" ]]; then
            echo "[submit] final already complete; skip sample=$sample"
            continue
        fi

        mkdir -p "${LOG_DIR}/${sample}" "${BAM_DIR}/${sample}"
        map_job_ids=()
        short_sample="${sample:0:28}"

        for db in "${DBS[@]}"; do
            bam="$(map_bam_path "$sample" "$db")"
            done="$(map_done_path "$sample" "$db")"
            if [[ -s "$bam" && -f "$done" ]]; then
                echo "[submit] mapping already complete; skip sample=$sample db=$db"
                continue
            fi

            map_log="${LOG_DIR}/${sample}/${sample}.${db}.map.%j.log"
            job_id="$(sbatch \
                --parsable \
                "${SBATCH_COMMON[@]}" \
                --job-name="orymap.${short_sample}.${db}" \
                --cpus-per-task="$MAP_CPUS" \
                --mem="${DB_MEM_MB[$db]}M" \
                --time="$MAP_TIME" \
                --output="$map_log" \
                --error="$map_log" \
                "$SCRIPT_PATH" map "$sample" "$db" "$fq")"
            clean_job_id="${job_id%%;*}"
            map_job_ids+=("$clean_job_id")
            printf 'map\t%s\t%s\t%s\t%s\t\n' \
                "$sample" "$db" "$clean_job_id" "${DB_MEM_MB[$db]}" >> "$submission_log"
            echo "[submit] map job=$clean_job_id sample=$sample db=$db mem=${DB_MEM_MB[$db]}MiB"
        done

        dependency=""
        if [[ "${#map_job_ids[@]}" -gt 0 ]]; then
            dependency="afterok:$(IFS=:; echo "${map_job_ids[*]}")"
        fi

        merge_log="${LOG_DIR}/${sample}/${sample}.merge_name_sort.%j.log"
        if [[ -n "$dependency" ]]; then
            job_id="$(sbatch \
                --parsable \
                "${SBATCH_COMMON[@]}" \
                --dependency="$dependency" \
                --kill-on-invalid-dep=yes \
                --job-name="orymerge.${short_sample}" \
                --cpus-per-task="$MERGE_CPUS" \
                --mem="${MERGE_MEM_MB}M" \
                --time="$MERGE_TIME" \
                --output="$merge_log" \
                --error="$merge_log" \
                "$SCRIPT_PATH" merge "$sample")"
        else
            job_id="$(sbatch \
                --parsable \
                "${SBATCH_COMMON[@]}" \
                --job-name="orymerge.${short_sample}" \
                --cpus-per-task="$MERGE_CPUS" \
                --mem="${MERGE_MEM_MB}M" \
                --time="$MERGE_TIME" \
                --output="$merge_log" \
                --error="$merge_log" \
                "$SCRIPT_PATH" merge "$sample")"
        fi
        clean_job_id="${job_id%%;*}"
        printf 'merge_sort\t%s\tALL\t%s\t%s\t%s\n' \
            "$sample" "$clean_job_id" "$MERGE_MEM_MB" "$dependency" >> "$submission_log"
        echo "[submit] merge job=$clean_job_id sample=$sample dependency_jobs=${#map_job_ids[@]}"
    done

    echo "[submit] job record=$submission_log"
    echo "[submit] finished=$(timestamp)"
}

mode="${1:-}"
case "$mode" in
    check)
        shift
        [[ "$#" -eq 0 ]] || { usage >&2; exit 2; }
        run_check
        ;;
    test)
        shift
        [[ "$#" -eq 0 ]] || { usage >&2; exit 2; }
        submit_test
        ;;
    submit)
        shift
        [[ "$#" -eq 0 ]] || { usage >&2; exit 2; }
        submit_all
        ;;
    map)
        shift
        run_map "$@"
        ;;
    merge)
        shift
        run_merge "$@"
        ;;
    -h|--help|help|"")
        usage
        ;;
    *)
        echo "ERROR: unknown mode: $mode" >&2
        usage >&2
        exit 2
        ;;
esac
