#!/bin/bash
#$ -N castor_kmer_sweep
#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs
#$ -S /bin/bash
#
# Repeated run of CASTOR-KRFE over every organism.
#
# The method selects the length of k itself inside the k_min..k_max range
# predefined in configuration_{organism}.ini, so there is no loop over k here:
# the loop is over the REPEATS repetitions, each one re-splitting the dataset so
# the repetitions are independent train/test partitions. The classification
# report of every run is appended to $OUTPUTS/reports.csv and aggregated into
# mean and standard deviation by aggregate_reports.py - the same metrics, the
# same files and the same format as the KEVOLVE sweep.
#
# Submission (all the organisms in a single job):
#     qsub doall_castor_kmer_sweep.sh
# Submission (one job per organism, using a SGE array job):
#     qsub -t 1-5 doall_castor_kmer_sweep.sh
# Submission of a subset of the organisms:
#     qsub doall_castor_kmer_sweep.sh denv hbv
#
# The job is resumable: an organism/k combination already completed is skipped,
# so the same job can be resubmitted after a wall clock kill.

set -u

FEHOME=/home/a61491
CASTOR=$FEHOME/CASTOR_KRFE
DATASETS=$FEHOME/datasets/original
BALANCEDATASET=$FEHOME/Fasta-splitter/FastaSplitter
TEMPROOT=/tmp2/felipe/castor

# Number of repetitions, and index of the first one (REPEAT_START lets a second
# job continue where a first one stopped)
REPEATS=${REPEATS:-100}
REPEAT_START=${REPEAT_START:-1}

# Range of the lengths of k. Empty means "use the k_min/k_max predefined in
# configuration_{organism}.ini", which is the point of this sweep.
K_MIN=${K_MIN:-1}
K_MAX=${K_MAX:-10}

# Destination of the results (kept in the home, so the job can be resumed)
OUTPUTS=${OUTPUTS:-$CASTOR/outputs}
CONFIGS_TMP=${CONFIGS_TMP:-$CASTOR/configs_tmp}
LOGS=${LOGS:-$CASTOR/logs}

# Rebuild the train/test split even if it already exists. Kept at 1 so every
# repetition gets an independent partition, which is what the standard deviation
# is measuring; set to 0 to repeat on a single frozen split instead.
RESPLIT=${RESPLIT:-1}
# Rerun the organism/k combinations already completed (FORCE=1)
FORCE=${FORCE:-0}
# Maximum duration of a single organism/k run, in seconds (empty = no limit)
RUN_TIMEOUT=${RUN_TIMEOUT:-}

# Organisms handled by this script (the name matches configuration_{organism}.ini)
ALL_ORGANISMS=(denv hbv hiv mkpx sars)

# Function to describe an organism: dataset directory, temporary directory and
# prefix of the concatenated train/test fasta files
describe_organism() {
	case "$1" in
		denv)      DATASET=$DATASETS/denv/data_clean;     TEMPDATASET=$TEMPROOT/denv;     PREFIX=denv ;;
		hbv)       DATASET=$DATASETS/HBV/data_clean;        TEMPDATASET=$TEMPROOT/HBV;        PREFIX=hbv ;;
		hiv)       DATASET=$DATASETS/hiv/data_clean;        TEMPDATASET=$TEMPROOT/hiv;        PREFIX=hiv ;;
		mkpx) DATASET=$DATASETS/mkpx/data_clean;       TEMPDATASET=$TEMPROOT/mkpx;       PREFIX=mkpx ;;
		sars)      DATASET=$DATASETS/sars_cov2/data_clean;  TEMPDATASET=$TEMPROOT/sars_cov2;  PREFIX=sars ;;
		*) echo "Unknown organism: $1" >&2; return 1 ;;
	esac
	# Directory holding the variant sub directories, split in train/ and test/
	SPLITROOT=$TEMPDATASET
	return 0
}

# Select the organisms to process: command line, then SGE array task, then all of them
if [ "$#" -gt 0 ]; then
	ORGANISMS=("$@")
elif [ -n "${SGE_TASK_ID:-}" ] && [ "${SGE_TASK_ID}" != "undefined" ]; then
	ORGANISMS=("${ALL_ORGANISMS[$((SGE_TASK_ID - 1))]}")
else
	ORGANISMS=("${ALL_ORGANISMS[@]}")
fi

mkdir -p "$OUTPUTS" "$CONFIGS_TMP" "$LOGS" "$TEMPROOT"

source $FEHOME/.venv/bin/activate

export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

REPORTS=${REPORTS:-$OUTPUTS/reports.csv}

# Pass the range of k only when it was explicitly given
K_RANGE=()
[ -n "$K_MIN" ] && K_RANGE+=(--k-min "$K_MIN")
[ -n "$K_MAX" ] && K_RANGE+=(--k-max "$K_MAX")

echo "[$(date '+%F %T')] $REPEATS repetitions over: ${ORGANISMS[*]} (k range: ${K_MIN:-config}..${K_MAX:-config})"

for REP in $(seq "$REPEAT_START" "$((REPEAT_START + REPEATS - 1))"); do
echo "########## [$(date '+%F %T')] repetition $REP ##########"

for ORGANISM in "${ORGANISMS[@]}"; do
	# Get the description of the organism
	if ! describe_organism "$ORGANISM"; then
		echo "[$(date '+%F %T')] $ORGANISM: unknown organism, skipped" >&2
		continue
	fi

	echo "[$(date '+%F %T')] $ORGANISM: staging $DATASET -> $TEMPDATASET"
	# Copy the dataset to the local scratch of the node (only once)
	if [ ! -d "$TEMPDATASET" ]; then
		cp -r "$DATASET" "$TEMPDATASET" || { echo "[$(date '+%F %T')] $ORGANISM: copy failed" >&2; continue; }
	fi

	TRAIN=$SPLITROOT/train/${PREFIX}_train.fasta
	TEST=$SPLITROOT/test/${PREFIX}_test.fasta

	# Build the balanced train/test split (kept between resumptions, so every
	# length of k is evaluated on exactly the same sequences)
	if [ "$RESPLIT" = "1" ] || [ ! -s "$TRAIN" ] || [ ! -s "$TEST" ]; then
		echo "[$(date '+%F %T')] $ORGANISM: building the balanced split in $SPLITROOT"
		if ! $BALANCEDATASET/testcl.sh "$SPLITROOT"; then
			echo "[$(date '+%F %T')] $ORGANISM: split failed, skipped" >&2
			echo -e "$(date '+%F %T')\t$ORGANISM\tk=-\tsplit failed" >> "$LOGS/errors.log"
			continue
		fi
		# Concatenate the per variant fasta files produced by the splitter
		( cd "$SPLITROOT/train" && cat *.fasta > "$TRAIN" ) || continue
		( cd "$SPLITROOT/test"  && cat *.fasta > "$TEST" )  || continue
	else
		echo "[$(date '+%F %T')] $ORGANISM: reusing the existing split ($TRAIN)"
	fi

	# Build the optional flags of the sweep
	EXTRA=()
	[ "$FORCE" = "1" ] && EXTRA+=(--force)
	[ -n "$RUN_TIMEOUT" ] && EXTRA+=(--timeout "$RUN_TIMEOUT")

	# Run the sweep of this organism: one temporary .ini and one output
	# directory per length of k, failures logged and never fatal
	python3 "$CASTOR/kmer_sweep.py" \
		--organisms "$ORGANISM" \
		--repeats 1 --repeat-start "$REP" \
		--config-dir "$CASTOR" \
		--configs-tmp "$CONFIGS_TMP" \
		--outputs "$OUTPUTS" \
		--logs "$LOGS" \
		--reports "$REPORTS" \
		--override "$ORGANISM:training_fasta=$TRAIN" \
		--override "$ORGANISM:testing_fasta=$TEST" \
		${K_RANGE[@]+"${K_RANGE[@]}"} \
		${EXTRA[@]+"${EXTRA[@]}"}

	echo "[$(date '+%F %T')] $ORGANISM: done"
done

# Refresh the mean/std summary after every repetition, so partial results stay
# usable if the job is killed before the last one
python3 "$CASTOR/aggregate_reports.py" --outputs "$OUTPUTS" --reports "$REPORTS"
done

echo "[$(date '+%F %T')] sweep finished, results in $OUTPUTS (errors in $LOGS/errors.log)"
echo "[$(date '+%F %T')] raw tables: $REPORTS | summary: $OUTPUTS/reports_summary.csv"
