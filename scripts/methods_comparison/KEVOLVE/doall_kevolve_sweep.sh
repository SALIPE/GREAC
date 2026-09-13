#!/bin/bash
#$ -S /bin/bash
#$ -N kevolve_sweep
#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs
#$ -cwd

# Automated k-mer length sweep (k = K_MIN..K_MAX) for several organisms,
# repeated REPEATS times so mean and standard deviation can be measured.
#
# Submit with:  qsub /home/a61491/KEVOLVE/doall_kevolve_sweep.sh
# Options through the environment, e.g.:
#   qsub -v ORGANISMS="sars hbv",K_MIN=4,K_MAX=8,REPEATS=100 doall_kevolve_sweep.sh
#
# The repetition is the OUTER loop: each repetition re-splits the datasets, so
# a repetition is an independent train/test partition and an independent run of
# the stochastic solution search, while every k inside one repetition still
# shares exactly the same partition (which is what makes the k values
# comparable). The classification report of every run is appended to
# $OUTPUTS/reports.csv, and aggregate_reports.py turns those tables into mean
# and standard deviation.
#
# Outputs land in /home/a61491/KEVOLVE/outputs/{organism}/k{k}/rep{r}/ .
#
# WARNING: the cost is REPEATS x organisms x (K_MAX - K_MIN + 1) runs. With the
# defaults that is 100 x 9 runs per organism, and a single large k already takes
# tens of minutes - size REPEATS/K_MAX to the time you actually have, or split
# the repetitions over several jobs with REPEAT_START/REPEATS.

KEVOLVE=/home/a61491/KEVOLVE
TEMPROOT=/tmp2/felipe
TEMPDIR=$TEMPROOT/kevolve
DATASETS=/home/a61491/datasets/original
BALANCEDATASET=/home/a61491/Fasta-splitter/FastaSplitter

#ORGANISMS=${ORGANISMS:-"sars denv hbv hiv mkpx"}
ORGANISMS=${ORGANISMS:-"hbv"}
K_MIN=${K_MIN:-2}
K_MAX=${K_MAX:-10}
# Number of repetitions of the whole sweep, and index of the first one
# (REPEAT_START lets a second job continue where a first one stopped)
REPEATS=${REPEATS:-100}
REPEAT_START=${REPEAT_START:-1}
OUTPUTS=${OUTPUTS:-$KEVOLVE/outputs}
REPORTS=${REPORTS:-$OUTPUTS/reports.csv}
# Per organism/k timeout in seconds (0 = no limit)
TIMEOUT=${TIMEOUT:-0}
# Set FORCE=1 to reprocess organism/k combinations already completed
FORCE_FLAG=""
[ "${FORCE:-0}" = "1" ] && FORCE_FLAG="--force"

mkdir -p "$TEMPROOT" "$TEMPDIR" "$OUTPUTS" "$KEVOLVE/logs"

source /home/a61491/.venv/bin/activate

# Working directory under /tmp2 holding the split of one organism.
temp_dataset_dir() {
	case $1 in
		sars) echo "$TEMPROOT/sars_cov2" ;;
		denv) echo "$TEMPROOT/dengue" ;;
		*)    echo "$TEMPROOT/$1" ;;
	esac
}

# Prepare the dataset of one organism: copy it to /tmp2, split it in balanced
# train/test sets and build the concatenated fasta files used by the .ini files.
# The copy is removed first so that every repetition gets a genuinely new
# random train/test partition instead of reusing the previous split.
prepare_dataset() {
	organism=$1
	rm -rf "$(temp_dataset_dir "$organism")"
	case $organism in
		sars)
			[ -d "$TEMPROOT/sars_cov2" ] || cp -r "$DATASETS/sars_cov2/data" "$TEMPROOT/sars_cov2"
			"$BALANCEDATASET/testcl_sars.sh" || return 1
			rm -f "$TEMPROOT/sars_cov2/train/sars_train.fasta" "$TEMPROOT/sars_cov2/test/sars_test.fasta"
			cat "$TEMPROOT"/sars_cov2/train/*.fasta > "$TEMPROOT/sars_cov2/train/sars_train.fasta"
			cat "$TEMPROOT"/sars_cov2/test/*.fasta  > "$TEMPROOT/sars_cov2/test/sars_test.fasta"
			;;
		denv)
			[ -d "$TEMPROOT/dengue" ] || cp -r "$DATASETS/dengue/data" "$TEMPROOT/dengue"
			"$BALANCEDATASET/testcl.sh" "$TEMPROOT/dengue" || return 1
			rm -f "$TEMPROOT/dengue/train/denv_train.fasta" "$TEMPROOT/dengue/test/denv_test.fasta"
			cat "$TEMPROOT"/dengue/train/*.fasta > "$TEMPROOT/dengue/train/denv_train.fasta"
			cat "$TEMPROOT"/dengue/test/*.fasta  > "$TEMPROOT/dengue/test/denv_test.fasta"
			;;
		hbv)
			[ -d "$TEMPROOT/hbv" ] || cp -r "$DATASETS/HBV/data" "$TEMPROOT/hbv"
			"$BALANCEDATASET/testcl.sh" "$TEMPROOT/hbv" || return 1
			rm -f "$TEMPROOT/hbv/train/hbv_train.fasta" "$TEMPROOT/hbv/test/hbv_test.fasta"
			cat "$TEMPROOT"/hbv/train/*.fasta > "$TEMPROOT/hbv/train/hbv_train.fasta"
			cat "$TEMPROOT"/hbv/test/*.fasta  > "$TEMPROOT/hbv/test/hbv_test.fasta"
			;;
		hiv)
			[ -d "$TEMPROOT/hiv" ] || cp -r "$DATASETS/hiv/data" "$TEMPROOT/hiv"
			"$BALANCEDATASET/testcl.sh" "$TEMPROOT/hiv" || return 1
			rm -f "$TEMPROOT/hiv/train/hiv_train.fasta" "$TEMPROOT/hiv/test/hiv_test.fasta"
			cat "$TEMPROOT"/hiv/train/*.fasta > "$TEMPROOT/hiv/train/hiv_train.fasta"
			cat "$TEMPROOT"/hiv/test/*.fasta  > "$TEMPROOT/hiv/test/hiv_test.fasta"
			;;
		mkpx)
			[ -d "$TEMPROOT/mkpx" ] || cp -r "$DATASETS/mkpx/data" "$TEMPROOT/mkpx"
			"$BALANCEDATASET/testcl.sh" "$TEMPROOT/mkpx" || return 1
			rm -f "$TEMPROOT/mkpx/train/mkpx_train.fasta" "$TEMPROOT/mkpx/test/mkpx_test.fasta"
			cat "$TEMPROOT"/mkpx/train/*.fasta > "$TEMPROOT/mkpx/train/mkpx_train.fasta"
			cat "$TEMPROOT"/mkpx/test/*.fasta  > "$TEMPROOT/mkpx/test/mkpx_test.fasta"
			;;
		*)
			echo "Unknown organism: $organism" >&2
			return 1
			;;
	esac
}

REPEAT_END=$((REPEAT_START + REPEATS - 1))

for rep in $(seq "$REPEAT_START" "$REPEAT_END"); do
	echo "########## [$(date '+%F %T')] repetition $rep / $REPEAT_END ##########"

	for organism in $ORGANISMS; do
		echo "=== [$(date '+%F %T')] rep $rep | preparing dataset: $organism ==="
		if ! prepare_dataset "$organism"; then
			echo "$(date '+%F %T')	organism=$organism	k=-	rep=$rep	dataset preparation failed" >> "$KEVOLVE/logs/errors.log"
			echo "!!! dataset preparation failed for $organism, skipping it"
			continue
		fi

		echo "=== [$(date '+%F %T')] rep $rep | k-mer sweep: $organism (k=$K_MIN..$K_MAX) ==="
		python3 -u "$KEVOLVE/run_kmer_sweep.py" \
			--organisms "$organism" \
			--k-min "$K_MIN" --k-max "$K_MAX" \
			--outputs "$OUTPUTS" \
			--configs "$KEVOLVE/configs_tmp" \
			--logs "$KEVOLVE/logs" \
			--timeout "$TIMEOUT" \
			--repeat "$rep" \
			--reports "$REPORTS" \
			$FORCE_FLAG
	done

	# Refresh the mean/std summary after every repetition, so partial results
	# are usable even if the job is killed before the last one.
	python3 -u "$KEVOLVE/aggregate_reports.py" --outputs "$OUTPUTS" --reports "$REPORTS"
done

echo "=== [$(date '+%F %T')] sweep finished, outputs in $OUTPUTS ==="
echo "=== raw tables: $REPORTS | summary: $OUTPUTS/reports_summary.csv ==="
# The datasets in /tmp2 are kept so an interrupted job can be resumed cheaply.
# Uncomment to clean the node:
rm -r $TEMPROOT
