#!/usr/bin/env python3
"""
Aggregation of the classification reports collected over the sweep repetitions.

Reads every run.log under the outputs tree (and/or an existing reports.csv
written by run_kmer_sweep.py) and produces two files:

  reports.csv          one row per organism/k/repetition/label - the raw tables
  reports_summary.csv  one row per organism/k/label - n, mean and standard
                       deviation of precision, recall, f1-score and runtime

The standard deviation is the sample one (ddof = 1), which is what a mean over
repetitions calls for; with a single repetition it is left empty.

Usage:
	python3 aggregate_reports.py --outputs /path/to/outputs
"""

import argparse
import csv
import glob
import os
import statistics
import sys

import report_parser

# Metrics aggregated over the repetitions
METRICS = ("precision", "recall", "f1_score", "seconds")

SUMMARY_FIELDS = (["organism", "k", "label", "kind", "n", "support"] +
	[metric + suffix for metric in METRICS for suffix in ("_mean", "_std")])


def runIdentity(path, outputs):
	"""Infer organism, k and repetition from the path of a run.log.

	Handles both layouts: outputs/{organism}/k{k}/rep{r}/run.log and the flat
	outputs/{organism}/k{k}/run.log, which is treated as repetition 1.
	"""
	parts = os.path.relpath(os.path.dirname(path), outputs).split(os.sep)
	repeat = 1
	if parts and parts[-1].startswith("rep"):
		try: repeat = int(parts[-1][3:])
		except ValueError: pass
		parts = parts[:-1]
	k = parts[-1][1:] if parts and parts[-1].startswith("k") else ""
	organism = parts[-2] if len(parts) >= 2 else (parts[0] if parts else "")
	return organism, k, repeat


def collectFromRunLogs(outputs):
	"""Parse every run.log below the outputs directory."""
	rows = []
	for path in sorted(glob.glob(os.path.join(outputs, "**", "run.log"), recursive = True)):
		organism, k, repeat = runIdentity(path, outputs)
		for row in report_parser.parseRunLog(path):
			row["organism"], row["k"], row["repeat"] = organism, k, repeat
			rows.append(row)
	return rows


def readReportsFile(path):
	"""Read an existing reports.csv, if there is one."""
	if not os.path.isfile(path): return []
	with open(path, newline = "") as f:
		return list(csv.DictReader(f))


def deduplicate(rows):
	"""Keep one row per organism/k/repetition/label (the last one wins)."""
	unique = {}
	for row in rows:
		unique[(str(row.get("organism", "")), str(row.get("k", "")),
			str(row.get("repeat", "")), str(row.get("label", "")))] = row
	return [unique[key] for key in sorted(unique)]


def number(value):
	"""Parse a csv cell as a float, or None when it is empty/not a number."""
	try: return float(value)
	except (TypeError, ValueError): return None


def summarise(rows):
	"""Mean and standard deviation of every metric over the repetitions."""
	groups = {}
	for row in rows:
		key = (str(row.get("organism", "")), str(row.get("k", "")), str(row.get("label", "")))
		groups.setdefault(key, []).append(row)

	summary = []
	for key in sorted(groups, key = lambda key: (key[0], int(key[1]) if key[1].isdigit() else 0, key[2])):
		organism, k, label = key
		group = groups[key]
		entry = {"organism": organism, "k": k, "label": label,
			"kind": group[0].get("kind", ""), "n": len(group),
			"support": group[0].get("support", "")}
		for metric in METRICS:
			values = [value for value in (number(row.get(metric)) for row in group) if value is not None]
			entry[metric + "_mean"] = round(statistics.fmean(values), 4) if values else ""
			entry[metric + "_std"] = round(statistics.stdev(values), 4) if len(values) > 1 else ""
		summary.append(entry)
	return summary


def write(path, fieldnames, rows):
	os.makedirs(os.path.dirname(path) or ".", exist_ok = True)
	with open(path, "w", newline = "") as f:
		writer = csv.DictWriter(f, fieldnames = list(fieldnames))
		writer.writeheader()
		for row in rows:
			writer.writerow({field: row.get(field, "") for field in fieldnames})


def main():
	argumentParser = argparse.ArgumentParser(
		description = "Aggregate the classification reports of the k-mer sweep repetitions")
	argumentParser.add_argument("--outputs", default = "outputs",
		help = "root directory of the outputs (default: ./outputs)")
	argumentParser.add_argument("--reports", default = None,
		help = "collected reports csv (default: {outputs}/reports.csv)")
	argumentParser.add_argument("--summary", default = None,
		help = "summary csv to write (default: {outputs}/reports_summary.csv)")
	argumentParser.add_argument("--no-rescan", action = "store_true",
		help = "summarise the existing reports csv without re-reading the run.log files")
	arguments = argumentParser.parse_args()

	reportsPath = arguments.reports or os.path.join(arguments.outputs, "reports.csv")
	summaryPath = arguments.summary or os.path.join(arguments.outputs, "reports_summary.csv")

	rows = readReportsFile(reportsPath)
	if not arguments.no_rescan:
		rows += collectFromRunLogs(arguments.outputs)
	rows = deduplicate(rows)

	if not rows:
		print("No classification report found under " + arguments.outputs, file = sys.stderr)
		return 1

	if not arguments.no_rescan:
		write(reportsPath, report_parser.FIELDS, rows)
	summary = summarise(rows)
	write(summaryPath, SUMMARY_FIELDS, summary)

	repetitions = len({(row.get("organism"), row.get("k"), row.get("repeat")) for row in rows})
	print(str(len(rows)) + " report rows from " + str(repetitions) + " runs")
	print("Raw tables : " + reportsPath)
	print("Summary    : " + summaryPath)

	worst = min((entry["n"] for entry in summary), default = 0)
	if worst < 2:
		print("Note: at least one organism/k/label has a single repetition, "
			"so its standard deviation is empty.")
	return 0


if __name__ == "__main__":
	sys.exit(main())
