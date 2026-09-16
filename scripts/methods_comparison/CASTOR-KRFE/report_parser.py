#!/usr/bin/env python3
"""
Parsing of the "Classification report" table printed at the end of a run.log.

The report looks like

	Classification report
	                 precision    recall  f1-score   support

	type1_reducted       0.97      0.56      0.71      2247
	...
	      accuracy                           0.50      4765
	     macro avg       0.76      0.46      0.42      4765
	  weighted avg       0.77      0.50      0.50      4765

A class row carries four numbers, the accuracy row only two (accuracy and
support) and the average rows four again. The label may contain spaces
("macro avg"), so the numbers are read from the end of the line.
"""

import os
import re

# Columns of the collected reports csv. k_length/n_k_mers are filled by
# CASTOR-KRFE, which picks the length of k itself inside its k_min..k_max range.
FIELDS = ("organism", "k", "repeat", "label", "kind",
          "precision", "recall", "f1_score", "support",
          "k_length", "n_k_mers", "seconds", "run_log")

# Scalars printed once per run, outside the classification report table. The
# LAST occurrence is kept: CASTOR-KRFE prints one of these per evaluated k while
# it sweeps its k_min..k_max range and the selected one last.
RUN_METRICS = {
	"k_length": r"Length of k\s*=\s*(\S+)",
	"n_k_mers": r"Number of k-mers\s*=\s*(\S+)",
	"seconds": r"^Time:\s*([\d.]+)",
}

# Evaluation of a single k, printed by the internal k sweep of CASTOR-KRFE
EVALUATION_METRICS = {
	"k": r"Length of k\s*=\s*(\d+)",
	"f1_score": r"Evaluation score \(F1 score\)\s*=\s*([\d.]+)",
	"n_k_mers": r"Number of k-mers\s*=\s*(\d+)",
}

NUMBER = re.compile(r"^\d+(?:\.\d+)?$")


def parseReportTable(lines):
	"""Convert the lines of a classification report table into rows."""
	rows = []
	for line in lines:
		if not line.strip(): continue
		tokens = line.split()
		# Header of the table and whatever follows the report
		if tokens[0] in ("precision", "Predictions", "KEVOLVE:"): continue
		numbers = []
		while tokens and NUMBER.match(tokens[-1]):
			numbers.append(float(tokens.pop()))
		numbers.reverse()
		label = " ".join(tokens)
		if not label or not numbers: continue
		if label == "accuracy" and len(numbers) == 2:
			rows.append({"label": "accuracy", "kind": "overall", "precision": "",
				"recall": "", "f1_score": numbers[0], "support": numbers[1]})
		elif len(numbers) == 4:
			rows.append({"label": label,
				"kind": "average" if label.endswith("avg") else "class",
				"precision": numbers[0], "recall": numbers[1],
				"f1_score": numbers[2], "support": numbers[3]})
	return rows


def parseEvaluations(path):
	"""Extract one row per k evaluated by the internal sweep of the method.

	The n-th length of k is paired with the n-th score and the n-th number of
	k-mers, in the order they appear in the log. A method that evaluates a
	single k (KEVOLVE) simply yields a single row, or none at all when it does
	not print these lines.
	"""
	if not os.path.isfile(path): return []
	with open(path, errors = "replace") as f:
		text = f.read()

	columns = {key: re.findall(pattern, text, flags = re.M)
		for key, pattern in EVALUATION_METRICS.items()}

	rows = []
	for values in zip(*columns.values()):
		row = dict(zip(columns, values))
		row["label"] = ""
		row["kind"] = "evaluation"
		row["run_log"] = path
		rows.append(row)
	return rows


def parseRunLog(path):
	"""Extract the classification report of a run.log (last one in the file)."""
	if not os.path.isfile(path): return []
	with open(path, errors = "replace") as f:
		text = f.read()

	start = text.rfind("Classification report")
	if start == -1: return []
	block = text[start:].splitlines()[1:]
	# The report ends where KEVOLVE reports where the predictions were saved
	end = len(block)
	for i, line in enumerate(block):
		if line.startswith("Predictions saved"):
			end = i
			break

	extra = {}
	for key, pattern in RUN_METRICS.items():
		matches = re.findall(pattern, text, flags = re.M)
		extra[key] = matches[-1] if matches else ""

	rows = parseReportTable(block[:end])
	for row in rows:
		row.update(extra)
		row["run_log"] = path
	return rows


# Self check: python3 report_parser.py
if __name__ == "__main__":
	import tempfile
	SAMPLE = """Length of k = 5
Number of k-mers = 10
Evaluation score (F1 score) = 0.80

Length of k = 6
Number of k-mers = 12
Evaluation score (F1 score) = 0.91

Time:  120

Classification report 
                 precision    recall  f1-score   support

type1_reducted       0.97      0.56      0.71      2247
type2_reducted       0.76      0.13      0.22      1432

      accuracy                           0.50      4765
     macro avg       0.76      0.46      0.42      4765
  weighted avg       0.77      0.50      0.50      4765

Predictions saved at the path: /x
"""
	with tempfile.NamedTemporaryFile("w", suffix = ".log", delete = False) as f:
		f.write(SAMPLE)
		path = f.name

	rows = parseRunLog(path)
	assert len(rows) == 5, rows
	assert [row["label"] for row in rows] == ["type1_reducted", "type2_reducted",
		"accuracy", "macro avg", "weighted avg"]
	# Label with a space, and the accuracy row with only two numbers
	assert rows[3]["precision"] == 0.76 and rows[3]["kind"] == "average"
	assert rows[2]["f1_score"] == 0.50 and rows[2]["precision"] == ""
	# Run scalars keep the LAST occurrence (the k the method selected)
	assert rows[0]["k_length"] == "6" and rows[0]["n_k_mers"] == "12"
	assert rows[0]["seconds"] == "120"

	evaluations = parseEvaluations(path)
	assert [(row["k"], row["f1_score"]) for row in evaluations] == [("5", "0.80"), ("6", "0.91")]
	assert all(row["kind"] == "evaluation" for row in evaluations)

	# A log without a report yields nothing instead of raising
	with tempfile.NamedTemporaryFile("w", suffix = ".log", delete = False) as f:
		f.write("crashed\n")
		empty = f.name
	assert parseRunLog(empty) == [] and parseEvaluations(empty) == []
	assert parseRunLog("/does/not/exist") == []

	os.remove(path)
	os.remove(empty)
	print("report_parser: all checks passed")
