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

# Columns of the collected reports csv
FIELDS = ("organism", "k", "repeat", "label", "kind",
          "precision", "recall", "f1_score", "support", "seconds", "run_log")

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

	match = re.search(r"^Time:\s*([\d.]+)", text, flags = re.M)
	seconds = float(match.group(1)) if match else ""

	rows = parseReportTable(block[:end])
	for row in rows:
		row["seconds"] = seconds
		row["run_log"] = path
	return rows
