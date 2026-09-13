#!/usr/bin/env python3
"""
Automated k-mer length sweep for KEVOLVE.

For each organism and each k in [K_MIN, K_MAX]:
  - a temporary copy of the original .ini is written with only the output paths
    and k changed (the original .ini files are never modified);
  - the analysis is executed with that configuration;
  - every output is written under the run directory
    (k_mers/, models/, predictions/, analysis/, run.log).

The sweep is repeated several times so that mean and standard deviation can be
measured over the repetitions. With --repeat R the run directory becomes
outputs/{organism}/k{k}/rep{R}/ and the configuration configs_tmp/{organism}/
rep{R}/k{k}.ini; without it the flat layout outputs/{organism}/k{k}/ is kept.
Each repetition is driven by the shell wrapper, which re-splits the dataset
first, so a repetition is an independent train/test partition (and an
independent run of the stochastic solution search) while every k inside one
repetition still shares exactly the same partition.

After each run the classification report found in run.log is parsed and appended
to outputs/reports.csv (one row per class/average), which is what
aggregate_reports.py summarises into mean and standard deviation.

Failures of one organism/k combination are logged to logs/errors.log and the
loop continues. Completed combinations are marked with a .done file, so the
script can be interrupted and resumed without reprocessing them.
"""

import argparse
import configparser
import csv
import datetime
import os
import shutil
import subprocess
import sys
import time

import report_parser

# Root of the KEVOLVE installation (directory holding this script)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Organism -> original configuration file (never modified)
ORGANISMS = {
	"sars": os.path.join(BASE_DIR, "configuration_file_sars.ini"),
	"denv": os.path.join(BASE_DIR, "configuration_file_denv.ini"),
	"hbv":  os.path.join(BASE_DIR, "configuration_file_hbv.ini"),
	"hiv":  os.path.join(BASE_DIR, "configuration_file_hiv.ini"),
	"mkpx": os.path.join(BASE_DIR, "configuration_file_mkpx.ini"),
}

# Sub-directories created for each organism/k run
RUN_SUBDIRS = ("k_mers", "models", "predictions", "analysis")


def runDirectory(outputs, organism, k, repeat):
	"""Run directory of one organism/k/repetition (flat layout when repeat is None)."""
	runDir = os.path.join(outputs, organism, "k" + str(k))
	if repeat is not None: runDir = os.path.join(runDir, "rep" + str(repeat))
	return runDir


def configurationPath(configs, organism, k, repeat):
	"""Temporary .ini path of one organism/k/repetition."""
	directory = os.path.join(configs, organism)
	if repeat is not None: directory = os.path.join(directory, "rep" + str(repeat))
	return os.path.join(directory, "k" + str(k) + ".ini")


def saveReport(reportsPath, organism, k, repeat, runDir):
	"""Append the classification report of one run to the collected reports file.

	One row per line of the report table (each class plus accuracy and the
	macro/weighted averages), so the repetitions can be aggregated afterwards.
	"""
	rows = report_parser.parseRunLog(os.path.join(runDir, "run.log"))
	if not rows:
		return 0
	os.makedirs(os.path.dirname(reportsPath) or ".", exist_ok = True)
	isNew = not os.path.isfile(reportsPath)
	with open(reportsPath, "a", newline = "") as f:
		writer = csv.DictWriter(f, fieldnames = report_parser.FIELDS)
		if isNew: writer.writeheader()
		for row in rows:
			row["organism"] = organism
			row["k"] = k
			row["repeat"] = repeat if repeat is not None else 1
			writer.writerow({field: row.get(field, "") for field in report_parser.FIELDS})
	return len(rows)


def log(message):
	"""Print a timestamped progress message (flushed for SGE log files)."""
	stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print("[" + stamp + "] " + message, flush = True)


def logError(errorLogPath, organism, k, message):
	"""Append an error entry to the error log file."""
	stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	os.makedirs(os.path.dirname(errorLogPath), exist_ok = True)
	with open(errorLogPath, "a") as f:
		f.write(stamp + "\torganism=" + organism + "\tk=" + str(k) + "\t" + message.replace("\n", " ") + "\n")


def prepareRunDirectory(runDir, force):
	"""Create a clean output directory tree for one organism/k run.

	The directories must be empty: ml.fit()/ml.predict() iterate over every
	file found in k_mers_path and model_path, so leftovers from a previous
	attempt would corrupt the results.
	"""
	if os.path.isdir(runDir) and force:
		shutil.rmtree(runDir)
	for subdir in RUN_SUBDIRS:
		path = os.path.join(runDir, subdir)
		if os.path.isdir(path): shutil.rmtree(path)
		os.makedirs(path)


def writeConfiguration(originalConfig, configPath, k, runDir):
	"""Write a copy of the original .ini with k and the output paths updated."""
	parser = configparser.ConfigParser()
	# Keep the comments/keys untouched apart from the values we set below
	if not parser.read(originalConfig):
		raise IOError("configuration file not readable: " + originalConfig)
	parser.set("parameters", "k", str(k))
	parser.set("parameters", "k_mers_path", os.path.join(runDir, "k_mers"))
	parser.set("parameters", "model_path", os.path.join(runDir, "models"))
	parser.set("parameters", "prediction_path", os.path.join(runDir, "predictions"))
	parser.set("parameters", "analysis_report_path", os.path.join(runDir, "analysis"))
	parser.set("parameters", "k_mers_to_analyze_path", os.path.join(runDir, "k_mers", "1.fasta"))
	os.makedirs(os.path.dirname(configPath), exist_ok = True)
	with open(configPath, "w") as f:
		parser.write(f)


def runCombination(organism, k, configPath, runDir, timeout):
	"""Execute one organism/k analysis, streaming its output to run.log."""
	logPath = os.path.join(runDir, "run.log")
	command = [sys.executable, "-u", os.path.join(BASE_DIR, "evaluate_kmer.py"), configPath]
	with open(logPath, "w") as logFile:
		logFile.write("command: " + " ".join(command) + "\n\n")
		logFile.flush()
		process = subprocess.Popen(command, cwd = BASE_DIR, stdout = logFile, stderr = subprocess.STDOUT)
		try:
			returnCode = process.wait(timeout = timeout)
		except subprocess.TimeoutExpired:
			process.kill()
			process.wait()
			raise RuntimeError("timeout after " + str(timeout) + "s (see " + logPath + ")")
	if returnCode != 0:
		raise RuntimeError("exit code " + str(returnCode) + " (see " + logPath + ")")


def tail(path, n = 15):
	"""Return the last n lines of a file, for the error log."""
	try:
		with open(path) as f: return " | ".join(f.read().splitlines()[-n:])
	except IOError:
		return ""


def main():
	argumentParser = argparse.ArgumentParser(description = "KEVOLVE k-mer length sweep")
	argumentParser.add_argument("--organisms", nargs = "+", default = sorted(ORGANISMS),
		choices = sorted(ORGANISMS), help = "organisms to process (default: all)")
	argumentParser.add_argument("--k-min", type = int, default = 1, help = "first k value (default: 1)")
	argumentParser.add_argument("--k-max", type = int, default = 10, help = "last k value (default: 10)")
	argumentParser.add_argument("--outputs", default = os.path.join(BASE_DIR, "outputs"),
		help = "root directory of the outputs (default: ./outputs)")
	argumentParser.add_argument("--configs", default = os.path.join(BASE_DIR, "configs_tmp"),
		help = "root directory of the temporary configuration files (default: ./configs_tmp)")
	argumentParser.add_argument("--logs", default = os.path.join(BASE_DIR, "logs"),
		help = "directory of the error log (default: ./logs)")
	argumentParser.add_argument("--timeout", type = int, default = 0,
		help = "per run timeout in seconds, 0 to disable (default: 0)")
	argumentParser.add_argument("--repeat", type = int, default = None,
		help = "index of the repetition; outputs go to outputs/{organism}/k{k}/rep{R} "
		       "(default: none, flat outputs/{organism}/k{k})")
	argumentParser.add_argument("--reports", default = None,
		help = "csv file collecting the classification reports of every run "
		       "(default: {outputs}/reports.csv)")
	argumentParser.add_argument("--force", action = "store_true",
		help = "reprocess organism/k combinations already completed")
	arguments = argumentParser.parse_args()

	errorLogPath = os.path.join(arguments.logs, "errors.log")
	timeout = arguments.timeout if arguments.timeout > 0 else None
	reportsPath = arguments.reports or os.path.join(arguments.outputs, "reports.csv")
	repeat = arguments.repeat

	nSuccess, nFailure, nSkipped = 0, 0, 0

	for organism in arguments.organisms:
		originalConfig = ORGANISMS[organism]
		if not os.path.isfile(originalConfig):
			logError(errorLogPath, organism, "-", "missing configuration file: " + originalConfig)
			log("ORGANISM " + organism + ": FAILED (missing configuration file " + originalConfig + ")")
			nFailure += 1
			continue

		for k in range(arguments.k_min, arguments.k_max + 1):
			runDir = runDirectory(arguments.outputs, organism, k, repeat)
			doneMarker = os.path.join(runDir, ".done")
			configPath = configurationPath(arguments.configs, organism, k, repeat)
			label = organism + " | k=" + str(k)
			if repeat is not None: label += " | rep=" + str(repeat)

			if os.path.isfile(doneMarker) and not arguments.force:
				log(label + " | SKIPPED (already completed)")
				nSkipped += 1
				continue

			log(label + " | START")
			start = time.time()
			try:
				prepareRunDirectory(runDir, arguments.force)
				writeConfiguration(originalConfig, configPath, k, runDir)
				runCombination(organism, k, configPath, runDir, timeout)
			except Exception as exception:
				nFailure += 1
				message = str(exception) + " || " + tail(os.path.join(runDir, "run.log"))
				logError(errorLogPath, organism, k, message)
				log(label + " | FAILED after " +
					str(round(time.time() - start, 1)) + "s: " + str(exception))
				# Continue with the next combination without stopping the loop
				continue
			duration = round(time.time() - start, 1)

			# Collect the classification report before moving on: a failure here
			# must not invalidate a run that actually completed.
			try:
				nRows = saveReport(reportsPath, organism, k, repeat, runDir)
				if not nRows:
					logError(errorLogPath, organism, k, "no classification report found in run.log")
					log(label + " | WARNING: no classification report found in run.log")
			except Exception as exception:
				logError(errorLogPath, organism, k, "report not collected: " + str(exception))
				log(label + " | WARNING: report not collected: " + str(exception))

			with open(doneMarker, "w") as f:
				f.write(datetime.datetime.now().isoformat() + "\tduration_s=" + str(duration) + "\n")
			nSuccess += 1
			log(label + " | SUCCESS in " + str(duration) + "s -> " + runDir)

	log("Sweep finished: " + str(nSuccess) + " success, " + str(nFailure) +
		" failure, " + str(nSkipped) + " skipped")
	log("Classification reports collected in: " + reportsPath)
	if nFailure: log("Errors logged in: " + errorLogPath)
	return 0


if __name__ == "__main__":
	sys.exit(main())
