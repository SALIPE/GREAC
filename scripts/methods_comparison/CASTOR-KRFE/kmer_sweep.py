# Imports
import os
import re
import csv
import sys
import json
import time
import shutil
import logging
import argparse
import datetime
import subprocess
import configparser

import report_parser

# Root directory of the project (this file location)
ROOT = os.path.dirname(os.path.abspath(__file__))

# Section holding the parameters inside the configuration files
SECTION = "parameters"

# Fields holding the range of k. CASTOR-KRFE selects the best length of k
# itself inside k_min..k_max, so the range of the original configuration is
# kept as it is and only overridden when --k-min/--k-max are given.
K_FIELDS = ("k_min", "k_max")

# Output fields to redirect, mapped to their destination inside the run directory.
# A destination ending with a separator is created as a directory, otherwise as a file.
OUTPUT_FIELDS = {
	"k_mers_path": "k_mers.fasta",
	"model_path": "model.pkl",
	"prediction_path": "prediction.csv",
	"analysis_report_path": "analysis/",
	"k_mers_to_analyze_path": "k_mers/1.fasta",
}


# Function to discover the organisms owning a configuration file
def discoverOrganisms(config_dir):
	# Initialize the organism dictionary
	organisms = dict()
	# Iterate through the configuration files of the directory
	for name in sorted(os.listdir(config_dir)):
		# Keep only the configuration_{organism}.ini files
		match = re.fullmatch(r"configuration_(.+)\.ini", name)
		# Skip anything else (configuration_file.ini is the template, not an organism)
		if not match or match.group(1) == "file": continue
		# Save the organism and the path of its configuration file
		organisms[match.group(1)] = os.path.join(config_dir, name)
	# Return the organism dictionary
	return organisms


# Function to parse the "organism:field=value" overrides given on the command line
def parseOverrides(values):
	# Initialize the override dictionary
	overrides = dict()
	# Iterate through the given overrides
	for value in values or []:
		# Split the organism from the assignment
		if ":" not in value or "=" not in value:
			raise ValueError("Invalid override (expected organism:field=value): " + value)
		organism, assignment = value.split(":", 1)
		field, field_value = assignment.split("=", 1)
		# Save the override for this organism
		overrides.setdefault(organism.strip(), dict())[field.strip()] = field_value.strip()
	# Return the override dictionary
	return overrides


# Function to redirect an output field to the run directory
def redirectOutput(field, original_value, run_directory):
	# Get the destination associated to this field
	destination = OUTPUT_FIELDS[field]
	# If the original value is a directory, keep a directory as destination
	if not os.path.splitext(original_value)[1]: destination = os.path.splitext(destination)[0] + "/"
	# Build the absolute destination path
	path = os.path.join(run_directory, destination)
	# Create the destination directory (the directory itself or the parent of the file)
	os.makedirs(path if destination.endswith("/") else os.path.dirname(path), exist_ok = True)
	# Return the destination path without the trailing separator
	return path.rstrip("/")


# Function to write the temporary configuration file of one run
def writeConfiguration(source_config, destination_config, k_range, run_directory, overrides):
	# Initialize the parser, keeping the case of the fields
	parser = configparser.ConfigParser()
	parser.optionxform = str
	# Read the original configuration file (never modified)
	if not parser.read(source_config):
		raise FileNotFoundError("Unreadable configuration file: " + source_config)
	# The parameter section is mandatory
	if not parser.has_section(SECTION):
		raise KeyError("Missing section [" + SECTION + "] in " + source_config)
	# The configuration must define the range of k explored by the method
	for field in K_FIELDS:
		if not parser.has_option(SECTION, field):
			raise KeyError("Missing field " + field + " in " + source_config)
	# Override the range only when the command line asks for it, so by default
	# the method runs on its own predefined k_min..k_max
	for field, value in zip(K_FIELDS, k_range):
		if value is not None: parser.set(SECTION, field, str(value))
	# Redirect every output field to the run directory
	for field in OUTPUT_FIELDS:
		if parser.has_option(SECTION, field):
			parser.set(SECTION, field, redirectOutput(field, parser.get(SECTION, field), run_directory))
	# Apply the user overrides (dataset paths staged by the job script, ...)
	for field, value in overrides.items(): parser.set(SECTION, field, value)
	# Save the temporary configuration file
	os.makedirs(os.path.dirname(destination_config), exist_ok = True)
	with open(destination_config, "w") as f: parser.write(f)
	# Return the path of the temporary configuration file
	return destination_config


# Function to extract the summary metrics from the log of a run
def parseMetrics(log_path):
	# Initialize the metrics
	metrics = {"f1_score": "", "k_length": "", "n_k_mers": ""}
	# Patterns displayed by the extraction step
	patterns = {
		"f1_score": r"Evaluation score \(F1 score\) =\s*(\S+)",
		"k_length": r"Length of k =\s*(\S+)",
		"n_k_mers": r"Number of k-mers =\s*(\S+)",
	}
	# Read the log of the run
	try:
		with open(log_path, errors = "replace") as f: content = f.read()
	except OSError: return metrics
	# Search each pattern inside the log
	for key, pattern in patterns.items():
		match = re.search(pattern, content)
		if match: metrics[key] = match.group(1)
	# Return the metrics
	return metrics


# Function to append the classification report of a run to the collected reports.
# Same file and same columns as the KEVOLVE sweep, so both methods are compared
# with exactly the same metrics (per class precision/recall/f1, accuracy and the
# macro/weighted averages, plus the length of k the method selected itself).
def saveReport(reports_path, organism, repeat, log_path):
	# Parse the classification report printed at the end of the log, plus one row
	# per k evaluated by the internal sweep of the method (kind = "evaluation"),
	# so the k sweep itself gets a mean and a standard deviation too
	rows = report_parser.parseRunLog(log_path) + report_parser.parseEvaluations(log_path)
	# Nothing to save when the run did not print anything parsable
	if not rows: return 0
	# Create the destination directory and remember if the header is needed
	os.makedirs(os.path.dirname(reports_path) or ".", exist_ok = True)
	is_new = not os.path.isfile(reports_path)
	# Append one row per line of the report table
	with open(reports_path, "a", newline = "") as f:
		writer = csv.DictWriter(f, fieldnames = list(report_parser.FIELDS))
		if is_new: writer.writeheader()
		for row in rows:
			row["organism"] = organism
			row["repeat"] = repeat
			# An evaluation row carries the k it evaluated; on a classification
			# report row k stays empty, the selected k being reported in k_length
			row.setdefault("k", "")
			writer.writerow({field: row.get(field, "") for field in report_parser.FIELDS})
	# Return the number of saved rows
	return len(rows)


# Function to run a single organism/k combination
def runCombination(run_directory, config_file, timeout):
	# Create the run directory
	os.makedirs(run_directory, exist_ok = True)
	# Path of the log of this run
	log_path = os.path.join(run_directory, "run.log")
	# Build the environment of the run (non interactive matplotlib backend)
	env = dict(os.environ)
	env["MPLBACKEND"] = "Agg"
	env["PYTHONUNBUFFERED"] = "1"
	# Build the command to execute
	command = [sys.executable, os.path.join(ROOT, "run_analysis.py"), config_file]
	# Start the chronometer
	started = time.time()
	# Run the analysis in a subprocess, so a crash never interrupts the sweep
	with open(log_path, "w") as log:
		log.write("# command: " + " ".join(command) + "\n")
		log.write("# configuration: " + config_file + "\n")
		log.flush()
		process = subprocess.run(command, cwd = ROOT, env = env, stdout = log,
			stderr = subprocess.STDOUT, timeout = timeout)
	# Compute the duration of the run
	duration = time.time() - started
	# Return the result of the run
	return process.returncode, duration, log_path


# Function to run the whole sweep
def sweep(args):
	# Discover the organisms to process
	available = discoverOrganisms(args.config_dir)
	organisms = args.organisms or sorted(available)
	# Parse the configuration overrides
	overrides = parseOverrides(args.override)
	# Build the absolute output directories
	configs_tmp = os.path.abspath(args.configs_tmp)
	outputs = os.path.abspath(args.outputs)
	logs = os.path.abspath(args.logs)
	os.makedirs(logs, exist_ok = True)
	# Configure the progress logging (console) and the error logging (file)
	logging.basicConfig(level = logging.INFO, format = "%(asctime)s [%(levelname)s] %(message)s",
		datefmt = "%Y-%m-%d %H:%M:%S", stream = sys.stdout)
	error_log = logging.getLogger("errors")
	error_log.propagate = False
	handler = logging.FileHandler(os.path.join(logs, "errors.log"))
	handler.setFormatter(logging.Formatter("%(asctime)s\t%(message)s", datefmt = "%Y-%m-%d %H:%M:%S"))
	error_log.addHandler(handler)
	# Path of the consolidated summary and of the collected classification reports
	summary_path = os.path.join(outputs, "summary.csv")
	reports_path = os.path.abspath(args.reports) if args.reports else os.path.join(outputs, "reports.csv")
	os.makedirs(outputs, exist_ok = True)
	# Initialize the counters
	succeeded, failed, skipped = 0, 0, 0
	# Iterate through the organisms
	for organism in organisms:
		# Skip the organisms without a configuration file
		if organism not in available:
			logging.error("Organism %s has no configuration file in %s, skipped", organism, args.config_dir)
			continue
		# One run of the whole method per repetition: CASTOR-KRFE explores its own
		# k_min..k_max range internally, so the repetition is the only loop here
		for repeat in range(args.repeat_start, args.repeat_start + args.repeats):
			# Build the paths of this organism/repetition
			run_directory = os.path.join(outputs, organism, "rep" + str(repeat))
			config_file = os.path.join(configs_tmp, organism, "rep" + str(repeat) + ".ini")
			done_marker = os.path.join(run_directory, ".done")
			# Skip the repetitions already completed (resumable sweep)
			if os.path.exists(done_marker) and not args.force:
				logging.info("[%s rep=%d] already done, skipped", organism, repeat)
				skipped += 1
				continue
			# Clean any partial result of a previous attempt
			if os.path.isdir(run_directory): shutil.rmtree(run_directory)
			# Announce the beginning of the run
			logging.info("[%s rep=%d] starting", organism, repeat)
			# Run the repetition, catching every failure to keep the sweep alive
			try:
				# Write the temporary configuration file
				writeConfiguration(available[organism], config_file, (args.k_min, args.k_max),
					run_directory, overrides.get(organism, dict()))
				# Run the analysis
				returncode, duration, log_path = runCombination(run_directory, config_file, args.timeout)
				# A non zero return code is a failure
				if returncode != 0:
					raise RuntimeError("run_analysis.py exited with code " + str(returncode) +
						" (see " + log_path + ")")
				# Get the metrics of the run
				metrics = parseMetrics(log_path)
				# Collect the classification report, the warning is not fatal for
				# a run that actually completed
				if not saveReport(reports_path, organism, repeat, log_path):
					logging.warning("[%s rep=%d] no classification report found in %s",
						organism, repeat, log_path)
				# Mark the repetition as completed
				with open(done_marker, "w") as f:
					json.dump({"organism": organism, "repeat": repeat, "duration_s": round(duration, 2),
						"finished_at": datetime.datetime.now().isoformat(timespec = "seconds"),
						**metrics}, f, indent = 2)
				# Append the run to the consolidated summary
				appendSummary(summary_path, organism, repeat, "success", duration, metrics, "")
				# Announce the success
				logging.info("[%s rep=%d] success in %.1fs (k = %s, F1 = %s, k-mers = %s)",
					organism, repeat, duration, metrics["k_length"] or "n/a",
					metrics["f1_score"] or "n/a", metrics["n_k_mers"] or "n/a")
				succeeded += 1
			except Exception as exception:
				# Register the error in the error log and continue with the next repetition
				message = type(exception).__name__ + ": " + str(exception)
				error_log.error("%s\trep=%d\t%s", organism, repeat, message.replace("\n", " "))
				logging.error("[%s rep=%d] FAILED - %s", organism, repeat, message)
				appendSummary(summary_path, organism, repeat, "failure", 0, parseMetrics(""), message)
				failed += 1
	# Display the final report of the sweep
	logging.info("Sweep finished: %d succeeded, %d failed, %d skipped (summary: %s)",
		succeeded, failed, skipped, summary_path)
	logging.info("Classification reports collected in: %s", reports_path)
	# A failure of a combination is not a failure of the sweep
	return 0


# Function to append a run to the consolidated summary
def appendSummary(summary_path, organism, repeat, status, duration, metrics, message):
	# Header of the summary file (k_length is the length of k the method selected)
	header = ["timestamp", "organism", "repeat", "status", "duration_s", "f1_score", "k_length",
		"n_k_mers", "message"]
	# Check if the header has to be written
	exists = os.path.exists(summary_path)
	# Append the run
	with open(summary_path, "a", newline = "") as f:
		writer = csv.writer(f)
		if not exists: writer.writerow(header)
		writer.writerow([datetime.datetime.now().isoformat(timespec = "seconds"), organism, repeat, status,
			round(duration, 2), metrics["f1_score"], metrics["k_length"], metrics["n_k_mers"], message])


# Function to parse the command line arguments
def parseArguments(argv):
	parser = argparse.ArgumentParser(description = "Repeat CASTOR-KRFE over several organisms, "
		"each run exploring the k_min..k_max range predefined in the configuration file")
	parser.add_argument("--organisms", nargs = "*", default = None,
		help = "Organisms to process (default: every configuration_{organism}.ini found)")
	parser.add_argument("--k-min", type = int, default = None,
		help = "Override k_min (default: the value predefined in the configuration file)")
	parser.add_argument("--k-max", type = int, default = None,
		help = "Override k_max (default: the value predefined in the configuration file)")
	parser.add_argument("--repeats", type = int, default = 100,
		help = "Number of repetitions of the method per organism (default: 100)")
	parser.add_argument("--repeat-start", type = int, default = 1,
		help = "Index of the first repetition (default: 1)")
	parser.add_argument("--reports", default = None,
		help = "Csv collecting the classification reports (default: {outputs}/reports.csv)")
	parser.add_argument("--config-dir", default = ROOT, help = "Directory of the original .ini files")
	parser.add_argument("--configs-tmp", default = os.path.join(ROOT, "configs_tmp"),
		help = "Directory of the temporary .ini files")
	parser.add_argument("--outputs", default = os.path.join(ROOT, "outputs"),
		help = "Directory of the outputs")
	parser.add_argument("--logs", default = os.path.join(ROOT, "logs"), help = "Directory of the logs")
	parser.add_argument("--timeout", type = float, default = None,
		help = "Maximum duration of a single organism/k run, in seconds")
	parser.add_argument("--force", action = "store_true",
		help = "Rerun the organism/k combinations already completed")
	parser.add_argument("--override", action = "append", default = [], metavar = "ORG:FIELD=VALUE",
		help = "Override a field of the temporary configuration of an organism (repeatable)")
	return parser.parse_args(argv)


# Entry point
if __name__ == "__main__":
	sys.exit(sweep(parseArguments(sys.argv[1:])))
