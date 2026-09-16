# Imports
import os
import sys

# Use a non interactive backend (the cluster nodes have no display)
os.environ.setdefault("MPLBACKEND", "Agg")

import ml
import krfe
import configuration


# Function to run the complete analysis (extraction, training, prediction)
def run_analysis(parameters):
	print("\nCASTOR-KRFE: extraction mode\n", flush = True)
	krfe.extract(parameters)
	print("\nCASTOR-KRFE: training mode\n", flush = True)
	ml.fit(parameters)
	print("\nCASTOR-KRFE: testing mode\n", flush = True)
	ml.predict(parameters)


# Entry point: python3 run_analysis.py <configuration_file.ini>
if __name__ == "__main__":
	# The configuration file is mandatory
	if len(sys.argv) != 2:
		print("Usage: python3 run_analysis.py <configuration_file.ini>", file = sys.stderr)
		sys.exit(2)
	# Get the parameters of the given configuration file
	parameters = configuration.getParameters(configuration_file = sys.argv[1])
	# Run the analysis
	run_analysis(parameters)
