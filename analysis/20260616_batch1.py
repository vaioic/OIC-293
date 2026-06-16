from pathlib import Path

import tissue_analyzer

input_dir = "../data/batch1/"
output_dir = "../processed/2026-06-03 Dev"

tissue_analyzer.process_directory(input_dir, output_dir)
