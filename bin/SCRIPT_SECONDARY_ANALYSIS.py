#!/usr/bin/env python3
# Import dependencies
import os
from datetime import datetime

# Find the latest RESULTS folder
results_folders = [d for d in os.listdir('.') if d.startswith('RESULTS_') and os.path.isdir(d)]
if not results_folders:
    print("ERROR: No RESULTS folder found!")
    exit(1)

latest_results = sorted(results_folders)[-1]
print(f"Working on: {latest_results}")

# Generate Summary File
summary_file = os.path.join(latest_results, 'SECONDARY_ANALYSIS_SUMMARY.txt')
with open(summary_file, 'w') as f:
    f.write("=" * 50 + "\n")
    f.write("SECONDARY ANALYSIS\n")
    f.write("=" * 50 + "\n")
    f.write(f"Execution time: {datetime.now()}\n")
    f.write(f"Results folder: {latest_results}\n")
    f.write("\n")
    f.write("Contents of results folder:\n")
    
    # List all files in the results folder
    for item in sorted(os.listdir(latest_results)):
        item_path = os.path.join(latest_results, item)
        if os.path.isfile(item_path):
            size = os.path.getsize(item_path)
            f.write(f"  [FILE] {item} ({size} bytes)\n")
        elif os.path.isdir(item_path):
            f.write(f"  [DIR]  {item}/\n")

# Function 1

# Finish script
with open(summary_file, 'a') as f:
    f.write("\n" + "=" * 50 + "\n")
    f.write("Secondary Analysis completed successfully!\n")
    f.write("=" * 50 + "\n")
print("Secondary analysis test completed!")
