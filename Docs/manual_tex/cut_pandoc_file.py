"""
This python script modifies a file by removing all 
lines after a line that contains the string "___PANDOC_IGNORE___"

Usage: 
    python cut_pandoc_file.py filename.tex

to modify the file filename.tex. The file is modified in place.
"""

import sys
from pathlib import Path

filename = sys.argv[1]

assert Path(filename).exists() and filename.endswith(".tex")

lines = []
with open(filename, "r") as f:
    line = f.readline()
    while line and "___PANDOC_IGNORE___" not in line:
        lines.append(line)
        line = f.readline()
    
with open(filename, "w") as f:
    f.writelines(lines)
