#!/usr/bin/env python3
import os
import subprocess
import sys
from pathlib import Path

def run_cmd(cmd):
    """Running shell and print log."""
    print(f"\n>>> Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    print(result.stderr)
        # Chỉ in ra STDERR nếu lệnh trả về mã lỗi (khác 0)
    if result.returncode != 0:
        print("--- AGAT ERROR LOG (STDERR) ---")
        print(result.stderr.strip())
    return result.returncode


afilename = "nmd_2_ecoli26.gff"
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PARENT_DIR, "data")
input_file = os.path.join(DATA_DIR, afilename)
output_file = afilename.split(".")[0] + ".gff3"
OUTPUT_DIR = os.path.join(DATA_DIR, output_file)

convert_cmd = ["conda", "run", "-n", "agat", "agat_convert_sp_gxf2gxf.pl", "-g", str(input_file), "-o", str(OUTPUT_DIR)]
ret = run_cmd(convert_cmd)

if ret == 0:
    print(f"Đã chuyển đổi xong -> {OUTPUT_DIR}")
else:
    print("Lỗi khi chuyển đổi.")
