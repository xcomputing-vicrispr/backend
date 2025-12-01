import subprocess
import os

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAR_PARENT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORK_DIR = os.path.join(PAR_PARENT_DIR, "Lindel")
input_file = os.path.join(WORK_DIR, "lindel_test_seq.txt")

DATA_DIR = os.path.join(PARENT_DIR, "data")

output_file = os.path.join(DATA_DIR, "lindel_output.txt")

def calLindelScore(data):

    with open(input_file, "w") as f:
        for i, seq in enumerate(data, start=1):
            f.write(f"{seq}")

    command = [
        "python3",
        "Lindel_prediction.py",
        "-f",
        "lindel_test_seq.txt",
        "-o",
        "/dev/stdout",
    ]

    process = subprocess.Popen(
        command,
        cwd=WORK_DIR,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    stdout_data, stderr_data = process.communicate()

    if process.returncode != 0:
        print("Return code:", process.returncode)
        print("Lindel error336:", stderr_data)

    return stdout_data

print(calLindelScore(["CACCGGCATCTTTGTAGCTAAGAGAGGTTTTATCGGTCACTGCTTGGGTCCCCACGCGTT"]))
