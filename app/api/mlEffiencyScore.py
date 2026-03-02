import json
import subprocess
import warnings, os
from rs3.seq import predict_seq
warnings.filterwarnings("ignore")

def get_ml_score(seqlist):
    seq_str = json.dumps(seqlist)
    
    conda_env = os.getenv("ML_CONDA_ENV", "env2_conda")
    
    
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
    TEST_FILE = os.path.abspath(os.path.join(CURRENT_DIR, "worker", "getrs2.py"))
    TEST_DIR = os.path.dirname(TEST_FILE)

    result = subprocess.Popen(
        ["conda", "run", "-n", conda_env, "python", TEST_FILE, seq_str],
        cwd=TEST_DIR,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    # ===========================

    stdout, stderr = result.communicate()
    
    print("=== ML Script STDOUT ===")
    print(seqlist)
    print(stdout.decode())
    print("=== ML Script STDERR ===")
    print(stderr.decode())
    
    try:
       # return json.loads(stdout.decode().strip())
        stdout_lines = stdout.decode().splitlines()
        for line in reversed(stdout_lines):
            line = line.strip()
            if line.startswith("[") and line.endswith("]"):
                return json.loads(line)
    except json.JSONDecodeError as e:
        print(f"JSON decode error: {e}")
        print("Raw output that caused error:")
        print(stdout.decode())
        raise  # Có thể bỏ dòng này nếu bạn không muốn exception tiếp tục

def get_ml_score_azi3(seqlist):
    res = predict_seq(seqlist, sequence_tracr='Hsu2013')
    return res.tolist()
