import os
import subprocess
import sys

def get_file_path(filename):
    return os.path.join("app", "data", filename)

def get_gff3_canonical_hash(filepath):
    """
    1. bỏ header/comment.
    2.  bỏ dòng trống hoặc dòng chỉ có dấu cách/tab.
    3. sort theo mã ASCII.
    4. md5sum: hashing nhi phan
    """
    try:
        cmd = (
            f"export LC_ALL=C; "
            f"grep -v '^#' {filepath} | "
            f"grep -v '^[[:space:]]*$' | "
            f"sed 's/[[:space:]]*$//' | "
            f"sort -S 25% --parallel=4 | "
            f"md5sum"
        )
        
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT).decode()
        return output.split()[0]
    except subprocess.CalledProcessError as e:
        print(f"Error {filepath}: {e.output.decode()}")
        return None

def check_gff3_binary(file_a, file_b):
    path1 = get_file_path(file_a)
    path2 = get_file_path(file_b)

    if not os.path.exists(path1) or not os.path.exists(path2):
        return "NO"

    hash_a = get_gff3_canonical_hash(path1)
    hash_b = get_gff3_canonical_hash(path2)

    if hash_a and hash_b and hash_a == hash_b:
        return "YES"
    else:
        return "NO"

if __name__ == "__main__":
    f1 = "atest.gff3"
    f2 = "btest.gff3"
    
    if len(sys.argv) >= 3:
        f1 = sys.argv[1]
        f2 = sys.argv[2]
        
    print(check_gff3_binary(f1, f2))