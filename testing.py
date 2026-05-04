import os
import subprocess

def get_file_path(filename):
    return os.path.join("app", "data", filename)

def get_multi_fasta_hash(filepath):
    try:
        cmd = f"seqkit sort -n {filepath} 2> /dev/null | seqkit seq -w 0 -u | md5sum"
        
        output = subprocess.check_output(cmd, shell=True).decode()
        return output.split()[0]
    except Exception as e:
        print(f"Error: {e}")
        return None

def main():
    file_a = "atest.fa"
    file_b = "btest.fa"
    path1 = get_file_path(file_a)
    path2 = get_file_path(file_b)

    if not os.path.exists(path1) or not os.path.exists(path2):
        print("File not found.")
        return

    print(f"comparing: {file_a} VS {file_b}")
    
    hash_a = get_multi_fasta_hash(path1)
    hash_b = get_multi_fasta_hash(path2)

    if hash_a and hash_b:
        print(f" Hash A (Sorted): {hash_a}")
        print(f" Hash B (Sorted): {hash_b}")

        if hash_a == hash_b:
            print("\SAME NHAU")
        else:
            print("\n DIFF")

if __name__ == "__main__":
    main()