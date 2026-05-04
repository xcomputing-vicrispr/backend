import subprocess
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq

def get_fasta_from_twobit(twobit_file: str, chromosome: str, chrstart: int, chrstop: int) -> str:
    try: 
        PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        DATA_DIR = os.path.join(PARENT_DIR, "data")

        twoBitToFa_path = os.path.join(DATA_DIR, "twoBitToFa")

        if(int(chrstart) > int(chrstop)):
            tmp = chrstart
            chrstart = chrstop
            chrstop = tmp

        twobit_file = os.path.basename(twobit_file)

        NEW_DATA_DIR = os.path.join(DATA_DIR, f"fasta_{twobit_file.split('.')[0]}")

        print(NEW_DATA_DIR)

        command_list = [
            twoBitToFa_path,
            f"{twobit_file}:{chromosome}:{chrstart}-{chrstop}",
            "/dev/stdout"
        ]

        result = subprocess.run(
            command_list,
            shell=False,
            capture_output=True,
            text=True,
            cwd=NEW_DATA_DIR,
        )
        parts = result.stdout.split('\n', 1)
        return parts[1].replace('\n', '')
    except Exception as e:
        raise Exception(f"twoBitfromFasta running failed: {str(e)}, please check again query and genome") from e


def extract_exon_by_gene(keyw: str):
    # Chuyển đến thư mục chứa dữ liệu
    os.chdir("/mnt/d/UET/CHOCHO/backend/data")

    cmd = f'grep "\\"{keyw}\\";" hg38.gtf > ext1.gtf'

    subprocess.run(cmd, shell=True)

    cmd = f"awk '$3 == \"exon\"' ext1.gtf > ext1_exon.gtf"

    subprocess.run(cmd, shell=True)

    cmd = f"awk '{{print $1\"\\t\"$4\"\\t\"$5\"\\t\"$7}}' ext1_exon.gtf > vitri.bed"

    subprocess.run(cmd, shell=True)

    cmd = f"./twoBitToFa osativa70.2bit -bed=vitri.bed output.fa"

    subprocess.run(cmd, shell=True)
    results = []

    for record in SeqIO.parse("output.fa", "fasta"):
        seq = str(record.seq)
        l = len(seq)
        for i in range(l - 2):
            if (seq[i] == 'C' and seq[i + 1] == 'C') or  (seq[i] == 'C' and seq[i + 1] == 'c') or (seq[i] == 'c' and seq[i + 1] == 'C') or (seq[i] == 'c' and seq[i + 1] == 'c'):
                pam_seq = seq[i:min(i + 23, l)]
                if len(pam_seq) != 23: continue
                pam_seq = Seq(pam_seq)
                rev_comp = str(pam_seq.reverse_complement())
                results.append(("-", rev_comp, i + 1))

        for i in range(l - 2):
            if (seq[i+2] == 'G' and seq[i + 1] == 'G') or (seq[i+2] == 'g' and seq[i + 1] == 'G') or (seq[i+2] == 'g' and seq[i + 1] == 'g') or (seq[i+2] == 'G' and seq[i + 1] == 'g'):
                pam_seq = seq[max(0, i - 20):i + 3]
                if len(pam_seq) != 23: continue
                results.append(("+", pam_seq, i + 1))

    return results

def fold_rna(seq: str):
    result = subprocess.run(
        ["RNAfold"],
        input=seq.encode(),
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL
    )
    output = result.stdout.decode().splitlines()
    if len(output) < 2:
        return None, None
    structure = output[1].split(" ")[0]
    mfe = float(output[1].split("(")[-1].replace(")", ""))
    return structure, mfe




