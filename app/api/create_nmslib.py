from Bio import SeqIO
import numpy as np
import nmslib, os, pickle
from Bio.Seq import Seq

genome_name = "ecoli"
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PARENT_DIR, "data")
fasta_path = os.path.join(DATA_DIR, f"{genome_name}.fa")
pkl_path = os.path.join(DATA_DIR, f"{genome_name}.pkl")
nms_path = os.path.join(DATA_DIR, f"{genome_name}.nmslib")

def encode_seq(seq):
    mapping = {"A":[1,0,0,0], "C":[0,1,0,0], "G":[0,0,1,0], "T":[0,0,0,1]}
    return np.array([mapping[b] for b in seq], dtype=np.float32).flatten()


def find_sgRNAs(seq, pam="NGG", guide_len=20):
    sgRNAs = []

    # --- Forward strand ---
    for i in range(len(seq) - len(pam)):
        if seq[i] == "N":
            continue
        pam_seq = seq[i:i+len(pam)]
        if pam_seq[1:] == "GG":   # PAM = NGG
            guide_start = i - guide_len
            if guide_start >= 0:
                guide = seq[guide_start:i]
                if(guide.count('N') > 0):
                    continue
                sgRNAs.append((guide, guide_start, "+"))

    # --- Reverse strand ---
    rc_seq = str(Seq(seq).reverse_complement())
    for i in range(len(rc_seq) - len(pam)):
        if rc_seq[i] == "N":
            continue
        pam_seq = rc_seq[i:i+len(pam)]
        if pam_seq[1:] == "GG":   # PAM = NGG
            guide_start = i - guide_len
            if guide_start >= 0:
                guide = rc_seq[guide_start:i]
                if(guide.count('N') > 0):
                    continue
                # Quy đổi ngược lại tọa độ gốc trên forward
                pos = len(seq) - (i + len(pam))
                sgRNAs.append((guide, pos, "-"))
    return sgRNAs


def sgRNA_generator(fasta_path):
    for record in SeqIO.parse(fasta_path, "fasta"):
        chrom = record.id
        seq = str(record.seq)
        guide_pos = set()
        for guide, pos, strand in find_sgRNAs(seq):
            if pos in guide_pos:
                continue
            guide_pos.add(pos)
            print(len(guide_pos), end="\r")
            vec = encode_seq(guide)
            yield vec, f"{chrom}:{pos}:{strand}"
        print(chrom, len(seq), "done")


index = nmslib.init(method='hnsw', space='l2')

batch_size = 100000
batch_vecs = []
batch_ids = []
all_ids = []

for vec, id_str in sgRNA_generator(fasta_path):
    vec /= np.linalg.norm(vec)
    batch_vecs.append(vec)
    batch_ids.append(id_str)
    all_ids.append(id_str)
    if len(batch_vecs) >= batch_size:
        index.addDataPointBatch(batch_vecs)
        with open("ids.txt", "a") as f:
            f.write("\n".join(batch_ids) + "\n")
        batch_vecs.clear()
        batch_ids.clear()

if len(batch_vecs) >= batch_size:
        index.addDataPointBatch(batch_vecs)
        with open("ids.txt", "a") as f:
            f.write("\n".join(batch_ids) + "\n")
        batch_vecs.clear()
        batch_ids.clear()

index.createIndex({'post': 2}, print_progress=True)
index.saveIndex(nms_path, save_data=False)

with open(pkl_path, "wb") as f:
    pickle.dump(all_ids, f)

print("Indexing complete")