import nmslib
from Bio import SeqIO
import os
import numpy as np
import time
import pickle
from Bio.Seq import Seq


def encode_seq(seq):
    mapping = {"A":[1,0,0,0], "C":[0,1,0,0], "G":[0,0,1,0], "T":[0,0,0,1]}
    return np.array([mapping[b] for b in seq], dtype=np.float32).flatten()


genome_name = "ecoli"
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PARENT_DIR, "data")
pkl_path = os.path.join(DATA_DIR, f"{genome_name}.pkl")
nms_path = os.path.join(DATA_DIR, f"{genome_name}.nmslib")

index = nmslib.init(method='hnsw', space='l2')
index.loadIndex(nms_path, load_data=False)

with open(pkl_path, "rb") as f:
    all_ids = pickle.load(f)

query_vec = encode_seq("CACGATCGTACGTCGATGCAAACG")
query_vec /= np.linalg.norm(query_vec)
ids, distances = index.knnQuery(query_vec, k=40)

print(ids, distances)

for i, dist in zip(ids, distances):
    print(all_ids[i], dist)
