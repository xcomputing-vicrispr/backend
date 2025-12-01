import gffutils, re, os, faiss, pickle, smtplib
from Bio.Seq import Seq
import numpy as np
from Bio import SeqIO
from collections import defaultdict
import pandas as pd
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication



genome_name = "nmd_2_ccruddi"
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PARENT_DIR, "data")
fasta_path = os.path.join(DATA_DIR, f"{genome_name}.fa")
fasta_path_test = os.path.join(DATA_DIR, f"faku.fa")
pkl_path = os.path.join(DATA_DIR, f"{genome_name}.pkl")

def find_sgRNAs_with_PAM(seq: str, chrom: str, pam: str):
    res = []
    pam_pattern_xuoi = pam_to_regex(pam)
    pattern_xuoi = re.compile(f"(?={pam_pattern_xuoi})")
    seq_upper = seq.upper()
    for m in pattern_xuoi.finditer(seq_upper):
        pam_start = m.start()
        res.append({
            "chrom": f'{chrom}:{pam_start-20}',
            "strand": "+",
            "seq": seq_upper[pam_start - 20 : pam_start + len(pam)],
            "seq_no_pam": seq_upper[pam_start - 20 : pam_start],
            "seed": seq_upper[m.start() - seed_len : m.start()]
        })

    pam_pattern_nguoc = pam_to_regex(str(Seq(pam).reverse_complement()))
    pattern_nguoc = re.compile(f"(?={pam_pattern_nguoc})")
    
    for m in pattern_nguoc.finditer(seq_upper):
        pam_start = m.start()
        seed = seq_upper[m.start() + len(pam): m.start() + len(pam) + seed_len]
        seed = str(Seq(seed).reverse_complement())

        res.append({
            "chrom": f'{chrom}:{pam_start}',
            "strand": "-",
            "seq": seq_upper[pam_start : pam_start + len(pam) + 20],
            "seq_no_pam": seq_upper[pam_start + len(pam): pam_start + len(pam) + 20],
            "seed": seed
        })
    return res

IUPAC_MAP = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "W": "[AT]",
    "R": "[AG]",
    "M": "[AC]",
    "Y": "[CT]",
    "S": "[GC]",
    "K": "[GT]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]"
}

def pam_to_regex(pam):
    return "".join(IUPAC_MAP.get(base, base) for base in pam.upper())

def one_hot_encode(seq):
    return np.concatenate([mapping[base] for base in seq.upper()])

mapping = {
    'A': np.array([1,0,0,0], dtype=np.float32),
    'C': np.array([0,1,0,0], dtype=np.float32),
    'G': np.array([0,0,1,0], dtype=np.float32),
    'T': np.array([0,0,0,1], dtype=np.float32)
}

faiss_index_path = os.path.join(DATA_DIR, f"{genome_name}.faiss")
index_loaded = faiss.read_index(faiss_index_path)

with open(pkl_path, 'rb') as f:
    sgRNAs_loaded = pickle.load(f)

PAM = "NGG"
seed_len = 9
unique_seed = []
all_sgRNAs = []
for record in SeqIO.parse(fasta_path, "fasta"):
        chrom = record.id
        seq = str(record.seq)

        sgRNAs_for_chrom = find_sgRNAs_with_PAM(seq, chrom, PAM)
        all_sgRNAs.extend(sgRNAs_for_chrom)

seed_to_indices = defaultdict(list)
for i, g in enumerate(all_sgRNAs):
    seed_to_indices[g["seed"]].append(i)

filtered = []

for seed, idxs in seed_to_indices.items():
    if len(idxs) == 1:
        only_index = idxs[0]          
        filtered.append(all_sgRNAs[only_index])

seq_list = [
    {
        "seq_no_pam": g["seq_no_pam"],
        "chrom": g["chrom"],
        "strand": g["strand"],
        "seq": g["seq"]    
    }
    for g in filtered
]

print(seq_list)
i = 0
encoded_queries = []
for query in seq_list:
    sequence = query["seq_no_pam"]
    if any(nuc not in "ACGT" for nuc in sequence):
        continue
    if len(sequence) != 20:
        continue
    vector = one_hot_encode(sequence)
    encoded_queries.append(vector)
    i = i + 1
    print(i)
query_vectors = np.array(encoded_queries, dtype=np.float32)
print("encoded xong")
sgRNAs_final = []
dl = []
vt = set()
D, I = index_loaded.search(query_vectors, k=2)
hammingDistance = 3

print("tim xong")
for q, idx_row, dist_row in zip(seq_list, I, D):
    idx = idx_row[1]
    dist = dist_row[1] / 2
    if dist < hammingDistance:
        continue
    dl.append(dist)
    print(q, sgRNAs_loaded[idx], dist)
    sgRNAs_final.append(q)
    print(sgRNAs_loaded[idx])

print(dl)
print("So luong sgRNA sau loc:", len(sgRNAs_final))
output_file = os.path.join(DATA_DIR, 'testing.csv')
print("Đã uả")

db_path = os.path.join(DATA_DIR, f"{genome_name}.db")
db = gffutils.FeatureDB(db_path, keep_order=True)
print(len(sgRNAs_final))
length_sgRNA = 20
flank = 100
results = []
j = 0 
for sgRNA in sgRNAs_final:
    j = j + 1
    loc_str = sgRNA
    position = int(sgRNA["chrom"].split(':')[1])
    sg_start = position
    sg_end = position + length_sgRNA
    
    chrom = sgRNA["chrom"].split(':')[0]
    
    genes = db.region(seqid=chrom, featuretype='gene')
    
    skip = True
    sgRNA_valid = False
    for gene in genes:
        gene_start = gene.start - flank
        gene_end = gene.end + flank
        
        if not (sg_end < gene_start or sg_start > gene_end):
            sgRNA_valid = True
            print(sgRNA, gene.id, gene.start, gene.end, j, len(sgRNAs_final))
            results.append({
                'sgRNA_seq': sgRNA["seq"],
                'sgRNA_loc': sgRNA["chrom"],
                'Strand': sgRNA["strand"],
                'gene_id': gene.id,
                'gene_start': gene.start,
                'gene_end': gene.end,
                'gene_data': '; '.join(f'{k}={",".join(v)}' for k, v in gene.attributes.items())
            })
            break
df = pd.DataFrame(results)
output_file = os.path.join(DATA_DIR, 'gw.csv')
df.to_csv(output_file, index=False)
print("Đã lưu kết quả")

sender = "vicrispr@gmail.com"
password = "yghnybppmeaehpxs"
maillist = ['linhprono01@gmail.com']

message = MIMEMultipart()
message["From"] = sender
message["To"] = ', '.join(maillist)
message["Subject"] = "Data từ ViCRISPR"

file_path = os.path.join(DATA_DIR, "gw.csv")
tsv_part = []
with open(file_path, "rb") as f:
    tsv_part = MIMEApplication(f.read(), Name="gw.csv")
tsv_part['Content-Disposition'] = 'attachment; filename="gw.csv"'
message.attach(tsv_part)

with smtplib.SMTP("smtp.gmail.com", 587) as server:
    server.starttls()
    server.login(sender, password)
    server.sendmail(sender, maillist, message.as_string())
    print("Email sent successfully")