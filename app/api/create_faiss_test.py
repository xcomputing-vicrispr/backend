import numpy as np
from pydantic import BaseModel
from app.api.lookUpsgRNA import IUPAC_MAP
from Bio.Seq import Seq
from fastapi import APIRouter
import gffutils, re, os, faiss, pickle, smtplib
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
import pandas as pd
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
import gffutils, re, os, faiss, pickle, smtplib
from .nonModel import GenomeUpdate, update_genome_status



router = APIRouter()

def updatePath():
    global PARENT_DIR, DATA_DIR
    PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    DATA_DIR = os.path.join(PARENT_DIR, "data")
    return

def get_paths(user_id: int, genome_name: str):
    base_name = f"nmd_{user_id}_{genome_name}"
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(parent_dir, "data")
    
    paths = {
        "parent_dir": parent_dir,
        "data_dir": data_dir,
        "fasta_path": os.path.join(data_dir, f"{base_name}.fa"),
        "pkl_path": os.path.join(data_dir, f"{base_name}.pkl"),
        "ori_pkl_path": os.path.join(data_dir, f"{base_name}ori.pkl"),
    }
    return paths


def find_sgRNAs_with_PAM_v2(seq: str, chrom: str, pam: str, seed_len=9):
    res = []
    pam_pattern_xuoi = pam_to_regex(pam)
    pattern_xuoi = re.compile(f"(?={pam_pattern_xuoi})")
    seq_upper = seq.upper()
    for m in pattern_xuoi.finditer(seq_upper):
        pam_start = m.start() + 1
        kx = seq_upper[m.start() - 20 : m.start() + len(pam)]
        if kx.startswith("CC") and kx.endswith("GG"):
            continue
        res.append({
            "chrom": f'{chrom}:{pam_start-20}-{pam_start-1}',
            "strand": "+",
            "seq": seq_upper[m.start() - 20 : m.start() + len(pam)],
            "seq_no_pam": seq_upper[m.start() - 20 : m.start()],
            "seed": seq_upper[m.start() - seed_len : m.start()]
        })
    print(len(res))
    seq_upper_nguoc = str(Seq(seq_upper).reverse_complement())
    
    for m in pattern_xuoi.finditer(seq_upper_nguoc):
        pam_start = m.start() + 1
        res.append({
            "chrom": f'{chrom}:{pam_start-20}-{pam_start-1}',
            "strand": "-",
            "seq": seq_upper_nguoc [m.start() - 20 : m.start() + len(pam)],
            "seq_no_pam": seq_upper_nguoc [m.start() - 20 : m.start()],
            "seed": seq_upper_nguoc [m.start() - seed_len : m.start()]
        })
    print(len(res))
    return res
def find_sgRNAs_with_PAM(seq: str, chrom: str, pam: str, seed_len=9):
    res = []
    pam_pattern_xuoi = pam_to_regex(pam)
    pattern_xuoi = re.compile(f"(?={pam_pattern_xuoi})")
    seq_upper = seq.upper()
    for m in pattern_xuoi.finditer(seq_upper):
        pam_start = m.start() + 1
        kx = seq_upper[m.start() - 20 : m.start() + len(pam)]
        if kx.startswith("CC") and kx.endswith("GG"):
           continue
        res.append({
            "chrom": f'{chrom}:{pam_start-20}-{pam_start-1}',
            "strand": "+",
            "seq": seq_upper[m.start() - 20 : m.start() + len(pam)],
            "seq_no_pam": seq_upper[m.start() - 20 : m.start()],
            "seed": seq_upper[m.start() - seed_len : m.start()]
        })
    print(len(res))
    pam_pattern_nguoc = pam_to_regex(str(Seq(pam).reverse_complement()))
    pattern_nguoc = re.compile(f"(?={pam_pattern_nguoc})")
    
    for m in pattern_nguoc.finditer(seq_upper):
        pam_start = m.start() + 1
        seed = seq_upper[m.start() + len(pam): m.start() + len(pam) + seed_len]
        seed = str(Seq(seed).reverse_complement())

        seq = str(Seq(seq_upper[m.start() : m.start() + len(pam) + 20]).reverse_complement())
        seq_no_pam = str(Seq(seq_upper[m.start() + len(pam): m.start() + len(pam) + 20]).reverse_complement())
        

        res.append({
            "chrom": f'{chrom}:{pam_start + len(pam)}-{pam_start + len(pam) + 20 - 1}',
            "strand": "-",
            "seq": seq,
            "seq_no_pam": seq_no_pam,
            "seed": seed
        })
    print(len(res))
    return res

NUC2ONEHOT = {
    'A': 0b1000,
    'C': 0b0100,
    'G': 0b0010,
    'T': 0b0001
}

def seq_to_bits(seq):
    """
    Mã hóa sequence thành one-hot binary vector
    Mỗi nucleotide = 4 bits
    """
    bits = 0
    for base in seq:
        bits = (bits << 4) | NUC2ONEHOT[base]
    
    # Convert to uint8 array
    bit_length = len(seq) * 4
    byte_length = (bit_length + 7) // 8
    return np.frombuffer(bits.to_bytes(byte_length, byteorder='big'), dtype=np.uint8)

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

mapping = {
    'A': np.array([1,0,0,0], dtype=np.float32),
    'C': np.array([0,1,0,0], dtype=np.float32),
    'G': np.array([0,0,1,0], dtype=np.float32),
    'T': np.array([0,0,0,1], dtype=np.float32)
}

def hamming(s1: str, s2: str):
    x = 0
    return sum(char1 != char2 for char1, char2 in zip(s1, s2))

def pam_to_regex(pam):
    return "".join(IUPAC_MAP.get(base, base) for base in pam.upper())

def one_hot_encode(seq):
    return np.concatenate([mapping[base] for base in seq.upper()])

def find_sgRNAs(seq: str, chrom: str, pam="NGG", guide_len=20):
    res = []
    pam_pattern_xuoi = pam_to_regex(pam)
    pattern_xuoi = re.compile(f"(?={pam_pattern_xuoi})")
    seq_upper = seq.upper()
    for m in pattern_xuoi.finditer(seq_upper):
        pam_start = m.start() + 1
        res.append({
            "chrom": f'{chrom}:{pam_start-20}-{pam_start-1}',
            "strand": "+",
            "seq_no_pam": seq_upper[m.start() - 20 : m.start()],
        })

    print(len(res))

    pam_pattern_nguoc = pam_to_regex(str(Seq(pam).reverse_complement()))
    pattern_nguoc = re.compile(f"(?={pam_pattern_nguoc})")
    
    for m in pattern_nguoc.finditer(seq_upper):
        pam_start = m.start() + 1
        seq_no_pam = str(Seq(seq_upper[m.start() + len(pam): m.start() + len(pam) + guide_len]).reverse_complement())
        res.append({
            "chrom": f'{chrom}:{pam_start + len(pam)}-{pam_start + len(pam) + guide_len - 1}',
            "strand": "-",
            "seq_no_pam": seq_no_pam,
        })
    print(len(res))
    return res

def load_all_metadata_from_pkl(file_path):
    all_metadata = []
    with open(file_path, 'rb') as f:
        while True:
            try:
                # Đọc từng đối tượng pickle (từng lô)
                batch_metadata = pickle.load(f)
                # Nối lô vào danh sách tổng
                all_metadata.extend(batch_metadata)
            except EOFError:
                break
    return all_metadata

class BuildFaissIndexRequest(BaseModel):
    owner_id: int
    genome_name: str
    PAM: str = "NGG"
    sgRNA_length: int = 20

class QueryFaissIndexRequest(BaseModel):
    owner_id: int
    genome_name: str
    PAM: str = "NGG"
    sgRNA_length: int = 20
    seed_length: int = 9
    hamming_distance: int = 3
    flank_before: int = 100
    flank_after: int = 100
    maillist: list[str] = []

@router.post("/runFaissPipeline")
async def runFaissPipeline(data: QueryFaissIndexRequest):

    #neu nhu availible thi ok ko thi cut

    update_data = GenomeUpdate(gname=data.genome_name,owner_id=data.owner_id, gw_state="navailible")
    update_genome_status(update_data)

    from .tasks import run_pipeline
    run_pipeline.delay(
        data.owner_id,
        data.genome_name,
        data.PAM,
        data.sgRNA_length,
        data.seed_length,
        data.hamming_distance,
        data.flank_before,
        data.flank_after,
        data.maillist
    )
    return {"status": "pending in queue"}