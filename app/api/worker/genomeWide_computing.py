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
from app.configs import get_settings

settings = get_settings()



router = APIRouter()

def updatePath():
    global PARENT_DIR, DATA_DIR
    PARENT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    DATA_DIR = os.path.join(PARENT_DIR, "data")
    return

def get_paths(user_id: int, genome_name: str):
    base_name = f"nmd_{user_id}_{genome_name}"
    parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    data_dir = os.path.join(parent_dir, "data")
    
    paths = {
        "parent_dir": parent_dir,
        "data_dir": data_dir,
        "fasta_path": os.path.join(data_dir, f"{base_name}.fa"),
        "anno_path": os.path.join(data_dir, f"{base_name}.gff3"),
        "filtered_anno_path": os.path.join(data_dir, f"nmd_{user_id}_{genome_name}_only_genes.gff3"),
        "pkl_path": os.path.join(data_dir, f"{base_name}.pkl"),
        "ori_pkl_path": os.path.join(data_dir, f"{base_name}ori.pkl"),
        "name_file": os.path.join(data_dir, f"gw_{base_name}.csv"),
    }
    return paths

def find_sgRNAs_with_PAM(seq: str, chrom: str, pam: str, seed_len=9, seq_len=20):
    res = []
    pam_pattern_xuoi = pam_to_regex(pam)
    pattern_xuoi = re.compile(f"(?={pam_pattern_xuoi})")
    seq_upper = seq.upper()
    for m in pattern_xuoi.finditer(seq_upper):
        pam_start = m.start() + 1
        kx = seq_upper[m.start() - seq_len : m.start() + len(pam)]
        if kx.startswith("CC") and kx.endswith("GG"):
           continue
        res.append({
            "chrom": f'{chrom}:{pam_start-seq_len}-{pam_start-1}',
            "strand": "+",
            "seq": seq_upper[m.start() - seq_len : m.start() + len(pam)],
            "seq_no_pam": seq_upper[m.start() - seq_len : m.start()],
            "seed": seq_upper[m.start() - seed_len : m.start()]
        })
    print(len(res))
    pam_pattern_nguoc = pam_to_regex(str(Seq(pam).reverse_complement()))
    pattern_nguoc = re.compile(f"(?={pam_pattern_nguoc})")
    
    for m in pattern_nguoc.finditer(seq_upper):
        pam_start = m.start() + 1
        seed = seq_upper[m.start() + len(pam): m.start() + len(pam) + seed_len]
        seed = str(Seq(seed).reverse_complement())

        seq = str(Seq(seq_upper[m.start() : m.start() + len(pam) + seq_len]).reverse_complement())
        seq_no_pam = str(Seq(seq_upper[m.start() + len(pam): m.start() + len(pam) + seq_len]).reverse_complement())
        

        res.append({
            "chrom": f'{chrom}:{pam_start + len(pam)}-{pam_start + len(pam) + seq_len - 1}',
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
            "chrom": f'{chrom}:{pam_start-guide_len}-{pam_start-1}',
            "strand": "+",
            "seq_no_pam": seq_upper[m.start() - guide_len : m.start()],
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

def load_filtered_genes(gff_path):
    genes_list = []
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            genes_list.append({
                'chrom': parts[0],
                'start': int(parts[3]),
                'end': int(parts[4]),
                'strand': parts[6],
                'attributes': parts[8]
            })
    return genes_list

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


def send_and_cleanup_data(file_path, maillist, settings):

    MAX_SIZE_MB = 15
    max_size_bytes = MAX_SIZE_MB * 1024 * 1024
    
    files_to_attach = []
    
    try:
        file_size = os.path.getsize(file_path)
        if file_size > max_size_bytes:
            part_num = 1
            with open(file_path, 'rb') as f:
                while True:
                    chunk = f.read(max_size_bytes)
                    if not chunk:
                        break
                    part_name = f"{file_path}.part{part_num}"
                    with open(part_name, 'wb') as p:
                        p.write(chunk)
                    files_to_attach.append(part_name)
                    part_num += 1
            print(f"Split into {len(files_to_attach)} phần.")
        else:
            files_to_attach.append(file_path)

        message = MIMEMultipart()
        message["From"] = settings.SMTP_SENDER
        message["To"] = ', '.join(maillist)
        message["Subject"] = f"Data from ViCRISPR - {os.path.basename(file_path.split('_')[-1])}"

        for p_path in files_to_attach:
            with open(p_path, "rb") as f:
                part = MIMEApplication(f.read(), Name=os.path.basename(file_path.split('_')[-1]))
                part['Content-Disposition'] = f'attachment; filename="{os.path.basename(file_path.split('_')[-1])}"'
                message.attach(part)

        with smtplib.SMTP(settings.SMTP_HOST, settings.SMTP_PORT) as server:
            server.starttls()
            server.login(settings.SMTP_SENDER, settings.SMTP_PASSWORD)
            server.sendmail(settings.SMTP_SENDER, maillist, message.as_string())
            print("Email sent successfully")

    except Exception as e:
        print(f"Error: {e}")
    
    finally:
        if os.path.exists(file_path):
            os.remove(file_path)
        
        for p_path in files_to_attach:
            if os.path.exists(p_path) and p_path != file_path:
                os.remove(p_path)

def buildFaissIndex(owner_id: int, genome_name: str, PAM: str, sgRNA_length: int):

    updatePath()
    paths = get_paths(owner_id, genome_name)
    fasta_path = paths["fasta_path"]
    pkl_path = paths["pkl_path"]
    ori_pkl_path = paths["ori_pkl_path"]
    anno_path = paths["anno_path"]
    filtered_anno_path = paths["filtered_anno_path"]

    with open(anno_path, 'r') as f_in, open(filtered_anno_path, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) > 2 and parts[2] == 'gene':
                f_out.write(line)

    name = f"nmd_{owner_id}_{genome_name}"

    faiss_path = os.path.join(DATA_DIR, f"{name}.faiss")

    sgRNA_len = sgRNA_length
    dim = sgRNA_len * 4

    # tạo index cơ bản + bọc IndexIDMap
    dim_bits = sgRNA_length * 4  # 2 bits/nu (A,C,G,T)
    index_flat = faiss.IndexBinaryFlat(dim_bits)
    index = faiss.IndexBinaryIDMap(index_flat)

    all_sgRNAs_metadata = []
    current_id = 0

    # add theo lô
    print("bat dau xay faiss index")
    with open(ori_pkl_path, 'ab') as f_meta:
        for record in SeqIO.parse(fasta_path, "fasta"):
            chrom = record.id
            seq = str(record.seq)

            sgRNAs_for_chrom = find_sgRNAs(seq, chrom, PAM, sgRNA_length)
            valid_sgRNAs = [g for g in sgRNAs_for_chrom if not any(nu not in 'ACGT' for nu in g["seq_no_pam"]) and len(g["seq_no_pam"]) == sgRNA_length]
            
            if not valid_sgRNAs:
                print(f"Không tìm thấy sgRNA hợp lệ nào trên nhiễm sắc thể {chrom}")
                continue

            # Lưu metadata của lô hiện tại
            pickle.dump(valid_sgRNAs, f_meta) 
            
            # Mã hóa lô sgRNA thành vector NumPy
            vectors_for_chrom = np.vstack([seq_to_bits(g["seq_no_pam"]) for g in valid_sgRNAs]).astype(np.uint8)

            
            # Tạo ID duy nhất cho lô này
            ids_for_chrom = np.arange(current_id, current_id + len(valid_sgRNAs)).astype(np.int64)
            
            # Thêm lô vector và ID vào Faiss index
            index.add_with_ids(vectors_for_chrom, ids_for_chrom)
            
            # Cập nhật ID cho lô tiếp theo
            current_id += len(valid_sgRNAs)

            print(f"xong {chrom}. {len(valid_sgRNAs)}")
        
    # --- 3. Lưu Index và Metadata ---
    print(f"{index.ntotal}")
    faiss.write_index_binary(index, faiss_path)

    print(f"faiss index {faiss_path}")

    index = None

    #bắt đầu gộp metadata từ file cũ
    all_sgRNAs_merged = load_all_metadata_from_pkl(ori_pkl_path)
    print(f"gộp {len(all_sgRNAs_merged)} bản")

    #ghi sang pkl mới
    with open(pkl_path, 'wb') as f:
        pickle.dump(all_sgRNAs_merged, f)
    print(f"file new {pkl_path}")

    os.remove(ori_pkl_path)
    print(f"xoa file {ori_pkl_path}")

    return

def queryFaissIndex(owner_id: int, genome_name: str, PAM: str, sgrna_length: int,
                          seed_length: int, hamming_distance: int, flank_up: int, flank_down: int, maillist: list[str]):

    print(maillist)

    updatePath()
    paths = get_paths(owner_id, genome_name)
    fasta_path = paths["fasta_path"]
    pkl_path = paths["pkl_path"]
    output_file_path = paths["name_file"]
    filtered_anno_path = paths["filtered_anno_path"]

    name = f"nmd_{owner_id}_{genome_name}"
    faiss_path = os.path.join(DATA_DIR, f"{name}.faiss")
    index_loaded = faiss.read_index_binary(faiss_path)


    with open(pkl_path, 'rb') as f:
        sgRNAs_loaded = pickle.load(f)

    PAM = PAM
    seed_len = seed_length
    all_sgRNAs = []
    for record in SeqIO.parse(fasta_path, "fasta"):
            chrom = record.id
            seq = str(record.seq)

            sgRNAs_for_chrom = find_sgRNAs_with_PAM(seq, chrom, PAM, seed_len)
            all_sgRNAs.extend(sgRNAs_for_chrom)
    print(len(all_sgRNAs))
   # sgRNAs_for_chrom = find_sgRNAs_with_PAM_v2(seq, chrom, PAM, seed_len)
    seed_strand_to_indices = defaultdict(list)
    # Gom nhóm theo (seed, strand)
    for i, g in enumerate(all_sgRNAs):
        key = (g["seed"])
        seed_strand_to_indices[key].append(i)

    filtered = []
    # Giữ lại những sgRNA có seed + strand xuất hiện duy nhất
    for seed, idxs in seed_strand_to_indices.items():
        if len(idxs) == 1:
            only_index = idxs[0]
            filtered.append(all_sgRNAs[only_index])

    # Tạo danh sách kết quả
    seq_list = [
        {
            "seq_no_pam": g["seq_no_pam"],
            "chrom": g["chrom"],
            "strand": g["strand"],
            "seq": g["seq"],
            "kc": ""
        }
        for g in filtered
    ]
    print(len(seq_list))
    count_cc_gg = sum(
    1 for g in seq_list 
    if g["seq"].startswith("CC") and g["seq"].endswith("GG")
    )
    # df = pd.DataFrame(seq_list)
    # output_file = os.path.join(DATA_DIR, 'gw.csv')
    # df.to_csv(output_file, index=False)
    print(len(seq_list))
    i = 0
    encoded_queries = []
    total = len(seq_list)
    for query in seq_list:
        sequence = query["seq_no_pam"]
        if any(nuc not in "ACGT" for nuc in sequence):
            continue
        if len(sequence) != sgrna_length:
            continue
        vector = seq_to_bits(sequence)
        encoded_queries.append(vector)
        i = i + 1
        if i % 100 == 0:
            print(f"Progress embedding: {i}/{total} ({i/total:.2%})", end='\r')
    query_vectors = np.vstack(encoded_queries).astype(np.uint8)

    print("encoded xong")
    sgRNAs_final = []
    dl = []
    vt = set()
    D, I = index_loaded.search(query_vectors, k=3)
    hammingDistance = hamming_distance

    print("tim xong")
    newdt = []
    i = 0
    for q, idx_row, dist_row in zip(seq_list, I, D):
        idx = idx_row[1]

        s1 = sgRNAs_loaded[idx_row[0]]  # xâu gần nhất
        s2 = sgRNAs_loaded[idx_row[1]]  # xâu thứ 2
        s3 = sgRNAs_loaded[idx_row[2]]  # xâu thứ 3
        smss1 = s1["seq_no_pam"]
        smss2 = s2["seq_no_pam"]
        smss3 = s3["seq_no_pam"]
        dist0 = hamming(q["seq_no_pam"], smss1)
        dist1 = hamming(q["seq_no_pam"], smss2)
        dist2 = hamming(q["seq_no_pam"], smss3)

        if (dist1 < hammingDistance):
            continue
        sgRNAs_final.append(
            {
                "seq_no_pam": q["seq_no_pam"],
                "chrom": q["chrom"],
                "strand": q["strand"],
                "seq": q["seq"],
                "kc": str(dist0) + "/" + str(dist1) + "/" + str(dist2),
                "details": str(s1) + "," + str(s2) + ", " + str(s3)
            }
        )
        i = i + 1
        if i % 100 == 0:
            print(f"Progress loc theo distance: {i}/{total} ({i/total:.2%})", end='\r')
        #print(q, sgRNAs_loaded[idx], dist1)
        #print(sgRNAs_loaded[idx])

    # df = pd.DataFrame(newdt)
    # output_file = os.path.join(DATA_DIR, 'gw.csv')
    # df.to_csv(output_file, index=False)
    # return
    print("So luong sgRNA sau loc:", len(sgRNAs_final))
    output_file = os.path.join(DATA_DIR, 'testing.csv')
    print("Đã uả")

    genes_data = load_filtered_genes(filtered_anno_path)
    results = []
    for j, sgRNA in enumerate(sgRNAs_final, 1):

        print(j)
        chrom_info = sgRNA["chrom"].split(':')
        chrom = chrom_info[0]
        position = chrom_info[1]
        sg_start = int(position.split('-')[0])
        sg_end = int(position.split('-')[1])
        
        for gene in genes_data:
            if gene['chrom'] != chrom:
                continue
                
            if gene['strand'] == '+':
                tss_start = gene['start'] - flank_up
                tss_end = gene['start'] + flank_down
            else:
                tss_start = gene['end'] - flank_down
                tss_end = gene['end'] + flank_up

            if not (sg_end < tss_start or sg_start > tss_end):
                print(f"Match: {sgRNA['seq']} | Pos: {gene['start']}-{gene['end']} | {j}/{len(sgRNAs_final)}")
                
                results.append({
                    'sgRNA_seq': sgRNA["seq"],
                    'sgRNA_loc': sgRNA["chrom"],
                    'Strand': sgRNA["strand"],
                    'gene_start': gene['start'],
                    'gene_end': gene['end'],
                    'gene_strand': gene['strand'], # Lưu thêm strand của gene để đối chiếu
                    'kc': sgRNA["kc"],
                    'details': sgRNA["details"],
                    'gene_data': gene['attributes']
                })                
                break
    df = pd.DataFrame(results)
    df.to_csv(output_file_path, index=False)
    print("Đã lưu kết quả")
    print(count_cc_gg)

    send_and_cleanup_data(output_file_path, maillist, settings)
    return

def cleanFaissIndex(owner_id: int, genome_name: str):
    updatePath()
    paths = get_paths(owner_id, genome_name)
    pkl_path = paths["pkl_path"]
    output_file_path = paths["name_file"]

    name = f"nmd_{owner_id}_{genome_name}"
    faiss_path = os.path.join(DATA_DIR, f"{name}.faiss")

    if os.path.exists(faiss_path):
        os.remove(faiss_path)
        print(f"Deleted Faiss index file: {faiss_path}")

    if os.path.exists(pkl_path):
        os.remove(pkl_path)
        print(f"Deleted pickle file: {pkl_path}")

    if os.path.exists(output_file_path):
        os.remove(output_file_path)
        print(f"Deleted output CSV file: {output_file_path}")