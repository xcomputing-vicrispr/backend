from fastapi import APIRouter, HTTPException
from pydantic import BaseModel



import subprocess
import smtplib, numpy as np
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from itertools import product
from app.database import get_db, SessionLocal
from app.models import EmailQueue
from app.configs import get_settings

settings = get_settings()

import re, os, json, random, string
from Bio.Seq import Seq

router = APIRouter()

sender = "vicrispr@gmail.com"
password = "yghnybppmeaehpxs"

IUPAC_MAP = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "W": ["A", "T"],
    "R": ["A", "G"],
    "M": ["A", "C"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "K": ["G", "T"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"]
}

class MailSession(BaseModel):
    idfile: str
    mail_list: list


def sendMail(idfile, maillist):
    db = SessionLocal()
    header_data = db.query(TaskMetadata).filter(TaskMetadata.query_id == idfile).first()
    sgrna_list = db.query(Sgrna).filter(Sgrna.query_id == idfile).order_by(Sgrna.stt.asc()).all()
    sgrna_list_dict = []

    metadata_dict = {
        'name': header_data.query_name,
        'spec': header_data.spec,
        'pam': header_data.pam,
        'sgRNA_len': header_data.sgrna_len,
        'gene_strand': header_data.gene_strand,
        'type_task': header_data.type_task,
        'min_product_size': header_data.min_product_size,
        'max_product_size': header_data.max_product_size,
        'min_primer_size': header_data.min_primer_size,
        'max_primer_size': header_data.max_primer_size,
        'optimal_primer_size': header_data.optimal_primer_size,
        'min_tm': header_data.min_tm,
        'max_tm': header_data.max_tm,
        'optimal_tm': header_data.optimal_tm,
        'status': header_data.status,
        'log': header_data.log,
        'queue_task_id': header_data.queue_task_id
    }

    for sgrna in sgrna_list:
        sgrna_dict = {
            'sequence': sgrna.sequence,
            'location': sgrna.location,
            'strand': sgrna.strand,
            'GC Content': sgrna.gc_content,
            'Self-complementary': sgrna.self_complementary,
            'mm0': sgrna.mm0,
            'mm1': sgrna.mm1,
            'mm2': sgrna.mm2,
            'mm3': sgrna.mm3,
            'cfdScore': sgrna.cfd_score,
            'mlScore': sgrna.ml_score,
            'microScore': sgrna.micro_score,
            'mmejpre': sgrna.mmej_pre,
            'Secondary structure with scaffold': sgrna.sec_structure,
            'bowtie_details': sgrna.bowtie_details,
            'mismatch_region': sgrna.mismatch_region,
            'rs3': sgrna.rs3_score,
        }
        sgrna_list_dict.append(sgrna_dict)
    
    json_data = [metadata_dict] + sgrna_list_dict
    json_str = json.dumps(json_data, indent=2, ensure_ascii=False)

    genename = metadata_dict['name']
    pam = metadata_dict['pam']
    spec = metadata_dict['spec']

    if sgrna_list_dict:
        headers = list(sgrna_list_dict[0].keys())
        tsv_header = "\t".join(headers)
        tsv_body = "\n".join(
            "\t".join(str(row.get(h, "")) for h in headers)
            for row in sgrna_list_dict
        )
        tsv_str = f"# genename: {genename}\n# spec: {spec}\n{tsv_header}\n{tsv_body}"
    else:
        tsv_str = f"# genename: {genename}\n# spec: {spec}\n"


    message = MIMEMultipart()
    message["From"] = settings.SMTP_SENDER
    message["To"] = ', '.join(maillist)
    message["Subject"] = "Data từ ViCRISPR"

    json_part = MIMEApplication(json_str.encode('utf-8'), Name="result.json")
    json_part['Content-Disposition'] = 'attachment; filename="result.json"'
    message.attach(json_part)

    tsv_part = MIMEApplication(tsv_str.encode('utf-8'), Name="result.tsv")
    tsv_part['Content-Disposition'] = 'attachment; filename="result.tsv"'
    message.attach(tsv_part)

    with smtplib.SMTP(settings.SMTP_HOST, settings.SMTP_PORT) as server:
        server.starttls()
        server.login(settings.SMTP_SENDER, settings.SMTP_PASSWORD)
        server.sendmail(settings.SMTP_SENDER, maillist, message.as_string())
        print("Email sent successfully")


def write_sgrna_to_fasta2(sgrnas, filename="app/data/sgrna_output.fa"):
    with open(filename, "w") as f:
        index = 0
        for seq in sgrnas:
            prefix = seq[:-3]
            prefix=prefix.upper()
            f.write(f">{index}\n{prefix}TGG\n")
            index += 1
            f.write(f">{index}\n{prefix}AGG\n")
            index += 1
            f.write(f">{index}\n{prefix}GGG\n")
            index += 1
            f.write(f">{index}\n{prefix}CGG\n")
            index += 1

def write_sgrna_to_fasta_with_NNAGAAW(sgrnas, filename="app/data/sgrna_output.fa"):
    # Định nghĩa các nucleotide có thể cho N và W
    N_bases = ['A', 'T', 'C', 'G']
    W_bases = ['A', 'T']

    with open(filename, "w") as f:
        index = 0
        for seq in sgrnas:
            prefix = seq[:-7]  # loại bỏ 7 ký tự PAM ở cuối
            # Tạo tất cả biến thể PAM
            for n1 in N_bases:
                for n2 in N_bases:
                    for w in W_bases:
                        pam = f"{n1}{n2}AGAA{w}"
                        f.write(f">{index}\n{prefix}{pam}\n")
                        index += 1

def write_sgrna_to_fasta_with_IUPAC(sgrnas, pam, idfile):
    index = 0 
    bases_per_position = [IUPAC_MAP[base] for base in pam]
    all_variants = list(product(*bases_per_position))
    base_dir="./app/data"
    output_filename = f"{base_dir}/{idfile}_sgrna_output.fa"
    with open(output_filename, "w") as f:
        for seq in sgrnas:
            prefix = seq[:-len(pam)]
            for all_posi in all_variants:
                pam_sequence = "".join(all_posi)
                full_sequence = prefix + pam_sequence
                f.write(f">{index}\n{full_sequence}\n")
                index += 1  

def leve(s1: str, s2: str):
    l = len(s1)
    dp = []
    for i in range(l + 5):
        row = [0] * (l + 5)
        dp.append(row)

    #n voi m kieu j cung bang nhau
    for i in range(l + 1):
        dp[i][0] = i
        dp[0][i] = i

    for i in range(1, l + 1):
        for j in range(1, l + 1):
            if s1[i - 1] == s2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                dp[i][j] = min(dp[i - 1][j - 1], dp[i - 1][j], dp[i][j - 1]) + 1

    return dp[l][l]
    
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PARENT_DIR, "data")
OUTPUT_DIR = "app/data"


def count_permu_IUPAC(s: str) -> int:
    iupac_counts = {
        "A": 1, "C": 1, "G": 1, "T": 1,
        "R": 2, "Y": 2, "S": 2, "W": 2,
        "K": 2, "M": 2,
        "B": 3, "D": 3, "H": 3, "V": 3,
        "N": 4
    }
    res = 1
    for ch in s.upper():
        res *= iupac_counts.get(ch, 1)
    return res
    
def xuly(line: str, datafile, pam_name: str, ll: int, off_target: int = 0, num_of_mismatches: int = 0):

    res = []
    parts = line.strip().split('\t')

    permu_num = count_permu_IUPAC(pam_name)

    x = int(int(parts[0]) / permu_num)
    if len(parts) < 8:
        datafile[x]["mm0"] += 1
        return -1, -1

    mismatch_info = parts[7]
    mismatch_info_parse = mismatch_info.split(",")
    positions = []
    
    for item in mismatch_info_parse[0:]:
        pos_str = item.split(":")[0]  
        positions.append(int(pos_str))
    if any(pos >= ll for pos in positions):
        return -1, -1
    
    #Chi phuc vu cho critirea thu 2
    if off_target == 1:
        mismatch_count_in_seed = sum(1 for pos in positions if pos > 9)#distal region la tu 0-8, con lai la seed (ko tinh PAM)
        if mismatch_count_in_seed > num_of_mismatches:
            return -1, -1

    if not mismatch_info:
        mis_match_num = 0
    else:
        mis_match_num = mismatch_info.count(":")

    res = f"{parts[2]}:{parts[3]}", parts[4], mismatch_info
    print(f"{x},{parts[1]},{parts[2]},{parts[3]},{parts[4]},{mismatch_info}")

    leveshtein_dis = 0
    mutated_seq = ""
    # if parts[1] == "+":
    #     seq1 = parts[4][0:ll]
    #     mutations = mismatch_info.split(",")

    #     seq_list = list(seq1)   
    #     for db in mutations:
    #         pos, change = db.split(":")
    #         ref, alt = change.split(">")
    #         pos = int(pos)
    #         seq_list[pos] = ref  
    #     mutated_seq = "".join(seq_list)
    #     leveshtein_dis = leve(mutated_seq, seq1)
    # else:
    #     seq1 = str(Seq(parts[4]).reverse_complement())
    #     seq1 = seq1[0:ll]
    #     mutations = mismatch_info.split(",")

    #     seq_list = list(seq1)
    #     for db in mutations:
    #         pos, change = db.split(":")
    #         ref, alt = change.split(">")
    #         pos = int(pos)
    #         seq_list[pos] = ref  
    #     mutated_seq = "".join(seq_list)
    #     mutated_seq = str(Seq(mutated_seq).reverse_complement())
    #     leveshtein_dis = leve(mutated_seq, seq1)
    res = f"{parts[2]}:{parts[3]}", mutated_seq, mismatch_info
    if mis_match_num == 1:
        datafile[x]["mm1"] += 1
    elif mis_match_num == 2:
        datafile[x]["mm2"] += 1
    elif mis_match_num == 3:
        datafile[x]["mm3"] += 1

    return x, res

import numpy as np

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def gc_score(gc_percent, len_without_pam):

    score = 1 - ((gc_percent - 50) / len_without_pam) ** 2
    return max(0, round(score, 4))

@router.post("/submitSendMail")
async def submitSendMail(data: MailSession):

    # from .tasks import submitSendMail_celery
    # task = submitSendMail_celery.delay(data.dict())
    # return {"status": "submitted", "task_id": task.id}

    #check xem cai idfile nay da xong chua, neu status la ok roi thi gui di
    #Neu chua xong, neu nhu trong db chua ton tai idfile nay: tao mot bang ghi trong db bao gom: idfile + maillist
    #Neu chua xong, nhung db da ton tai idfile nay: tao list gom mailist, roi gop no voi mailist moi sau do unique
    #               roi lai luu quay tro lai db
    try:
        db = SessionLocal()

        header_data = db.query(TaskMetadata).filter(TaskMetadata.query_id == data.idfile).first()

        status = header_data.status
        
        if status == 'success':
            sendMail(data.idfile, data.mail_list)
            return {
                "status": "sent",
                "message": "Check your mailbox in some minutes"
            }
        else:
            existing_emails = {
                email[0] for email in 
                db.query(EmailQueue.email).filter(EmailQueue.idfile == data.idfile).all()
            }
            print(existing_emails)
            new_queue_items = [
                EmailQueue(idfile=data.idfile, email=email)
                for email in data.mail_list
                if email not in existing_emails
            ]
            print(2232323232323)

            if new_queue_items:
                db.bulk_save_objects(new_queue_items)
                db.commit()
            

            return {
                "status": "queued",
                "message": f"Đã thêm {len(new_queue_items)} email mới vào hàng đợi.",
            }
            #ghi vao trong database idfile va mailist de cho luc nao xong thi moi gui
    except Exception as e:
        return {"error": str(e)}




print(leve("GAAAACAAGAGACGAAGCGG", "GAAAATAAGAGACGAAGCGA"))  # 2
