from celery import Celery
import time
from pydantic import BaseModel
import httpx
import subprocess
import time, gffutils
import asyncio

import smtplib, numpy as np
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
import signal
from celery.exceptions import SoftTimeLimitExceeded, Terminated

from .cmdprocess import get_fasta_from_twobit
from Bio.Seq import Seq
import re, os, json, random, string
from .lookUpsgRNA import *
from .nonModel import *
from .cfdEffiencyScore import get_cfd_score
from .mlEffiencyScore import get_ml_score, get_ml_score_azi3
from .calLindel import calLindelScore
from .export import sendMail, xuly, MailSession, sigmoid, gc_score, write_sgrna_to_fasta2, write_sgrna_to_fasta_with_NNAGAAW
from .export import sender, password, OUTPUT_DIR, DATA_DIR, write_sgrna_to_fasta_with_IUPAC

from .worker.computing import GeneNameComputing, CoordinateComputing, FastaComputing, indexComputing, getMMRegion
from .worker.genomeWide_computing import updatePath, buildFaissIndex, queryFaissIndex, cleanFaissIndex

import redis
import os
from celery.signals import task_prerun
from app.configs import get_settings

settings = get_settings()

QUEUE_POSITION_PREFIX = "vicrispr_queue:"

celery = Celery(
    "crispr_app",
    broker=settings.RABBITMQ_BROKER_URL,
    backend=settings.REDIS_BACKEND_URL
)
#DB1
redis_client = redis.Redis.from_url(
    url=settings.REDIS_APP_URL,  
    decode_responses=True
)

#DB2
redis_client_fq = aioredis.from_url(
    url=settings.REDIS_FQ_URL,
    decode_responses=True
)

@celery.task(queue='send_mail')
def submitSendMail_celery(data):


    data = MailSession(**data)
    seq_list = []
    seq_list_ml = []
    seq_list_lindel = []


    filename = "vcp" + data.idfile + ".json"
    file_path = os.path.join(OUTPUT_DIR, filename)
    with open(file_path, 'r', encoding='utf-8') as f:
        datafile = json.load(f)


    pam_name = datafile[0]["pam"]
    bowtie_index_file = datafile[0]["name"] + "_index"
    datafile = datafile[1:]
    datafile = datafile[:-1]

    for i, row in enumerate(datafile):
        seq_list.append(row["sequence"])
        seq_list_ml.append(row["mlseq"])

    write_sgrna_to_fasta_with_IUPAC(seq_list, pam_name)


    if (pam_name == "NGG"):
        ml_score = get_ml_score(seq_list_ml)
    else:
        ml_score = ["N/S"] * len(datafile)

    command = [
        "bowtie",
        "-v", "3",
        "-a",
        "-f",
        "-x", bowtie_index_file,
        "sgrna_output.fa",
        "/dev/stdout"
    ]

    process = subprocess.Popen(
        command,
        cwd=DATA_DIR,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )



    for line in process.stdout:
        line = line.strip()
        idseq, ttmm = xuly(line, datafile)
        datafile[idseq]["bowtie_details"] += ttmm + "; "
        x, scr = get_cfd_score(line)
        if x == -1:
            continue
        datafile[x]["cfdScore"] += scr
        print(str(x) + " " + str(datafile[x]["cfdScore"]))


    process.stdout.close()
    process.wait()

    for i, row in enumerate(datafile):
        print(
            f"{i}: {row.get('sequence', '')}, "
            f"location={row.get('location', '')}, "
            f"mm0={row.get('mm0', 0)}, "
            f"mm1={row.get('mm1', 0)}, "
            f"mm2={row.get('mm2', 0)}, "
            f"mm3={row.get('mm3', 0)}"
        )

    for i in range(len(datafile)):
        datafile[i]["mlScore"] = ml_score[i]

    sendMail(data.mail_list, datafile, datafile[0]['name'])

    return


@celery.task(bind=True, queue='lookUpComputing')
def lookUpComputing_celery(self, redis_key, idd, request, generalSetting, casData, primerConfigData, type):
    
    try:
        current_task_id = self.request.id
        if type == "gene_name":
            GeneNameComputing(current_task_id, idd, request, generalSetting, casData, primerConfigData)
        elif type == "coordinate":
            CoordinateComputing(current_task_id, idd, request, generalSetting, casData, primerConfigData)   
        else:
            FastaComputing(current_task_id, idd, request, generalSetting, casData, primerConfigData)
    except Exception as e:
        print(f"Task {self.request.id} failed với lỗi: {str(e)}")
        
        try:
            gene_name = request.get('gene_name', 'unknown')
            spec = request.get('species', 'unknown')
            PAM = casData.get('pam', 'unknown')
            sgRNA_len = generalSetting.get('sgRNA_len', 20)
            
            save_sgRNA_list(
                idd, [], gene_name, spec, PAM, sgRNA_len, type,
                0, 0, 0, 0, 0, 0, 0, 0,
                queue_task_id=current_task_id,
                status="failed",
                log=f"Error: {str(e)}"
            )
        except Exception as save_error:
            print(f"Không thể save error status: {str(save_error)}")

        raise
    finally:
        final_value = asyncio.run(get_and_decr_redis(redis_client_fq, redis_key))
        print(f"Task hoàn thành. Key đã được giảm về {final_value}")
        
    return

@task_prerun.connect
def on_task_prerun(sender, task_id, **kwargs):
    queue_name = sender.queue
    if queue_name:
        zset_key = f"{QUEUE_POSITION_PREFIX}{queue_name}"
        redis_client.zrem(zset_key, task_id)

import signal
import os
from celery.exceptions import SoftTimeLimitExceeded, Terminated

@celery.task(bind=True, queue='nonModel')
def uploadNonModel_celery(self, redis_key, session_id, user_id, fa_name, anno_name):
    flag_key = f"task_decr_flag:{user_id}:{self.request.id}"
    processes = []
    cleanup_done = False
    task_cancelled = False
    
    def cleanup_all_processes(signum=None, frame=None):
        nonlocal cleanup_done, task_cancelled
        if cleanup_done:
            return
        cleanup_done = True
        task_cancelled = True  # ← Đánh dấu task bị cancel
        
        print(f"Task cancelled - cleaning up {len(processes)} processes...")
        for proc in processes:
            try:
                if proc.poll() is None:
                    print(f"Terminating subprocess PID {proc.pid}")
                    try:
                        pgid = os.getpgid(proc.pid)
                        os.killpg(pgid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    
                    try:
                        proc.wait(timeout=3)
                    except subprocess.TimeoutExpired:
                        try:
                            os.killpg(pgid, signal.SIGKILL)
                            proc.wait(timeout=1)
                        except:
                            pass
            except Exception as e:
                print(f"Error killing subprocess: {e}")
    
    def check_cancelled(): #check xem co bi cencel khong
        if task_cancelled:
            raise Terminated("Task was cancelled by user")
    
    original_sigterm = signal.signal(signal.SIGTERM, cleanup_all_processes)
    original_sigint = signal.signal(signal.SIGINT, cleanup_all_processes)
    
    try:
        start = time.time()
        print("sap ghi xong fasta")

        name_fa = fa_name.split(".")[0]
        name_anno = anno_name.split(".")[0]

        nonmdFA = f"nmd_{user_id}_{name_fa}.fa"
        nonmdAN = f"nmd_{user_id}_{name_anno}.gff3"

        print(name_fa, name_anno, nonmdFA, nonmdAN)
        update_data = GenomeUpdate(gname=name_fa,owner_id=user_id,status="prepare processing data", log="")
        update_genome_status(update_data)
        
        check_cancelled()
        
        temp_fasta = None
        temp_anno = None
        found = False
        for ext in [".fa", ".fna"]:
            temp_fasta = os.path.join(DATA_DIR, f"{session_id}_temp_nmd_{user_id}_{name_fa}{ext}")
            if os.path.exists(temp_fasta):
                found = True
                break
        if not found:
            raise FileNotFoundError("Không tìm thấy file temp")
        
        found = False
        for ext in [".gtf", ".gff", ".gff3"]:
            temp_anno  = os.path.join(DATA_DIR, f"{session_id}_temp_nmd_{user_id}_{name_anno}{ext}")
            if os.path.exists(temp_anno):
                found = True
                break
        if not found:
                raise FileNotFoundError("Ko thay file temp")

        check_cancelled()

        fasta_path = os.path.join(DATA_DIR, nonmdFA)
        anno_path = os.path.join(DATA_DIR, nonmdAN)

        print(fasta_path)
        os.rename(temp_fasta, fasta_path)

        check_cancelled()

        cmd = [
            "conda", "run", "-n", "agat",
            "agat_convert_sp_gxf2gxf.pl",
            "--gff", temp_anno,
            "-o", anno_path,
            "-v", "2"
        ]
        print("chuan bi convert annotaion file")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, preexec_fn=os.setpgrp)
        processes.append(proc)
        stdout, stderr = proc.communicate()
        
        check_cancelled()
        
        print("AGAT stdout:", stdout)
        print("AGAT stderr:", stderr)
        
        if proc.returncode != 0:
            raise Exception(f"AGAT conversion failed with code {proc.returncode}: {stderr}")
        
        if not os.path.exists(anno_path):
            raise FileNotFoundError(f"AGAT output file not created: {anno_path}")
        
        if os.path.getsize(anno_path) == 0:
            raise Exception(f"AGAT created empty file: {anno_path}")
        
        print(f"AGAT conversion successful. Output size: {os.path.getsize(anno_path)} bytes")

        for temp_file in [temp_fasta, temp_anno]:
            if os.path.exists(temp_file):
                os.remove(temp_file)
        
        check_cancelled()
        
        cmd_2bit = "./faToTwoBit" + " " + nonmdFA + " " + nonmdFA.split(".")[0] + ".2bit"
        cmd_bowtie_build = "bowtie-build" + " " + nonmdFA + " " + nonmdFA.split(".")[0] + "_index"

        update_data = GenomeUpdate(gname=name_fa,owner_id=user_id,status="processConverting fa to 2bit", log="")
        update_genome_status(update_data)
        proc = subprocess.Popen(cmd_2bit, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=DATA_DIR, preexec_fn=os.setpgrp)
        processes.append(proc)
        stdout, stderr = proc.communicate()
        proc.wait()
        
        check_cancelled()
        
        if proc.returncode != 0:
            raise Exception(f"2bit conversion failed: {stderr}")

        end_2bit = time.time()
        print("da chuyen fa sang 2bit", end_2bit - start)

        update_data = GenomeUpdate(gname=name_fa,owner_id=user_id,status="bowtie indexing", log="")
        update_genome_status(update_data)
        proc = subprocess.Popen(cmd_bowtie_build, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=DATA_DIR, preexec_fn=os.setpgrp)
        processes.append(proc)
        stdout, stderr = proc.communicate()
        proc.wait()
        
        check_cancelled()
        
        if proc.returncode != 0:
            raise Exception(f"Bowtie build failed: {stderr}")

        end_bowtie = time.time()
        print("da tao bowtie index", end_bowtie - start)

        update_data = GenomeUpdate(gname=name_fa,owner_id=user_id,status="gffutils indexing", log="")
        update_genome_status(update_data)
        
        check_cancelled()  # ← Check TRƯỚC gffutils
        
        db_path = os.path.join(DATA_DIR, f"nmd_{user_id}_{name_anno}.db")
        gff3_file = os.path.join(DATA_DIR, nonmdAN)
        
        print(f"Creating gffutils database from {gff3_file}")
        if not os.path.exists(gff3_file):
            raise FileNotFoundError(f"GFF3 file not found: {gff3_file}")
        
        file_size = os.path.getsize(gff3_file)
        print(f"GFF3 file size: {file_size} bytes")
        if file_size == 0:
            raise Exception(f"GFF3 file is empty: {gff3_file}")
        
        db = gffutils.create_db(
            gff3_file,
            dbfn=db_path,
            merge_strategy="create_unique",
            keep_order=True,
            disable_infer_transcripts=True,
            disable_infer_genes=True,
        )
        print("da tao db")
        end_gff = time.time()
        createMMRegion(nonmdAN)

        update_data = GenomeUpdate(gname=name_fa,owner_id=user_id,status="Success, genome ready", log="")
        update_genome_status(update_data)
        print("da tao gff3", end_gff - start)
    
    except (SoftTimeLimitExceeded, Terminated, KeyboardInterrupt) as e:
        print(f"Task terminated: {type(e).__name__}")
        cleanup_all_processes()
        update_data = GenomeUpdate(gname=fa_name.split(".")[0], owner_id=user_id, status="Cancelled", log="Task was terminated")
        update_genome_status(update_data)
        return {"status": "cancelled", "message": "Task was cancelled by user"}
    except Exception as e:
        print(f"Task failed: {e}")
        cleanup_all_processes()
        update_data = GenomeUpdate(gname=fa_name.split(".")[0], owner_id=user_id, status="Failed", log=str(e))
        update_genome_status(update_data)
        raise
    finally:
        signal.signal(signal.SIGTERM, original_sigterm)
        signal.signal(signal.SIGINT, original_sigint)
        cleanup_all_processes()
        
        cuop_duoc_co = redis_client_fq.set(
            flag_key,
            "task",
            ex=7200,
            nx=True
        )
        if cuop_duoc_co:
            final_value = asyncio.run(get_and_decr_redis(redis_client_fq, redis_key))
            print(f"Task hoàn thành. Key đã được giảm về {final_value}")
        


@celery.task(bind=True, queue='genomeWide')
def run_pipeline(self, uid, genome_name, pam, sgrna_length,
                 seed_region, hamming_distance, flank_up, flank_down, emails):
    try:
        # build faiss index
        buildFaissIndex(uid, genome_name, pam, sgrna_length)
        # goi qury
        queryFaissIndex(uid, genome_name, pam, sgrna_length,
                          seed_region, hamming_distance, flank_up, flank_down, emails)
    except Exception as e:
        print(f"Error in genome wide pipeline: {e}")
    finally:
        update_data = GenomeUpdate(gname=genome_name,owner_id=uid, gw_state="available")
        update_genome_status(update_data)
        cleanFaissIndex(uid, genome_name)
        print("tra availible genome wide cho genome")



    


