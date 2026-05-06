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
from sqlalchemy import func

from .worker.computing import GeneNameComputing, CoordinateComputing, FastaComputing
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
redis_client_fq_celery = redis.Redis.from_url(
    settings.REDIS_FQ_URL,
    decode_responses=True
)

def get_and_decr_redis(redis_client, key: str):    
    value_before_str = redis_client.get(key) 
    
    value_before = int(value_before_str) if value_before_str else 0
    print(f"Giá trị TRƯỚC khi giảm của {key}: {value_before}")
    
    value_after = redis_client.decr(key)
    
    print(f"Giá trị SAU khi giảm của {key}: {value_after}")
    
    return value_after 


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
            
            save_sgRNA_list_dbv(
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
        final_value = get_and_decr_redis(redis_client_fq_celery, redis_key)
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
def uploadNonModel_celery(self, redis_key, session_id, user_id, fa_name, anno_name, display_id):
    flag_key = f"task_decr_flag:{user_id}:{self.request.id}"
    processes = []
    cleanup_done = False
    task_cancelled = False
    
    def cleanup_all_processes(signum=None, frame=None):
        nonlocal cleanup_done, task_cancelled
        if cleanup_done:
            return
        cleanup_done = True
        task_cancelled = True

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
    
    def check_cancelled():
        if task_cancelled:
            raise Terminated("Task was cancelled by user")
    
    original_sigterm = signal.signal(signal.SIGTERM, cleanup_all_processes)
    original_sigint = signal.signal(signal.SIGINT, cleanup_all_processes)
    
    try:
        start = time.time()

        name_fa = fa_name.split(".")[0]
        name_anno = anno_name.split(".")[0]

        # Temporary working filenames (legacy pattern, used during processing only)
        nonmdFA = f"nmd_{user_id}_{name_fa}.fa"
        nonmdAN = f"nmd_{user_id}_{name_anno}.gff3"

        print(f"[DEDUP] Starting upload: fa={name_fa}, anno={name_anno}")
        update_data = GenomeUpdate(display_id=display_id, status="prepare processing data", log="")
        update_genome_status(update_data)
        
        check_cancelled()
        
        # ── Step 1: Locate temp files ──
        temp_fasta = None
        found = False
        for ext in [".fa", ".fna"]:
            temp_fasta = os.path.join(DATA_DIR, f"{session_id}_temp_nmd_{user_id}_{name_fa}{ext}")
            if os.path.exists(temp_fasta):
                found = True
                break
        if not found:
            raise FileNotFoundError("Không tìm thấy file temp fasta")
        
        temp_anno = None
        found = False
        for ext in [".gtf", ".gff", ".gff3"]:
            temp_anno = os.path.join(DATA_DIR, f"{session_id}_temp_nmd_{user_id}_{name_anno}{ext}")
            if os.path.exists(temp_anno):
                found = True
                break
        if not found:
            raise FileNotFoundError("Không tìm thấy file temp annotation")

        check_cancelled()

        # ── Step 2: Rename temp fasta to working path ──
        fasta_path = os.path.join(DATA_DIR, nonmdFA)
        anno_path = os.path.join(DATA_DIR, nonmdAN)

        os.rename(temp_fasta, fasta_path)

        check_cancelled()

        # ── Step 3: AGAT convert annotation ──
        cmd = [
            "agat_convert_sp_gxf2gxf.pl",
            "--gff", temp_anno,
            "-o", anno_path,
            "-v", "2"
        ]

        print("[DEDUP] Converting annotation file with AGAT...")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, preexec_fn=os.setpgrp)
        processes.append(proc)
        stdout, stderr = proc.communicate()
        
        check_cancelled()
        
        if proc.returncode != 0:
            raise Exception(f"AGAT conversion failed with code {proc.returncode}: {stderr}")
        
        if not os.path.exists(anno_path):
            raise FileNotFoundError(f"AGAT output file not created: {anno_path}")
        
        if os.path.getsize(anno_path) == 0:
            raise Exception(f"AGAT created empty file: {anno_path}")
        
        print(f"[DEDUP] AGAT conversion successful. Output size: {os.path.getsize(anno_path)} bytes")

        # Clean temp annotation source (temp_fasta already renamed)
        if os.path.exists(temp_anno):
            os.remove(temp_anno)
        
        check_cancelled()

        # ── Step 4: Two-phase FASTA deduplication ──
        from .hashing import find_existing_fasta, find_existing_gff3

        update_data = GenomeUpdate(display_id=display_id, status="hashing fasta", log="")
        update_genome_status(update_data)

        fasta_size = os.path.getsize(fasta_path)

        db_session = SessionLocal()
        try:
            existing_fasta_id, canonical_fasta_hash = find_existing_fasta(db_session, fasta_path)
        finally:
            db_session.close()

        fasta_is_duplicate = existing_fasta_id is not None
        fasta_canonical = canonical_fasta_hash or existing_fasta_id

        if not fasta_canonical:
            raise Exception("Failed to compute fasta hash")

        check_cancelled()

        if fasta_is_duplicate:
            # ── FASTA DUPLICATE: skip 2bit + bowtie build ──
            print(f"[DEDUP] Fasta content already exists (hash={fasta_canonical}). Skipping 2bit/bowtie.")
            # Remove the working fasta file since shared copy already exists
            if os.path.exists(fasta_path):
                os.remove(fasta_path)
        else:
            # ── FASTA NEW: create shared dir and process ──
            print(f"[DEDUP] New fasta content (hash={fasta_canonical}). Processing...")
            shared_fasta_dir = os.path.join(DATA_DIR, f"fasta_{fasta_canonical}")
            os.makedirs(shared_fasta_dir, exist_ok=True)

            shared_fasta_file = os.path.join(shared_fasta_dir, f"{fasta_canonical}.fa")
            shutil.move(fasta_path, shared_fasta_file)

            # faToTwoBit
            update_data = GenomeUpdate(display_id=display_id, status="converting fa to 2bit", log="")
            update_genome_status(update_data)
            cmd_2bit = f"{os.path.join(DATA_DIR, 'faToTwoBit')} {fasta_canonical}.fa {fasta_canonical}.2bit"
            proc = subprocess.Popen(cmd_2bit, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=shared_fasta_dir, preexec_fn=os.setpgrp)
            processes.append(proc)
            stdout, stderr = proc.communicate()
            proc.wait()
            check_cancelled()
            if proc.returncode != 0:
                raise Exception(f"2bit conversion failed: {stderr}")
            print(f"[DEDUP] 2bit created in {time.time() - start:.1f}s")

            # bowtie-build
            update_data = GenomeUpdate(display_id=display_id, status="bowtie indexing", log="")
            update_genome_status(update_data)
            cmd_bowtie = f"bowtie-build {fasta_canonical}.fa {fasta_canonical}_index"
            proc = subprocess.Popen(cmd_bowtie, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=shared_fasta_dir, preexec_fn=os.setpgrp)
            processes.append(proc)
            stdout, stderr = proc.communicate()
            proc.wait()
            check_cancelled()
            if proc.returncode != 0:
                raise Exception(f"Bowtie build failed: {stderr}")
            print(f"[DEDUP] Bowtie index created in {time.time() - start:.1f}s")

        # ── Step 5: Two-phase GFF3 deduplication ──
        update_data = GenomeUpdate(display_id=display_id, status="hashing annotation", log="")
        update_genome_status(update_data)

        anno_size = os.path.getsize(anno_path)

        db_session = SessionLocal()
        try:
            existing_gff3_id, canonical_gff3_hash = find_existing_gff3(db_session, anno_path)
        finally:
            db_session.close()

        gff3_is_duplicate = existing_gff3_id is not None
        gff3_canonical = canonical_gff3_hash or existing_gff3_id

        if not gff3_canonical:
            raise Exception("Failed to compute gff3 hash")

        check_cancelled()

        if gff3_is_duplicate:
            # ── GFF3 DUPLICATE: skip gffutils + sorted files ──
            print(f"[DEDUP] GFF3 content already exists (hash={gff3_canonical}). Skipping gffutils/sort.")
            if os.path.exists(anno_path):
                os.remove(anno_path)
        else:
            # ── GFF3 NEW: create shared dir and process ──
            print(f"[DEDUP] New GFF3 content (hash={gff3_canonical}). Processing...")
            shared_anno_dir = os.path.join(DATA_DIR, f"anno_{gff3_canonical}")
            os.makedirs(shared_anno_dir, exist_ok=True)

            shared_gff3_file = os.path.join(shared_anno_dir, f"{gff3_canonical}.gff3")
            shutil.move(anno_path, shared_gff3_file)

            # gffutils
            update_data = GenomeUpdate(display_id=display_id, status="gffutils indexing", log="")
            update_genome_status(update_data)
            check_cancelled()

            gffdb_path = os.path.join(shared_anno_dir, f"{gff3_canonical}.db")

            print(f"[DEDUP] Creating gffutils database from {shared_gff3_file}")
            gffutils.create_db(
                shared_gff3_file,
                dbfn=gffdb_path,
                merge_strategy="create_unique",
                keep_order=True,
                disable_infer_transcripts=True,
                disable_infer_genes=True,
                force=True
            )
            print(f"[DEDUP] gffutils db created in {time.time() - start:.1f}s")

            # createMMRegion — adapted for shared dir
            _createMMRegion_shared(shared_gff3_file, shared_anno_dir, gff3_canonical)

        # ── Step 6: Update DB with hashes ──
        update_data = GenomeUpdate(
            display_id=display_id,
            status="success",
            log="genome ready"
        )
        update_genome_status(update_data)

        # Update hash columns directly
        db_session = SessionLocal()
        try:
            genome = db_session.query(Genome).filter(Genome.id_for_user_display == display_id).first()
            if genome:
                genome.id_use_for_us_fasta = fasta_canonical
                genome.id_use_for_us_gff3 = gff3_canonical
                genome.fasta_size = fasta_size
                genome.anno_size = anno_size
                genome.upload_timestamp = func.now()
                db_session.commit()
                print(f"[DEDUP] Genome record updated: fasta_hash={fasta_canonical}, gff3_hash={gff3_canonical}")
                if fasta_is_duplicate:
                    print(f"[DEDUP] Fasta files SHARED with existing genome")
                if gff3_is_duplicate:
                    print(f"[DEDUP] GFF3 files SHARED with existing genome")
        finally:
            db_session.close()

        # Clean up any remaining working files
        for f in [fasta_path, anno_path]:
            if os.path.exists(f):
                os.remove(f)

        print(f"[DEDUP] Total time: {time.time() - start:.1f}s")
    
    except (SoftTimeLimitExceeded, Terminated, KeyboardInterrupt) as e:
        print(f"Task terminated: {type(e).__name__}")
        cleanup_all_processes()
        update_data = GenomeUpdate(display_id=display_id, status="Cancelled", log="Task was terminated")
        update_genome_status(update_data)
        return {"status": "cancelled", "message": "Task was cancelled by user"}
    except Exception as e:
        print(f"Task failed: {e}")
        cleanup_all_processes()
        update_data = GenomeUpdate(display_id=display_id, status="Failed", log=str(e))
        update_genome_status(update_data)
        raise
    finally:
        signal.signal(signal.SIGTERM, original_sigterm)
        signal.signal(signal.SIGINT, original_sigint)
        cleanup_all_processes()
        
        try:
            cleanup_temp_files_by_session(user_id, session_id)
        except Exception as cleanup_err:
            print(f"[CLEANUP] Failed to clean temp files: {cleanup_err}")
        
        cuop_duoc_co = redis_client_fq_celery.set(
            flag_key,
            "task",
            ex=7200,
            nx=True
        )
        if cuop_duoc_co:
            final_value = get_and_decr_redis(redis_client_fq_celery, redis_key)
            print(f"Task hoàn thành. Key đã được giảm về {final_value}")


def _createMMRegion_shared(gff3_file: str, output_dir: str, gff3_canonical: str):
    """
    Create sorted exon/gene region files in the shared annotation directory.
    Adapted from createMMRegion() for shared storage.
    """
    exon_out = os.path.join(output_dir, f"{gff3_canonical}_exons.sorted.gff3")
    gene_out = os.path.join(output_dir, f"{gff3_canonical}_genes.sorted.gff3")

    cmd_exon = f"awk 'BEGIN{{OFS=\"\\t\"}} $3 == \"exon\"' {gff3_file} | sort -k1,1V -k4,4n > {exon_out}"
    cmd_gene = f"awk 'BEGIN{{OFS=\"\\t\"}} $3 == \"gene\"' {gff3_file} | sort -k1,1V -k4,4n > {gene_out}"

    try:
        subprocess.run(cmd_exon, shell=True, check=True, executable='/bin/bash')
        subprocess.run(cmd_gene, shell=True, check=True, executable='/bin/bash')
        print(f"[DEDUP] Created sorted files:\n  {exon_out}\n  {gene_out}")
    except subprocess.CalledProcessError as e:
        print(f"[DEDUP] Error creating sorted region files: {e}")



@celery.task(bind=True, queue='genomeWide')
def run_pipeline(self, display_id, pam, sgrna_length,
                 seed_region, hamming_distance, flank_up, flank_down, emails):
    try:
        # build faiss index
        buildFaissIndex(display_id, pam, sgrna_length)
        # goi query
        queryFaissIndex(display_id, pam, sgrna_length,
                          seed_region, hamming_distance, flank_up, flank_down, emails)
    except Exception as e:
        print(f"Error in genome wide pipeline: {e}")
    finally:
        update_data = GenomeUpdate(display_id=display_id, gw_state="available")
        update_genome_status(update_data)
        cleanFaissIndex(display_id)
        print("tra availible genome wide cho genome")




    


