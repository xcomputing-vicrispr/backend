from fastapi import APIRouter, HTTPException, File, UploadFile, Depends
from sqlalchemy.orm import Session
from sqlalchemy import func
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from .cfdEffiencyScore import get_cfd_score
from .mlEffiencyScore import get_ml_score
from typing import Literal, Optional
from fastapi import Form, Body, Request
from app.database import get_db, SessionLocal
from app.models import Genome
from celery.result import AsyncResult
import subprocess, aiofiles, shutil, glob
import smtplib, numpy as np
import redis.asyncio as aioredis
from app.configs import get_settings
import redis, asyncio
import uuid

settings = get_settings() 


import re, os, json, random, string
import gffutils, time

redis_client_fq = aioredis.from_url(
    url=settings.REDIS_FQ_URL,
    decode_responses=True
)

redis_client_fq_celery = redis.Redis.from_url(
    settings.REDIS_FQ_URL,
    decode_responses=True
)
MAX_CONCURRENT_TASKS = 4

router = APIRouter()
class GenomeCreate(BaseModel):
    gname: str
    owner_id: int
    status: str
    log: str
    kbstorage: int
    session_id: str


PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PARENT_DIR, "data")
OUTPUT_DIR = "app/data"
TMP_DIR = os.path.join(DATA_DIR, "tmp")   # nơi giữ các chunk tạm
CHUNK_PREFIX = "chunk_"                    # tiền tố đặt tên chunk

def resolve_genome_paths(genome: Genome) -> dict:
    """
    Resolve actual file paths for a genome record.
    """
    result = {}
    if genome.id_use_for_us_fasta:
        h = genome.id_use_for_us_fasta
        fasta_dir = os.path.join(DATA_DIR, f"fasta_{h}")
        result["fasta"] = os.path.join(fasta_dir, f"{h}.fa")
        result["twobit"] = os.path.join(fasta_dir, f"{h}.2bit")
        result["bowtie_index"] = os.path.join(fasta_dir, f"{h}_index")
        result["fasta_dir"] = fasta_dir
    if genome.id_use_for_us_gff3:
        h = genome.id_use_for_us_gff3
        anno_dir = os.path.join(DATA_DIR, f"anno_{h}")
        result["gff3"] = os.path.join(anno_dir, f"{h}.gff3")
        result["gffdb"] = os.path.join(anno_dir, f"{h}.db")
        result["exon_sorted"] = os.path.join(anno_dir, f"{h}_exons.sorted.gff3")
        result["gene_sorted"] = os.path.join(anno_dir, f"{h}_genes.sorted.gff3")
        result["anno_dir"] = anno_dir
    return result


def generate_user_display_id() -> str:
    """Generate a UUID v4 for user display."""
    return str(uuid.uuid4())
class FileCheck(BaseModel):
    user_id: int
    name: str

class GenomeSignUp(BaseModel):
    gname: str
    file_size: int
    owner_id: int

class GenomeUpdate(BaseModel):
    display_id: Optional[str] = None
    gname: Optional[str] = None              
    owner_id: Optional[int] = None           
    status: Optional[str] = None
    log: Optional[str] = None
    task_queue_id: Optional[str] = None
    gw_state: Optional[str] = None
    
def get_and_decr_redis(redis_client, key: str):    
    value_before_str = redis_client.get(key) 
    
    value_before = int(value_before_str) if value_before_str else 0
    print(f"Giá trị TRƯỚC khi giảm của {key}: {value_before}")
    
    value_after = redis_client.decr(key)
    
    print(f"Giá trị SAU khi giảm của {key}: {value_after}")
    
    return value_after 


def update_genome_status(update: GenomeUpdate):
    db = SessionLocal()
    try:
        if update.display_id:
            genome = db.query(Genome).filter(Genome.id_for_user_display == update.display_id).first()
        elif update.owner_id is not None and update.gname:
            genome = db.query(Genome).filter(
                (Genome.owner_id == update.owner_id) & (Genome.gname == update.gname)
            ).first()
        else:
            raise HTTPException(status_code=400, detail="Missing identifiers for Genome update")
        if not genome:
            raise HTTPException(status_code=404, detail="Genome không tồn tại")

        # Cập nhật nếu có giá trị mới
        if update.status is not None:
            genome.status = update.status
        if update.log is not None:
            genome.log = update.log
        if update.task_queue_id is not None:
            genome.task_queue_id = update.task_queue_id
        if update.gw_state is not None:
            genome.gw_state = update.gw_state
        
        if update.status and update.status.lower() in ["success", "failed", "cancelled", "error", "no_result"]:
            genome.upload_timestamp = func.now()

        db.commit()
        db.refresh(genome)
    finally:
        db.close()
    return genome
# @router.post("/checkDuplicate")
# async def checkDuplicate(data: FileCheck):

#     #check trong database
#     user_id = data.user_id
#     filename = data.name.split(".")[0]
#     print(f"Received data: {data}")  # Debug log
#     print(f"Name: {data.name}")     

#     nonmdFA = f"nmd_{user_id}_{filename}.fa"
#     final_path = os.path.join(DATA_DIR, nonmdFA)

#     tempFA = f"temp_nmd_{user_id}_{filename}.fa"
#     temp_path = os.path.join(DATA_DIR, tempFA)
#     if os.path.exists(final_path) or os.path.exists(temp_path):
#         return JSONResponse(content={"duplicate": True})
#     return JSONResponse(content={"duplicate": False})

def createMMRegion(anno_file):
    # Lấy (.gff, .gtf, .gff3)
    ext = os.path.splitext(anno_file)[1]

    pt = anno_file.split(".")[0]

    # Tạo tên file đầu ra cho exon và gene
    exon_out = f"{pt}_exons.sorted{ext}"
    gene_out = f"{pt}_genes.sorted{ext}"


    exon_out = os.path.join(DATA_DIR, exon_out)
    gene_out = os.path.join(DATA_DIR, gene_out)

    anno_file = os.path.join(DATA_DIR, anno_file)

    # Lệnh trích xuất exon
    cmd_exon = f"awk 'BEGIN{{OFS=\"\\t\"}} $3 == \"exon\"' {anno_file} | sort -k1,1V -k4,4n > {exon_out}"
    # Lệnh trích xuất gene
    cmd_gene = f"awk 'BEGIN{{OFS=\"\\t\"}} $3 == \"gene\"' {anno_file} | sort -k1,1V -k4,4n > {gene_out}"

    try:
        # Chạy lần lượt hai lệnh
        subprocess.run(cmd_exon, shell=True, check=True, executable='/bin/bash')
        subprocess.run(cmd_gene, shell=True, check=True, executable='/bin/bash')
        print(f"Tạo thành công:\n  {exon_out}\n  {gene_out}")
    except subprocess.CalledProcessError as e:
        print(" Lỗi khi xử lý:", e)


class createData(BaseModel):
    session_id: str
    user_id: int
    fasta_file_name: str
    annotation_file_name: str
    display_id: str

@router.post("/createDataForNonModel")
async def createDataForNonModel(request_fe: Request, data: createData):

    ip = request_fe.client.host
    queue_name = 'nonModel'
    redis_key = f"running:{queue_name}:{ip}"

   #------------------------------------------------
    from .tasks import uploadNonModel_celery
    user_id = data.user_id
    fa_name = data.fasta_file_name
    anno_name = data.annotation_file_name   

    try:
        # Dùng apply_async để truyền kwargs và tên queue
        task = uploadNonModel_celery.apply_async(
            kwargs={
                'redis_key': redis_key,
                'session_id': data.session_id,
                'user_id': user_id,
                'fa_name': fa_name,
                'anno_name': anno_name,
                'display_id': data.display_id
            },
            queue=queue_name
        )
        task_id = task.id

    except Exception as e:
        await redis_client_fq.decr(redis_key)
        raise HTTPException(
            status_code=503, 
            detail=f"Không thể xếp hàng task: {e}"
        )

    update_data = GenomeUpdate(display_id=data.display_id, status="In queue", log="", task_queue_id=task_id)
    update_genome_status(update_data)
    return



def _chunk_dir(session_id: str, file_id: str, user_id: int, file_type: Literal["fasta", "annotation"]) -> str:
    type_folder = "fa" if file_type == "fasta" else "anno"
    return os.path.join(TMP_DIR, str(user_id), session_id, type_folder, file_id)

def _chunk_path(session_id: str, file_id: str, user_id: int, index: int, file_type: Literal["fasta", "annotation"]) -> str:
    return os.path.join(_chunk_dir(session_id, file_id, user_id, file_type), f"{CHUNK_PREFIX}{index:010d}")


def _final_paths(session_id: str, user_id: int, filename: str, file_type: Literal["fasta", "annotation"]) -> str:
    name_no_ext = filename.rsplit(".", 1)[0]
    if file_type == "fasta":
        out_name = f"{session_id}_temp_nmd_{user_id}_{name_no_ext}.fa"
    else:
        out_name = f"{session_id}_temp_nmd_{user_id}_{filename}"
    return os.path.join(DATA_DIR, out_name)

@router.post("/upload-chunk")
async def upload_chunk(
    session_id: str = Form(...),
    chunk: UploadFile = File(...,),#chunk data
    file_id: str = Form(...,), # id chunk
    filename: str = Form(...,),# ten file goc
    index: int = Form(...,),# chunk id
    total: int = Form(...,), # tong so chunk
    file_type: Literal["fasta", "annotation"] = Form(...), #file type
    user_id: int = Form(...,)
):
    
    print(1111111)
    if total <= 0 or index < 0 or index >= total:
        raise HTTPException(status_code=400, detail="chunk j day???")

    print(222222222)
    cdir = _chunk_dir(session_id, file_id, user_id, file_type)
    os.makedirs(cdir, exist_ok=True)

    #tai chunk xuong
    chunk_path = _chunk_path(session_id, file_id, user_id, index, file_type)
    async with aiofiles.open(chunk_path, "wb") as f:
        while True:
            data = await chunk.read(1024 * 1024)  # doc tu tu
            if not data:
                break
            await f.write(data)
    print(333333333)

    return {"ok": True, "received_index": index}

@router.post("/merge-chunks")
async def merge_chunks(
    session_id: str = Body(...),
    file_id: str = Body(...),
    filename: str = Body(...),
    total: int = Body(...),
    file_type: Literal["fasta", "annotation"] = Body(...),
    user_id: int = Body(...)
):
    cdir = _chunk_dir(session_id, file_id, user_id, file_type)
    if not os.path.isdir(cdir):
        raise HTTPException(status_code=404, detail="error1")

    out_path = _final_paths(session_id, user_id, filename, file_type)

    if os.path.exists(out_path):
        raise HTTPException(status_code=400, detail="error2")

    #gop file
    async with aiofiles.open(out_path, "wb") as out_f:
        for i in range(total):
            cpath = _chunk_path(session_id, file_id, user_id, i, file_type)
            if not os.path.exists(cpath):
                raise HTTPException(status_code=400, detail=f"Thiếu chunk {i}/{total}")
            async with aiofiles.open(cpath, "rb") as cf:
                while True:
                    data = await cf.read(1024 * 1024)
                    if not data:
                        break
                    await out_f.write(data)

    # Xóa thư mục chunk sau khi merge thành công
    shutil.rmtree(cdir, ignore_errors=True)
    return {"ok": True, "output": out_path}


def cleanup_temp_files_by_session(user_id, session_id):
    
    user_type_dir = os.path.join(TMP_DIR, str(user_id), session_id)
    if os.path.isdir(user_type_dir):
        shutil.rmtree(user_type_dir, ignore_errors=True)
        print(f" Cleaned folder for user {user_id}, session_id: {session_id}")
  
    cleanup_pattern = os.path.join(DATA_DIR, f"*{session_id}*")
    files_to_delete = glob.glob(cleanup_pattern)

    for item_path in files_to_delete:
        try:
            if os.path.isfile(item_path):
                os.remove(item_path)
                print(f"--- deleted: {os.path.basename(item_path)}")
        except Exception as e:
            print(f" Error when delete {item_path}: {e}")

    return
    

@router.post("/add_new_genome_to_db")
async def newGenomeSign(request_fe: Request, new: GenomeCreate, db: Session = Depends(get_db)):

    ip = request_fe.client.host
    queue_name = 'nonModel'
    redis_key = f"running:{queue_name}:{ip}"

    current = await redis_client_fq.incr(redis_key)
    await redis_client_fq.expire(redis_key, 3600 * 24)

    print(current)
    if current > MAX_CONCURRENT_TASKS:
        await redis_client_fq.decr(redis_key) 
        raise HTTPException(
            status_code=429, 
            detail=f"chỉ có thể có {MAX_CONCURRENT_TASKS} task non-model chạy cùng lúc từ cùng một IP, vui lòng thử lại hoặc cancel task khác đang chạy"
        )
    #------------------------------------------------
    try:
        display_id = generate_user_display_id()
        new_genome = Genome(gname=new.gname,
                        kbstorage=new.kbstorage,
                        status=new.status,
                        owner_id=new.owner_id,
                        log = new.log,
                        upload_id=new.session_id,
                        id_for_user_display=display_id
                    )
        db.add(new_genome)
        db.commit()
    except Exception as e:
        await redis_client_fq.decr(redis_key)
        cleanup_temp_files_by_session(new.owner_id, new.session_id)
        raise HTTPException(status_code=500, detail=f"Lỗi khi thêm genome mới: {e}")
    return {"msg": "genome dc tao", "display_id": display_id}

class queryAllGenome(BaseModel):
    owner_id: int

@router.post("/getAllGenomeForUser")
def getAllGenomeForUser(new: queryAllGenome, db: Session = Depends(get_db)):
    genomes = db.query(Genome).filter(Genome.owner_id == new.owner_id).all()
    return genomes

class cleanTempQuery(BaseModel):
    session_id: Optional[str] = None
    owner_id: int
    file_type: Literal["fasta", "annotation"]
    file_name: str

@router.post("/cleanTmp")
def clean_user_tmp(data: cleanTempQuery):

    user_id = data.owner_id
    file_type = data.file_type
    file_name = data.file_name.split('.')[0]
    session_id = data.session_id

    type_folder = "fa" if file_type == "fasta" else "anno"
    user_type_dir = os.path.join(TMP_DIR, str(user_id), session_id, type_folder)
    if os.path.isdir(user_type_dir):
        shutil.rmtree(user_type_dir, ignore_errors=True)
        print(f" Cleaned {type_folder} folder for user {user_id}: {user_type_dir}")
  
    file_type_map = {
        "fasta": [".fna", ".fa"],
        "annotation": [".gff", ".gtf", ".gff3"]
    }

    extensions = file_type_map.get(file_type, [])

    for ext in extensions:
        temp_pattern = os.path.join(DATA_DIR, f"{session_id}_temp_nmd_{user_id}_{file_name}*{ext}")
        temp_items = glob.glob(temp_pattern)
        
        for item_path in temp_items:
            try:
                if os.path.isfile(item_path):
                    os.remove(item_path)
                    print(f" Deleted temp file: {os.path.basename(item_path)}")
            except Exception as e:
                print(f" Error deleting {item_path}: {e}")
    return {"data": 200}

@router.post("/cleanTmplea")
def clean_when_leave(data: dict = Body(...)):
    user_id = data["owner_id"]
    file_type = data["file_type"]
    file_name = data["file_name"].split('.')[0]
    session_id = data.get("session_id", None)

    type_folder = "fa" if file_type == "fasta" else "anno"
    user_type_dir = os.path.join(TMP_DIR, str(user_id), session_id, type_folder)
    if os.path.isdir(user_type_dir):
        shutil.rmtree(user_type_dir, ignore_errors=True)
        print(f" Cleaned {type_folder} folder for user {user_id}: {user_type_dir}")
  
    file_type_map = {
        "fasta": [".fna", ".fa"],
        "annotation": [".gff", ".gtf", ".gff3"]
    }

    extensions = file_type_map.get(file_type, [])

    for ext in extensions:
        temp_pattern = os.path.join(DATA_DIR, f"{session_id}_temp_nmd_{user_id}_{file_name}*{ext}")
        temp_items = glob.glob(temp_pattern)
        
        for item_path in temp_items:
            print(item_path)
            try:
                if os.path.isfile(item_path):
                    os.remove(item_path)
                    print(f" Deleted temp file: {os.path.basename(item_path)}")
            except Exception as e:
                print(f" Error deleting {item_path}: {e}")
    return {"data": 200}
    
class DeleteTask(BaseModel):
    task_queue_id: str
    user_id: int
    faname: str
    upload_id: Optional[str] = None

@router.post("/deleteTask")
def deleteTask(data: DeleteTask, request: Request):

    print("aduvjp")
    from .tasks import celery
    tid = data.task_queue_id
    user_id = data.user_id
    faname = data.faname
    session_id = data.upload_id

    client_ip = request.client.host
    queue_name = 'nonModel'
    redis_key = f"running:{queue_name}:{client_ip}"

    flag_key = f"task_decr_flag:{user_id}:{tid}"

    task_result = AsyncResult(tid, app=celery)

    print(task_result.state)
    
    if task_result.state in ["SUCCESS", "FAILURE", "REVOKED"]:
        #xoa khoi database dong ma co taskqueueid = tid
        print("task nay se ko dc queue quan tam")
        return {"code": "done", "message": "ko the huy task"}
    
    if task_result.state == "PENDING":
        update_data = GenomeUpdate(gname=faname.split(".")[0], owner_id=user_id, status="Cancelled", log="Task was pending and revoked before execution")
        update_genome_status(update_data)
    
    celery.control.revoke(tid, terminate=True, signal='SIGINT')

    cuop_duoc_co = redis_client_fq.set(flag_key, "api", ex=7200, nx=True)
    prefixes = [ f"nmd_{user_id}_{faname.split('.')[0]}", f"{session_id}_temp_nmd_{user_id}_{faname.split('.')[0]}"]
    for prefix in prefixes:
        pattern = os.path.join(DATA_DIR, f"{prefix}*")
        for file_path in glob.glob(pattern):
            try:
                os.remove(file_path)
            except IsADirectoryError:
                shutil.rmtree(file_path, ignore_errors=True)
            except Exception as e:
                print(f"lỗi khi delete {file_path}: {e}")
    if cuop_duoc_co:
        final_value = get_and_decr_redis(redis_client_fq_celery, redis_key)
        print(f"Task hoàn thành. Key đã được giảm về {final_value}")
        return {"message": "API đã hủy và decr"}
    else:
        return {"message": "task se tự decr"}

class removeGenomeRequest(BaseModel):
    gname: str
    owner_id: int
    upload_id: Optional[str] = None

    
@router.post("/removeGenome")
def removeGenome(data: removeGenomeRequest):
    db = SessionLocal()
    try:
        genome = db.query(Genome).filter(
                (Genome.owner_id == data.owner_id) & (Genome.gname == data.gname)
            ).first()
        if not genome:
            raise HTTPException(status_code=404, detail="Genome không tồn tại")

        fasta_hash = genome.id_use_for_us_fasta
        gff3_hash = genome.id_use_for_us_gff3

        db.delete(genome)
        db.commit()

        # Reference counting: only delete physical files if no other genome references them
        if fasta_hash:
            ref_count = db.query(Genome).filter(Genome.id_use_for_us_fasta == fasta_hash).count()
            if ref_count == 0:
                shared_fasta_path = os.path.join(DATA_DIR, f"fasta_{fasta_hash}")
                if os.path.isdir(shared_fasta_path):
                    shutil.rmtree(shared_fasta_path, ignore_errors=True)
                    print(f"🗑️ Deleted shared fasta dir: {shared_fasta_path}")

        if gff3_hash:
            ref_count = db.query(Genome).filter(Genome.id_use_for_us_gff3 == gff3_hash).count()
            if ref_count == 0:
                shared_anno_path = os.path.join(DATA_DIR, f"anno_{gff3_hash}")
                if os.path.isdir(shared_anno_path):
                    shutil.rmtree(shared_anno_path, ignore_errors=True)
                    print(f"🗑️ Deleted shared anno dir: {shared_anno_path}")

        # Also clean up any legacy nmd_ files if they exist
        prefixes = [ f"nmd_{data.owner_id}_{data.gname.split('.')[0]}", f"{data.upload_id}_temp_nmd_{data.owner_id}_{data.gname.split('.')[0]}"]
        for prefix in prefixes:
            pattern = os.path.join(DATA_DIR, f"{prefix}*")
            for file_path in glob.glob(pattern):
                try:
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path, ignore_errors=True)
                except Exception:
                    pass

        return {"message": "genome da dc xoa"}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Lỗi khi xóa genome: {e}")
    finally:
        db.close()


@router.get("/genome-status/{display_id}/{owner_id}")
def get_genome_status(display_id: str, owner_id: int, db: Session = Depends(get_db)):

    try:
        genome = db.query(Genome).filter(
            (Genome.owner_id == owner_id) & (Genome.id_for_user_display == display_id)
        ).first()
        
        if not genome:
            return {
                "status": "not_found",
                "log": "Genome không tồn tại",
                "gw_state": "unknown",
                "gname": display_id,
                "owner_id": owner_id
            }
        
        return {
            "status": genome.status,
            "log": genome.log,
            "gw_state": genome.gw_state,
            "gname": genome.gname,
            "owner_id": genome.owner_id,
            "task_queue_id": genome.task_queue_id,
            "created_at": genome.kbstorage,
            "id_for_user_display": genome.id_for_user_display,
            "id_use_for_us_fasta": genome.id_use_for_us_fasta,
            "id_use_for_us_gff3": genome.id_use_for_us_gff3,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error while taking status: {str(e)}")

@router.get("/getGenomeMetadata/{display_id}")
def get_genome_metadata(display_id: str, db: Session = Depends(get_db)):
    """Lấy thông tin cơ bản của một genome dựa trên Display ID."""
    genome = db.query(Genome).filter(Genome.id_for_user_display == display_id).first()
    if not genome:
        raise HTTPException(status_code=404, detail="Genome không tồn tại")
    
    return {
        "gname": genome.gname,
        "status": genome.status,
        "gw_state": genome.gw_state,
        "kbstorage": genome.kbstorage,
        "upload_timestamp": genome.upload_timestamp,
        "display_id": genome.id_for_user_display
    }
