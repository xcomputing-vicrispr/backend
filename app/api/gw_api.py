"""
Genome-Wide API — Submit, Status, Download, History
"""
import json
import os
from fastapi import APIRouter, HTTPException, Depends
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import Optional
from sqlalchemy.orm import Session

from app.database import SessionLocal, get_db
from app.models import Genome, GWTaskMetadata
from app.api.hashing import generate_query_hash
from app.api.nonModel import GenomeUpdate, update_genome_status

router = APIRouter()

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PARENT_DIR, "data")


# ── Request / Response schemas ──

class GWSubmitRequest(BaseModel):
    display_id: str
    PAM: str = "NGG"
    sgRNA_length: int = 20
    seed_region: int = 9
    hamming_distance: int = 3
    flank_up: int = 100
    flank_down: int = 100
    emails: list[str] = []


# ── Helper ──

def update_gw_task(hashing_code: str, **fields):
    """Update gw_task_metadata fields by hashing_code."""
    from sqlalchemy import func as sa_func
    db = SessionLocal()
    try:
        task = db.query(GWTaskMetadata).filter(GWTaskMetadata.hashing_code == hashing_code).first()
        if task:
            for k, v in fields.items():
                setattr(task, k, v)
            if fields.get('state') in ('success', 'failed'):
                task.completed_at = sa_func.now()
            db.commit()
    finally:
        db.close()


# ── Endpoints ──

@router.post("/submit")
async def submit_gw_task(data: GWSubmitRequest):
    """
    Submit a genome-wide analysis task.
    Uses hashing_code for caching: same params → same hash → return cached result.
    """
    db = SessionLocal()
    try:
        # 1. Validate genome
        genome = db.query(Genome).filter(Genome.id_for_user_display == data.display_id).first()
        if not genome:
            raise HTTPException(status_code=404, detail="Genome not found")
        if genome.status != 'success':
            raise HTTPException(status_code=400, detail=f"Genome is not ready (status: {genome.status})")

        # 2. Check genome-wide lock
        if genome.gw_state == 'navailable':
            raise HTTPException(
                status_code=409,
                detail="A genome-wide task is already running on this genome. Please wait for it to finish."
            )

        # 3. Size check (same logic as original gw_faiss.py)
        if genome.kbstorage and genome.kbstorage > 60956679:
            raise HTTPException(status_code=400, detail="Genome size exceeds the limit for genome-wide analysis")

        # 4. Compute hashing_code
        query_params = {
            "type": "genome_wide",
            "display_id": data.display_id,
            "pam": data.PAM,
            "sgrna_length": data.sgRNA_length,
            "seed_region": data.seed_region,
            "hamming_distance": data.hamming_distance,
            "flank_up": data.flank_up,
            "flank_down": data.flank_down,
        }
        hashing_code = generate_query_hash(query_params)

        # 5. Cache check
        existing = db.query(GWTaskMetadata).filter(
            GWTaskMetadata.hashing_code == hashing_code
        ).first()

        if existing and existing.state not in ('failed',):
            print(f"[GW CACHE] Hit: {hashing_code}, state={existing.state}")
            return {"hashing_code": hashing_code, "cached": True, "state": existing.state}

        # 6. Cache miss — create new record
        # Clean old failed record if exists
        if existing and existing.state == 'failed':
            db.delete(existing)
            db.commit()

        # Filter valid emails
        valid_emails = [e.strip() for e in data.emails if e.strip()]

        new_task = GWTaskMetadata(
            hashing_code=hashing_code,
            genome_display_id=data.display_id,
            pam=data.PAM,
            sgrna_len=data.sgRNA_length,
            seed_region=data.seed_region,
            hamming_distance=data.hamming_distance,
            flank_up=data.flank_up,
            flank_down=data.flank_down,
            emails=json.dumps(valid_emails) if valid_emails else None,
            state="pending",
            current_phase="queued",
        )
        db.add(new_task)
        db.commit()

        # 7. Lock genome for GW processing
        update_data = GenomeUpdate(display_id=data.display_id, gw_state="navailable")
        update_genome_status(update_data)

        # 8. Dispatch Celery task
        from app.api.tasks import run_gw_pipeline
        celery_task = run_gw_pipeline.delay(
            hashing_code,
            data.display_id,
            data.PAM,
            data.sgRNA_length,
            data.seed_region,
            data.hamming_distance,
            data.flank_up,
            data.flank_down,
            valid_emails,
        )

        # 9. Save celery task id
        task_record = db.query(GWTaskMetadata).filter(
            GWTaskMetadata.hashing_code == hashing_code
        ).first()
        if task_record:
            task_record.queue_task_id = celery_task.id
            db.commit()

        return {"hashing_code": hashing_code, "cached": False, "state": "pending"}

    finally:
        db.close()


@router.get("/status/{hashing_code}")
async def get_gw_status(hashing_code: str):
    """Get current status of a genome-wide task."""
    db = SessionLocal()
    try:
        task = db.query(GWTaskMetadata).filter(
            GWTaskMetadata.hashing_code == hashing_code
        ).first()

        if not task:
            raise HTTPException(status_code=404, detail="Task not found")

        return {
            "hashing_code": task.hashing_code,
            "genome_display_id": task.genome_display_id,
            "state": task.state,
            "current_phase": task.current_phase,
            "log": task.log,
            "pam": task.pam,
            "sgrna_len": task.sgrna_len,
            "seed_region": task.seed_region,
            "hamming_distance": task.hamming_distance,
            "flank_up": task.flank_up,
            "flank_down": task.flank_down,
            "result_file": task.result_file,
            "result_count": task.result_count,
            "created_at": task.created_at.isoformat() if task.created_at else None,
            "completed_at": task.completed_at.isoformat() if task.completed_at else None,
        }
    finally:
        db.close()


@router.get("/download/{hashing_code}")
async def download_gw_result(hashing_code: str):
    """Download the result CSV file."""
    db = SessionLocal()
    try:
        task = db.query(GWTaskMetadata).filter(
            GWTaskMetadata.hashing_code == hashing_code
        ).first()

        if not task:
            raise HTTPException(status_code=404, detail="Task not found")
        if task.state != 'success':
            raise HTTPException(status_code=400, detail=f"Task is not complete (state: {task.state})")
        if not task.result_file:
            raise HTTPException(status_code=404, detail="Result file not available")

        file_path = os.path.join(DATA_DIR, task.result_file)
        if not os.path.exists(file_path):
            raise HTTPException(status_code=404, detail="Result file has been cleaned up")

        return FileResponse(
            path=file_path,
            filename=f"gw_result_{hashing_code[:12]}.csv",
            media_type="text/csv"
        )
    finally:
        db.close()


@router.get("/history/{genome_display_id}")
async def get_gw_history(genome_display_id: str):
    """Get all genome-wide tasks for a specific genome."""
    db = SessionLocal()
    try:
        tasks = db.query(GWTaskMetadata).filter(
            GWTaskMetadata.genome_display_id == genome_display_id
        ).order_by(GWTaskMetadata.created_at.desc()).all()

        return [
            {
                "hashing_code": t.hashing_code,
                "pam": t.pam,
                "sgrna_len": t.sgrna_len,
                "state": t.state,
                "current_phase": t.current_phase,
                "result_count": t.result_count,
                "created_at": t.created_at.isoformat() if t.created_at else None,
                "completed_at": t.completed_at.isoformat() if t.completed_at else None,
            }
            for t in tasks
        ]
    finally:
        db.close()
