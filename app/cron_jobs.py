from fastapi import FastAPI
import os, glob, time, shutil
import asyncio
from app.database import SessionLocal
from app.models import TaskMetadata, GWTaskMetadata
from datetime import datetime, timedelta
from sqlalchemy import delete

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PARENT_DIR, "app/data")
TMP_DIR = os.path.join(PARENT_DIR, "app/data/tmp")

PATTERNS = [
    (DATA_DIR, "*_sgrna_output.fa"),
    (DATA_DIR, "*_mm_annotation.bed"),
    (DATA_DIR, "*_raw.bed"),
    (DATA_DIR, "*temp*"),
    (PARENT_DIR, "*.agat.log"),
]
vcp_PATTERNS = [
    (DATA_DIR, "vcp*"),
    (DATA_DIR, "*_primerInp.txt"),
]

tmp_PATTERNS = []

MAX_AGE = 60 * 60  #1 tieng
MAX_AGE_vcp = 60 * 60 * 24# 1 ngay
MAX_AGE_tmp = 60 * 60 * 24 #1 ngay
MAX_AGE_GW_RESULT = 60 * 60 * 24 * 7  # 7 ngày

#da update thanh 1 task se ton tai trong 3 ngay

async def cleanup_files():

    db = SessionLocal()
    while True:
        now = time.time()
        print("bat dau truy quet")
        for folder, pattern in PATTERNS:
            path = os.path.join(folder, pattern)
            for file in glob.glob(path):
                try:
                    mtime = os.path.getmtime(file)
                    if now - mtime > MAX_AGE:
                        os.remove(file)
                        print(f"Deleted old file: {file}")
                except Exception as e:
                    print("Error deleting file:", file, e)

        for folder, pattern in vcp_PATTERNS:
            path = os.path.join(folder, pattern)
            for file in glob.glob(path):
                try:
                    mtime = os.path.getmtime(file)
                    if now - mtime > MAX_AGE_vcp:
                        os.remove(file)
                        print(f"Deleted old file: {file}")
                except Exception as e:
                    print("Error deleting file:", file, e)

            now = time.time()

            if os.path.exists(TMP_DIR):
                for name in os.listdir(TMP_DIR):
                    subfolder = os.path.join(TMP_DIR, name)
                    if os.path.isdir(subfolder):
                        try:
                            mtime = os.path.getmtime(subfolder)
                            if now - mtime > MAX_AGE_tmp:
                                shutil.rmtree(subfolder)
                                print(f"Deleted old temp folder: {subfolder}")
                        except Exception as e:
                            print(f"Error deleting folder {subfolder}: {e}")

        try: 
            limit_time = datetime.now() - timedelta(days=3)
        
            stmt = delete(TaskMetadata).where(TaskMetadata.created_at < limit_time)
            db.execute(stmt)
            db.commit()
            print(f"Deleted old task data")
        except Exception as e:
            print(f"Error deleting old task {e}")

        # --- Genome-Wide result cleanup (7-day retention) ---
        try:
            gw_limit = datetime.now() - timedelta(days=7)

            # Find expired GW tasks with result files
            expired_tasks = db.query(GWTaskMetadata).filter(
                GWTaskMetadata.created_at < gw_limit
            ).all()

            for task in expired_tasks:
                # Delete result CSV file if exists
                if task.result_file:
                    file_path = os.path.join(DATA_DIR, task.result_file)
                    if os.path.exists(file_path):
                        os.remove(file_path)
                        print(f"[GW Cleanup] Deleted result file: {file_path}")

            # Delete expired GW task records
            stmt_gw = delete(GWTaskMetadata).where(GWTaskMetadata.created_at < gw_limit)
            db.execute(stmt_gw)
            db.commit()
            print(f"[GW Cleanup] Deleted expired GW task records")

            # Also clean orphaned gw_result_*.csv files older than 7 days
            gw_pattern = os.path.join(DATA_DIR, "gw_result_*.csv")
            for f in glob.glob(gw_pattern):
                try:
                    mtime = os.path.getmtime(f)
                    if now - mtime > MAX_AGE_GW_RESULT:
                        os.remove(f)
                        print(f"[GW Cleanup] Deleted orphaned result file: {f}")
                except Exception as e:
                    print(f"[GW Cleanup] Error deleting {f}: {e}")

        except Exception as e:
            print(f"[GW Cleanup] Error: {e}")

        await asyncio.sleep(1000)
