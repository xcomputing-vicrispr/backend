from fastapi import FastAPI
import requests, os, glob, time, shutil
import asyncio


PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PARENT_DIR, "app/data")
TMP_DIR = os.path.join(PARENT_DIR, "app/data/tmp")



PATTERNS = [
    (DATA_DIR, "*_sgrna_output.fa"),
    (DATA_DIR, "*_mm_annotation.bed"),
    (DATA_DIR, "*_raw.bed"),
    (PARENT_DIR, "*.agat.log"),
]
vcp_PATTERNS = [
    (DATA_DIR, "vcp*"),
    (DATA_DIR, "*_primerINP.txt"),
]

MAX_AGE = 60 * 60  #1 tieng
MAX_AGE_vcp = 60 * 60 * 24 * 7# 7 ngay
MAX_AGE_tmp = 60 * 60 * 24 #1 ngay


async def cleanup_files():
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

        await asyncio.sleep(1000)
