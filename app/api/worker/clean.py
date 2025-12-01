import asyncio
import shutil
import os
from app.api.tasks import redis_client
from app.api.nonModel import TMP_DIR  # 🔹 import cấu hình path

async def cleanup_inactive_uploads():
    """Dọn rác các upload đã hết hạn trong Redis (mất mạng, đóng trình duyệt...)"""
    while True:
        for key in redis_client.keys("upload:*"):
            # Redis TTL = -2 nghĩa là key không tồn tại (đã hết hạn)
            if redis_client.ttl(key) == -2:
                try:
                    _, user_id, file_id = key.split(":")
                except ValueError:
                    continue

                path = os.path.join(TMP_DIR, f"{user_id}_{file_id}")

                if os.path.isdir(path):
                    shutil.rmtree(path, ignore_errors=True)
                    print(f"Cleaned up {path}")
        await asyncio.sleep(60)
