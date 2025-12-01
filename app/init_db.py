import os
import sys
import psycopg2
from sqlalchemy.engine.url import make_url
from sqlalchemy.exc import OperationalError
import time 

from app.configs import get_settings
from app.database import Base, engine 
from app.models import *

def initialize_database():
    """
    Thực hiện 2 bước: 
    1. Tạo Database trống (nếu chưa có).
    2. Tạo tất cả các Tables trong Database.
    """
    settings = get_settings()
    db_url = make_url(settings.DATABASE_URL)
    DB_NAME = db_url.database
    
    # URL kết nối đến database 'postgres' mặc định để tạo DB mới
    postgres_url_parts = db_url._replace(database='postgres')
    POSTGRES_URL = str(postgres_url_parts)
    
    #tạo db trống
    conn = None
    try:
        print("--- BƯỚC 1: TẠO DATABASE TRỐNG ---")
        conn_params = {
            'host': db_url.host,
            'port': db_url.port or 5432, 
            'user': db_url.username,
            'password': db_url.password,
            'dbname': 'postgres'
        }
        conn = psycopg2.connect(**conn_params)
        conn.autocommit = True
        cursor = conn.cursor()
        
        cursor.execute(f"SELECT 1 FROM pg_database WHERE datname = '{DB_NAME}'")
        exists = cursor.fetchone()
        
        if not exists:
            print(f"Database '{DB_NAME}' chưa tồn tại. Đang tạo...")
            owner_clause = f" OWNER {db_url.username}" if db_url.username else ""
            cursor.execute(f"CREATE DATABASE {DB_NAME}{owner_clause};")
            print(f"Database '{DB_NAME}' đã được tạo thành công.")
        else:
            print(f"Database '{DB_NAME}' đã tồn tại, bỏ qua bước tạo DB.")
        
        cursor.close()

    except psycopg2.errors.OperationalError as e:
        print("LỖI KẾT NỐI POSTGRES: Đảm bảo PostgreSQL Server đang chạy và User/Password trong .env là đúng.")
        print(f"Chi tiết lỗi: {e}")
        return
    except Exception as e:
        print(f"Lỗi không xác định khi tạo DB: {e}")
        return
    finally:
        if conn:
            conn.close()
    time.sleep(1) 

    # tạo bảng
    try:
        print("\n--- BƯỚC 2: TẠO CÁC BẢNG TRONG DATABASE ---")
        # Lệnh này sẽ kết nối đến DB vừa được tạo và tạo các bảng
        Base.metadata.create_all(bind=engine)
        print("Tất cả các bảng (tables) đã được tạo thành công!")
        
    except OperationalError as e:
        print("LỖI TẠO BẢNG: Không thể kết nối hoặc tạo bảng.")
        print(f"Chi tiết lỗi: {e}")
    except Exception as e:
        print(f"Lỗi không xác định khi tạo bảng: {e}")

if __name__ == "__main__":
    initialize_database()