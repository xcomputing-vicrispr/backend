import pandas as pd

# --- BẠN CẦN THAY ĐỔI CÁC GIÁ TRỊ NÀY ---

FILE1_NAME = 'mine24.csv'  # Tên file 1 của bạn
FILE2_NAME = 'chopchop24.csv'  # Tên file 2 của bạn

# Tên cột chứa 'seq' trong mỗi file
COLUMN_NAME_FILE1 = 'sgRNA_seq' 
COLUMN_NAME_FILE2 = 'sgRNA_seq' 

# --- HẾT PHẦN THAY ĐỔI ---

try:
    # 1. Đọc 2 file CSV
    df1 = pd.read_csv(FILE1_NAME)
    df2 = pd.read_csv(FILE2_NAME)

    # 2. Kiểm tra xem các cột bạn nhập có tồn tại trong file không
    if COLUMN_NAME_FILE1 not in df1.columns:
        print(f"Lỗi: Không tìm thấy cột '{COLUMN_NAME_FILE1}' trong file {FILE1_NAME}")
        print(f"Các cột có trong file 1 là: {list(df1.columns)}")
        
    elif COLUMN_NAME_FILE2 not in df2.columns:
        print(f"Lỗi: Không tìm thấy cột '{COLUMN_NAME_FILE2}' trong file {FILE2_NAME}")
        print(f"Các cột có trong file 2 là: {list(df2.columns)}")
        
    else:
        # 3. Lấy ra các chuỗi (sequence) và đưa vào set để tìm kiếm nhanh
        # (Giữ nguyên cấu trúc của bạn)
        seqs_file1 = set(df1[COLUMN_NAME_FILE1].dropna().unique())
        seqs_file2 = set(df2[COLUMN_NAME_FILE2].dropna().unique())

        # 4. Tìm các seq có trong file 1 NHƯNG KHÔNG có trong file 2
        # (Giữ nguyên cấu trúc của bạn)
        missing_seqs = seqs_file1.difference(seqs_file2)

        # 5. In kết quả
        if not missing_seqs:
            print("Chúc mừng! Tất cả sgRNA từ file 1 đều có mặt trong file 2.")
        else:
            print(f"Tìm thấy {len(missing_seqs)} sgRNA từ file 1 KHÔNG có trong file 2:")
            
            # --- PHẦN THÊM VÀO ---
            # Bây giờ, dùng set 'missing_seqs' để lọc lại df1 và lấy các dòng tương ứng
            # .isin() kiểm tra xem giá trị trong cột có nằm trong set 'missing_seqs' không
            missing_rows_df = df1[df1[COLUMN_NAME_FILE1].isin(missing_seqs)].copy()
            
            # Nếu file 1 có thể có nhiều seq trùng lặp, ta chỉ giữ lại dòng đầu tiên
            missing_rows_df = missing_rows_df.drop_duplicates(subset=[COLUMN_NAME_FILE1], keep='first')
            # --- HẾT PHẦN THÊM VÀO ---


            # Kiểm tra xem file 1 có đủ 3 cột không
            if len(df1.columns) < 3:
                print("--------------------------------------------------")
                print(f"(File 1 có ít hơn 3 cột, chỉ có thể in cột '{COLUMN_NAME_FILE1}')")
                # Lặp qua set gốc như cũ
                for seq in missing_seqs:
                    print(seq)
            else:
                # Lấy tên của cột thứ 2 và thứ 3 (theo index 1 và 2)
                col2_name = df1.columns[1]
                col3_name = df1.columns[2]
                
                print("--------------------------------------------------")
                # In tiêu đề
                print(f"{COLUMN_NAME_FILE1.ljust(25)} | {col2_name.ljust(15)} | {col3_name.ljust(15)}")
                print("-" * 60) # In dòng kẻ ngang

                # Lặp qua DataFrame đã lọc (thay vì lặp qua set)
                for index, row in missing_rows_df.iterrows():
                    print(f"{str(row[COLUMN_NAME_FILE1]).ljust(25)} | {str(row[col2_name]).ljust(15)} | {str(row[col3_name]).ljust(15)}")


except FileNotFoundError as e:
    print(f"Lỗi: Không tìm thấy file. Vui lòng kiểm tra lại tên file: {e.filename}")
except Exception as e:
    print(f"Đã xảy ra lỗi ngoài dự kiến: {e}")