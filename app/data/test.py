import os

def clean_sgrna_txt(input_path, output_path):
    # Kiểm tra file đầu vào
    if not os.path.exists(input_path):
        print(f"Lỗi: Không tìm thấy file {input_path}")
        return

    print(f"Bắt đầu lọc file: {input_path}")
    
    try:
        with open(input_path, 'r', encoding='utf-8') as fin, \
             open(output_path, 'w', encoding='utf-8') as fout:
            
            # 1. Đọc và giữ lại dòng tiêu đề (Header)
            header = fin.readline()
            if header:
                fout.write(header)
            
            count_removed = 0
            count_kept = 0
            
            # 2. Duyệt qua từng dòng dữ liệu còn lại
            for line in fin:
                # Nếu dòng không chứa dấu "?" thì mới ghi vào file đích
                if '?' not in line:
                    fout.write(line)
                    count_kept += 1
                else:
                    count_removed += 1
            
            print("-" * 30)
            print(f"Kết quả xử lý:")
            print(f" * Số dòng giữ lại: {count_kept}")
            print(f" * Số dòng chứa '?' đã bị loại bỏ: {count_removed}")
            print(f" * File sạch đã lưu tại: {output_path}")
            print("-" * 30)

    except Exception as e:
        print(f"Đã xảy ra lỗi khi xử lý: {e}")

# --- Cấu hình tên file ---
input_file = "ccruddi-guides.txt"  # Thay tên file .txt của bạn vào đây
output_file = "ccruddi-guides_final.txt"


import os

def extract_seq_only_to_csv(input_path, output_path):
    # Kiểm tra file đầu vào
    if not os.path.exists(input_path):
        print(f"Lỗi: Không tìm thấy file {input_path}")
        return

    print(f"Đang xử lý: {input_path}")
    
    try:
        with open(input_path, 'r', encoding='utf-8') as fin, \
             open(output_path, 'w', encoding='utf-8') as fout:
            
            # 1. Đọc và bỏ qua dòng tiêu đề (Header) của file gốc
            fin.readline()
            
            # 2. Ghi tiêu đề mới cho file CSV đầu ra
            fout.write("seq\n")
            
            count_kept = 0
            
            # 3. Duyệt qua từng dòng dữ liệu
            for line in fin:
                line = line.strip()
                if not line:
                    continue  # Bỏ qua dòng trống
                
                # Nếu dòng không chứa dấu "?"
                if '?' not in line:
                    # Tách các cột bằng dấu phẩy và lấy cột đầu tiên (index 0)
                    parts = line.split(',')
                    if len(parts) > 0:
                        seq = parts[0].strip()
                        # Loại bỏ dấu ngoặc kép nếu có (ví dụ "GTTG...")
                        seq = seq.replace('"', '').replace("'", "")
                        
                        fout.write(f"{seq}\n")
                        count_kept += 1
            
            print("-" * 30)
            print(f"Hoàn thành!")
            print(f" * Tổng số trình tự đã trích xuất: {count_kept}")
            print(f" * Kết quả lưu tại: {output_path}")
            print("-" * 30)

    except Exception as e:
        print(f"Đã xảy ra lỗi: {e}")

# --- Cấu hình ---
input_file = "ccruddi-guides_final.txt"  # File gốc của bạn
output_file = "sgrna_sequences_final.csv"

if __name__ == "__main__":
    extract_seq_only_to_csv(input_file, output_file)

# if __name__ == "__main__":
#     clean_sgrna_txt(input_file, output_file)