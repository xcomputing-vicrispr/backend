import pandas as pd

def compare_overlap(file1, file2):
    try:
        # Đọc file, lấy cột đầu tiên (index 0)
        # usecols=[0] giúp tiết kiệm bộ nhớ nếu file lớn
        df1 = pd.read_csv(file1, usecols=[0])
        df2 = pd.read_csv(file2, usecols=[0])
        
        # Chuyển đổi thành set để loại bỏ trùng lặp nội bộ và so sánh nhanh
        # Chúng ta lấy values của cột đầu tiên
        set1 = set(df1.iloc[:, 0].dropna().astype(str))
        set2 = set(df2.iloc[:, 0].dropna().astype(str))
        
        # Tính toán các tập hợp
        intersection = set1.intersection(set2)
        only_in_f1 = set1 - set2
        only_in_f2 = set2 - set1
        
        # Hiển thị kết quả
        print(f"--- Kết quả so sánh ---")
        print(f"File 1 ({file1}): {len(set1)} trình tự duy nhất.")
        print(f"File 2 ({file2}): {len(set2)} trình tự duy nhất.")
        print("-" * 30)
        print(f"Số lượng trình tự trùng nhau (Overlap): {len(intersection)}")
        
        if len(set1) > 0:
            percent_overlap_f1 = (len(intersection) / len(set1)) * 100
            print(f"Tỷ lệ trùng lặp so với file 1: {percent_overlap_f1:.2f}%")
            
        if len(set2) > 0:
            percent_overlap_f2 = (len(intersection) / len(set2)) * 100
            print(f"Tỷ lệ trùng lặp so với file 2: {percent_overlap_f2:.2f}%")
            
        print("-" * 30)
        print(f"Số lượng chỉ có trong {file1}: {len(only_in_f1)}")
        print(f"Số lượng chỉ có trong {file2}: {len(only_in_f2)}")
        
        # Tùy chọn: Lưu các trình tự trùng nhau ra file mới
        if len(intersection) > 0:
            pd.Series(list(intersection), name='overlapping_seq').to_csv("overlapping_sequences.csv", index=False)
            print(f"\nĐã lưu {len(intersection)} trình tự trùng nhau vào file 'overlapping_sequences.csv'")

    except FileNotFoundError as e:
        print(f"Lỗi: Không tìm thấy file. Vui lòng kiểm tra lại tên file.")
    except Exception as e:
        print(f"Đã xảy ra lỗi: {e}")

# Thực thi
compare_overlap("gw.csv", "sgrna_sequences_final.csv")