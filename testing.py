import os
import subprocess

def get_file_path(filename):
    return os.path.join("app", "data", filename)

def get_multi_fasta_hash(filepath):
    try:
        cmd = f"seqkit sort -n {filepath} 2> /dev/null | seqkit seq -w 0 | md5sum"
        
        output = subprocess.check_output(cmd, shell=True).decode()
        # Lấy phần tử đầu tiên (chính là mã MD5)
        return output.split()[0]
    except Exception as e:
        print(f"❌ Lỗi: {e}")
        return None

def main():
    file_a = "atest.fa"
    file_b = "btest.fa"
    path1 = get_file_path(file_a)
    path2 = get_file_path(file_b)

    if not os.path.exists(path1) or not os.path.exists(path2):
        print("❌ Lỗi: Không tìm thấy file.")
        return

    print(f"🧬 Đang kiểm tra Multi-FASTA: {file_a} VS {file_b}")
    
    hash_a = get_multi_fasta_hash(path1)
    hash_b = get_multi_fasta_hash(path2)

    if hash_a and hash_b:
        print(f"🔑 Hash A (Sorted): {hash_a}")
        print(f"🔑 Hash B (Sorted): {hash_b}")

        if hash_a == hash_b:
            print("\n✅ KẾT QUẢ: GIỐNG NHAU TUYỆT ĐỐI.")
            print("   (Đã kiểm tra tất cả Chromosome, bất kể thứ tự dòng hay ngắt dòng)")
        else:
            print("\n🚫 KẾT QUẢ: KHÁC NHAU.")
            print("   (Có thể do tên Chromosome khác nhau hoặc có biến dị nucleotide)")

if __name__ == "__main__":
    main()