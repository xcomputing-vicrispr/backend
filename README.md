# Project Setup Guide

## 0. Clone dự án

Clone từ Git về máy:

```bash
git clone https://github.com/xcomputing-vicrispr/Backend.git
```

Đi vào thư mục dự án:

```bash
cd Backend
```
## 1. Cài đặt Miniconda

Tải và chạy installer của Miniconda, ví dụ:

```bash
bash ~/Downloads/Miniconda3-latest-Linux-x86_64.sh
```

> Cần tìm một bản Miniconda; nếu tên file khác, hãy thay lại cho đúng.
> Bạn có thể chạy file `.sh` từ thư mục bất kỳ — thư mục cài đặt sẽ được hỏi trong quá trình setup.

---

## 2. Cài đặt System Dependencies

```bash
sudo apt update
sudo apt install -y \
    build-essential pkg-config python3-apt libffi-dev libcairo2-dev \
    python-dev-is-python3 gir1.2-gtk-3.0 meson ninja-build \
    libdbus-1-dev libgirepository1.0-dev
```

---

## 3. Tạo Các Môi Trường Conda

### ➤ Tạo environment `agat`

```bash
conda env create -f agat_enviroment.yml
```

### ➤ Tạo environment `env2`

```bash
conda env create -f env2_environment.yml
```

### ➤ Update environment `base`

```bash
conda activate base
conda env update -f base_enviroment.yml
```

---

## 4. Khởi Tạo Database

Đi đến thư mục dự án:

```bash
conda activate base
cd /app
python init_db.py
```

---

## 5. Khởi Động Server

```bash
conda activate base
python -m uvicorn app.main:app --reload
```
