## 0\. Clone Dự án

Thực hiện lệnh sau để nhân bản kho lưu trữ (repository) từ Git và điều hướng vào thư mục dự án:

```bash
git clone https://github.com/xcomputing-vicrispr/Backend.git
```

## 1\. Cấu hình Mạng Docker

Tạo một mạng Docker tùy chỉnh để cho phép giao tiếp nội bộ giữa các container dịch vụ (Backend, Frontend và các dịch vụ phụ trợ như Database, Celery,..).

```bash
# Tạo mạng Docker dùng chung
docker network create vicrispr-net
```

## 2\. Khởi động Dịch vụ với Docker Compose

Sử dụng `docker-compose` để xây dựng (build) các image cần thiết và khởi động tất cả các container dịch vụ được định nghĩa trong tệp cấu hình (`docker-compose.yml`).

```bash
docker-compose up --build 
```

> **Lưu ý:** Chạy **một lần** cho cả thiết lập backend và frontend.

## 3\. Xác minh Truy cập Dịch vụ

Sau khi khởi động, ứng dụng sẽ có thể truy cập được qua cổng đã cấu hình.

  * **Địa chỉ truy cập:** `http://localhost:8098`

> **Quan trọng:** Cần một khoảng thời gian nhất định (vài phút) để tất cả các thành phần dịch vụ (Web Server, Database, Celery worker) khởi động hoàn tất và ổn định. Vui lòng chờ đợi cho đến khi server hoàn tất quá trình khởi tạo trước khi truy cập.
