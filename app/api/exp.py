import nmslib

# Tạo index với khoảng cách Levenshtein cho chuỗi
index = nmslib.init(
    method='hnsw',
    space='leven',  # hoặc 'levenshtein' nếu build có
    data_type=nmslib.DataType.OBJECT_AS_STRING
)

# Dữ liệu mẫu (chuỗi DNA)
data = ["ACTG", "ACCG", "ATCG", "GCTA"]

# Thêm dữ liệu vào index
index.addDataPointBatch(data)

# Xây index
index.createIndex({'post': 2}, print_progress=True)

# Query thử
query = "ACGG"
ids, distances = index.knnQuery(query, k=2)
print(ids, distances)
