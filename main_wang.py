import pymysql.cursors
import pyarrow as pa
import pandas as pd

# 数据库连接参数
host = '数据库地址'
user = '用户名'
password = '密码'
database = '数据库名'
tables = []
# 初始化 Arrow Table Schema
# 需要根据实际查询结果的列类型来定义，这里只是一个示例
schema = pa.schema([
    ('gene1', pa.string()),
    ('gene2', pa.string()),
    ('cor', pa.float64()),
    ('p1', pa.float64()),
    ('p2', pa.float64()),
    # 根据你的数据添加更多列
])

# 连接数据库
connection = pymysql.connect(host=host, user=user, password=password, database=database,
                             cursorclass=pymysql.cursors.SSCursor)

try:
    cursor = connection.cursor()
    for table_name in tables:
        cursor.execute(
            f"""
           SELECT a.gene1,
       a.`log10p`      AS p1,
       a.gene2,
       log10p.`log10p` AS p2,
       a.`cor_pearson.binMean`
FROM (SELECT s12g2t14.`gene1`,
             s12g2t14.`gene2`,
             s12g2t14.`cor_pearson.binMean`,
             log10p.log10p

      FROM s12g2t14
               INNER JOIN log10p ON (s12g2t14.gene1 = log10p.`gene` and log10p.`abbr_id` = 's12g2t14')) AS a
         INNER JOIN log10p ON (a.gene2 = log10p.`gene` and log10p.`abbr_id` = 's12g2t14')
            """
        )
        # 初始化 Arrow 文件写入器
        with pa.OSFile(f'{table_name}.arrow', 'wb') as sink:
            with pa.ipc.new_file(sink, schema) as writer:

                # 分批读取和处理数据
                while True:
                    # 读取一批数据，大小可根据内存调整
                    batch_data = cursor.fetchmany(size=10000)
                    if not batch_data:
                        break  # 数据读取完毕

                    # 转换为 Pandas DataFrame
                    df = pd.DataFrame(batch_data, columns=[desc[0] for desc in cursor.description])

                    # 转换 DataFrame 为 Arrow RecordBatch
                    batch = pa.RecordBatch.from_pandas(df, schema=schema, preserve_index=False)

                    # 写入 Arrow 文件
                    writer.write_batch(batch)
                    print(f"{table_name}成功保存为 Arrow 文件")

finally:
    cursor.close()
    connection.close()
