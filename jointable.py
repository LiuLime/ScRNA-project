"""join the correlation table and log10p value table together
@ Wang Zekun 2024/4/11
"""


import pymysql.cursors
import pyarrow as pa
import pandas as pd

# 数据库连接参数
host = 'localhost'
user = 'root'

database = 'scrna_mt'
# database = 'scrna_mc'
# database = 'scrna_4060t'
# database = 'scrna_4060c'

tables = ['s10g1t10', 's10g1t30', 's11g2t10', 's11g2t11', 's11g2t17', 's11g2t32', 's11g2t33', 's11g2t34', 's11g2t35',
          's11g2t36', 's11g2t37', 's11g2t5', 's11g2t6', 's12g2t14', 's12g2t21', 's12g2t38', 's12g2t40', 's12g2t41',
          's12g2t5', 's1g1t1', 's2g2t2', 's2g2t3', 's3g2t4', 's3g2t5', 's3g2t6', 's4g2t7', 's5g1t5', 's6g2t8', 's7g1t8',
          's8g2t10', 's8g2t11', 's8g2t16', 's8g2t17', 's8g2t18', 's8g2t20', 's8g2t22', 's8g2t24', 's8g2t25', 's8g2t26',
          's8g2t29', 's8g2t5', 's8g2t6', 's8g2t9', 's9g2t10', 's9g2t30']
# tables = ['s10g1t10c10', 's10g1t10c14', 's10g1t10c23', 's10g1t10c6', 's10g1t10c7', 's10g1t10c8', 's10g1t10c9', 's10g1t30c10', 's10g1t30c20', 's10g1t30c6', 's10g1t30c7', 's10g1t30c8', 's10g1t30c9', 's11g2t10c6', 's11g2t10c7', 's11g2t10c9', 's11g2t11c10', 's11g2t11c17', 's11g2t11c19', 's11g2t11c20', 's11g2t11c22', 's11g2t11c23', 's11g2t11c28', 's11g2t11c29', 's11g2t11c7', 's11g2t11c9', 's11g2t17c11', 's11g2t17c20', 's11g2t17c23', 's11g2t17c7', 's11g2t17c9', 's11g2t33c27', 's11g2t33c8', 's11g2t34c20', 's11g2t34c27', 's11g2t34c6', 's11g2t34c7', 's11g2t34c9', 's11g2t35c10', 's11g2t35c23', 's11g2t35c26', 's11g2t35c27', 's11g2t35c28', 's11g2t35c29', 's11g2t35c38', 's11g2t35c7', 's11g2t35c8', 's11g2t35c9', 's11g2t36c10', 's11g2t36c20', 's11g2t36c23', 's11g2t36c27', 's11g2t36c6', 's11g2t36c7', 's11g2t36c9', 's11g2t37c9', 's11g2t5c10', 's11g2t5c2', 's11g2t5c20', 's11g2t5c23', 's11g2t5c45', 's11g2t5c9', 's11g2t6c10', 's11g2t6c20', 's11g2t6c26', 's11g2t6c27', 's11g2t6c29', 's11g2t6c7', 's11g2t6c8', 's11g2t6c9', 's12g2t14c13', 's12g2t14c16', 's12g2t14c20', 's12g2t14c25', 's12g2t14c6', 's12g2t21c13', 's12g2t21c16', 's12g2t21c2', 's12g2t21c20', 's12g2t21c24', 's12g2t21c3', 's12g2t21c4', 's12g2t38c25', 's12g2t38c31', 's12g2t38c4', 's12g2t40c13', 's12g2t40c14', 's12g2t40c17', 's12g2t40c20', 's12g2t40c25', 's12g2t40c4', 's12g2t41c13', 's12g2t41c19', 's12g2t41c25', 's12g2t41c36', 's12g2t5c2', 's12g2t5c38', 's12g2t5c4', 's12g2t5c6', 's1g1t1c1', 's1g1t1c2', 's1g1t1c3', 's1g1t1c6', 's2g2t2c10', 's2g2t2c11', 's2g2t2c2', 's2g2t2c7', 's2g2t2c9', 's2g2t3c10', 's2g2t3c12', 's2g2t3c7', 's2g2t3c9', 's3g2t4c13', 's3g2t4c14', 's3g2t4c15', 's3g2t4c16', 's3g2t4c17', 's3g2t4c2', 's3g2t4c4', 's3g2t4c8', 's3g2t5c10', 's3g2t5c13', 's3g2t5c14', 's3g2t5c17', 's3g2t5c19', 's3g2t5c20', 's3g2t5c21', 's3g2t5c22', 's3g2t5c23', 's3g2t5c4', 's3g2t5c7', 's3g2t5c8', 's3g2t5c9', 's3g2t6c10', 's3g2t6c17', 's3g2t6c20', 's3g2t6c23', 's3g2t6c26', 's3g2t6c27', 's3g2t6c7', 's3g2t6c8', 's3g2t6c9', 's4g2t7c19', 's4g2t7c28', 's4g2t7c29', 's4g2t7c3', 's4g2t7c4', 's4g2t7c8', 's5g1t5c10', 's5g1t5c13', 's5g1t5c15', 's5g1t5c19', 's5g1t5c2', 's5g1t5c20', 's5g1t5c22', 's5g1t5c24', 's5g1t5c30', 's5g1t5c4', 's5g1t5c6', 's5g1t5c7', 's5g1t5c9', 's6g2t8c10', 's6g2t8c14', 's6g2t8c27', 's6g2t8c31', 's6g2t8c6', 's6g2t8c7', 's6g2t8c8', 's6g2t8c9', 's7g1t8c10', 's7g1t8c27', 's7g1t8c4', 's7g1t8c6', 's7g1t8c7', 's7g1t8c8', 's7g1t8c9', 's8g2t10c23', 's8g2t10c29', 's8g2t10c31', 's8g2t10c32', 's8g2t10c7', 's8g2t10c8', 's8g2t10c9', 's8g2t11c10', 's8g2t11c16', 's8g2t11c17', 's8g2t11c19', 's8g2t11c20', 's8g2t11c22', 's8g2t11c23', 's8g2t11c7', 's8g2t11c9', 's8g2t18c10', 's8g2t18c12', 's8g2t18c23', 's8g2t18c7', 's8g2t18c8', 's8g2t18c9', 's8g2t20c10', 's8g2t20c13', 's8g2t20c14', 's8g2t20c16', 's8g2t20c24', 's8g2t20c4', 's8g2t20c43', 's8g2t20c8', 's8g2t22c13', 's8g2t22c19', 's8g2t22c4', 's8g2t22c7', 's8g2t22c8', 's8g2t25c14', 's8g2t25c23', 's8g2t25c24', 's8g2t25c26', 's8g2t25c4', 's8g2t25c7', 's8g2t25c8', 's8g2t26c13', 's8g2t26c19', 's8g2t26c27', 's8g2t26c8', 's8g2t29c13', 's8g2t29c2', 's8g2t29c24', 's8g2t29c4', 's8g2t29c43', 's8g2t29c8', 's8g2t29c9', 's8g2t5c13', 's8g2t5c16', 's8g2t5c19', 's8g2t5c21', 's8g2t5c22', 's8g2t5c24', 's8g2t5c28', 's8g2t5c4', 's8g2t5c6', 's8g2t5c9', 's8g2t6c10', 's8g2t6c20', 's8g2t6c21', 's8g2t6c22', 's8g2t6c23', 's8g2t6c31', 's8g2t6c45', 's8g2t6c7', 's8g2t6c8', 's8g2t6c9', 's8g2t9c13', 's8g2t9c22', 's8g2t9c24', 's8g2t9c4', 's8g2t9c7', 's8g2t9c9', 's9g2t10c10', 's9g2t10c14', 's9g2t10c23', 's9g2t10c44', 's9g2t10c6', 's9g2t10c7', 's9g2t10c8', 's9g2t10c9', 's9g2t30c10', 's9g2t30c20', 's9g2t30c22', 's9g2t30c8', 's9g2t30c9']
# tables = ['s12g2t41', 's3g2t4', 's3g2t5', 's5g1t5', 's6g2t8', 's7g1t8', 's8g2t10', 's8g2t11', 's8g2t20', 's8g2t26']
# tables = ['s12g2t41c25', 's3g2t4c13', 's3g2t4c16', 's3g2t4c17', 's3g2t4c2', 's3g2t4c8', 's3g2t5c10', 's3g2t5c14', 's3g2t5c19', 's3g2t5c20', 's3g2t5c21', 's3g2t5c23', 's3g2t5c4', 's3g2t5c7', 's3g2t5c8', 's3g2t5c9', 's5g1t5c10', 's5g1t5c13', 's5g1t5c15', 's5g1t5c19', 's5g1t5c2', 's5g1t5c20', 's5g1t5c22', 's5g1t5c24', 's5g1t5c30', 's5g1t5c4', 's5g1t5c6', 's5g1t5c7', 's5g1t5c9', 's6g2t8c10', 's6g2t8c14', 's6g2t8c27', 's6g2t8c31', 's6g2t8c6', 's6g2t8c7', 's6g2t8c8', 's6g2t8c9', 's7g1t8c10', 's7g1t8c27', 's7g1t8c4', 's7g1t8c6', 's7g1t8c7', 's7g1t8c8', 's7g1t8c9', 's8g2t10c20', 's8g2t10c23', 's8g2t10c29', 's8g2t10c7', 's8g2t10c9', 's8g2t11c16', 's8g2t11c7', 's8g2t20c10', 's8g2t20c13', 's8g2t20c14', 's8g2t20c16', 's8g2t20c24', 's8g2t20c4', 's8g2t20c43', 's8g2t20c8']

# 初始化 Arrow Table Schema
# 需要根据实际查询结果的列类型来定义，这里只是一个示例
schema = pa.schema([
    ('gene1', pa.string()),
    ('log10p1', pa.float64()),
    ('gene2', pa.string()),
    ('log10p2', pa.float64()),
    ('cor_pearson.binMean', pa.float64()),

    # 根据你的数据添加更多列
])

# 连接数据库
connection = pymysql.connect(host=host, user=user, database=database,
                             cursorclass=pymysql.cursors.SSCursor)

try:
    cursor = connection.cursor()
    for table_name in tables:
        cursor.execute(
            f"""
           SELECT a.gene1,
               a.`log10p`      AS log10p1,
               a.gene2,
               log10p.`log10p` AS log10p2,
               a.`cor_pearson.binMean`
            FROM (SELECT {table_name}.`gene1`,
                         {table_name}.`gene2`,
                         {table_name}.`cor_pearson.binMean`,
                         log10p.log10p

                FROM {table_name}
               INNER JOIN log10p ON ({table_name}.gene1 = log10p.`gene` and log10p.`abbr_id` = '{table_name}')) AS a
            INNER JOIN log10p ON (a.gene2 = log10p.`gene` and log10p.`abbr_id` = '{table_name}')
            """
        )
        # 初始化 Arrow 文件写入器
        with pa.OSFile(f'./{database}/{table_name}.arrow', 'wb') as sink:
            with pa.ipc.new_file(sink, schema) as writer:

                # 分批读取和处理数据
                while True:
                    # 读取一批数据，大小可根据内存调整
                    batch_data = cursor.fetchmany(size=100000)
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
