# # %%
# from sqlalchemy import create_engine, Column, Integer, String, ForeignKey
# from sqlalchemy.orm import relationship, sessionmaker, declarative_base
#
# Base = declarative_base()
#
#
# # 定义子表模型
# class Child(Base):
#     __tablename__ = 'child'
#
#     id = Column(Integer, primary_key=True)
#     name = Column(String(30))
#     # 假设有其他字段
#
#
# # 定义总表模型，其中包含指向子表的外键
# class Parent(Base):
#     __tablename__ = 'parent'
#
#     id = Column(Integer, primary_key=True)
#     child_id = Column(Integer, ForeignKey('child.id'))
#
#     # 通过relationship建立与Child模型的联系，这将允许我们方便地访问子表的实例
#     child = relationship("Child")
#
#
# # 创建数据库引擎
# engine = create_engine('mysql+mysqlconnector://root@localhost:3306/test')
#
# # 创建所有表
# Base.metadata.create_all(engine)
#
# # 创建会话
# Session = sessionmaker(bind=engine)
# session = Session()
#
# # 示例：创建一个Child实例和一个与之关联的Parent实例
# # new_child = Child(name ="child 1")
# # new_parent = Parent(child=new_child)
# # session.add(new_parent)
#
# new_child = [Child("Child 1"), Child("Child2")]
# new_parent=Parent(child=new_child)
# session.add_all(new_parent)
#
#
# session.commit()
#
# # 查询示例：获取所有Parent对象及其关联的Child对象
# parents = session.query(Parent).all()
# for parent in parents:
#     print(f'Parent ID: {parent.id}, Child Name: {parent.child.name}')

# %%
from database import db
from sqlalchemy import text, create_engine


# def add_id(table):
#     with engine.connect() as conn:
#         conn.execute(text(f"ALTER TABLE {table} DROP PRIMARY KEY"))
#         conn.execute(text(f"ALTER TABLE {table} ADD COLUMN id INT AUTO_INCREMENT PRIMARY KEY FIRST"))
#         conn.commit()
#
#
# db_list = ["scrna_mt", "scrna_mc", "scrna_4060t", "scrna_4060c"]
# for db_path in db_list:
#     engine = create_engine(f'mysql+mysqlconnector://root@localhost:3306/{db_path}')
#     db = db.databaseHelper(db_path)
#     db.init()
#     _, tasks = db.get_available_tasks()
#     print("available_tasks", tasks)
#     # print(db.get_table_cutoff('s8g2t20',0.7,1))
#     for task in tasks:
#         table = task["Table"]
#         add_id(table)

# %%
# from utils.log import logger
# log=logger()
# from database import db
# db_list = ["scrna_mt", "scrna_mc", "scrna_4060t", "scrna_4060c"]
# db = db.databaseHelper(db_list[2])
# db.init()
#
# log.debug(db.get_correlation_records_above_threshold('s8g2t20',0.8))

# %%
# from utils.log import logger
#
# log = logger()
# log.error("test error")
# log.debug("test debug")
#
# #%%
# import pandas as pd
# test = [{'id': 191983, 'gene': 'A2ML1', 'log10p': '4.52676230992907', 'abbr_id': 's3g2t4'},
#         {'id': 192095, 'gene': 'AC024293.1', 'log10p': '25.3048788411024', 'abbr_id': 's3g2t4'}]
#
# data=pd.DataFrame(test)

#%%
# test=[{'id': 2, 'Table': 's12g2t41c25'}, {'id': 3, 'Table': 's3g2t4c13'}]
# print(list(zip(*test[0].items())))
#%%
# import pandas as pd
#
# path = "/Users/liuyuting/WorkSpace/ScRNA_project/mysql/20230328_median_10humanStudies-20230410_40-60_5humanStudies.1ageGrpGE50cells/ageGrp_byMedian/tissue-cellType_level/log10p.csv"
# data = pd.read_csv(path, sep=",",header=0)
#
# count=0
# for idx, row in data.iterrows():
#     if row["log10p"] =='#NAME?':
#         count+=1
# print(count)

#%%
from MySQL import sql
import pandas as pd
test_p=pd.DataFrame([{"gene":"a","log10p":3,"abbr_id":"s1"},
                     {"gene":"b","log10p":33,"abbr_id":"s2"}])

s=sql('test')
with s:
    # tasks=s.get_avalible_studies()
    # print(tasks)
    # s.create_log10p_table(test_p)
    table=s.screen_by_corr_log10p('s12g2t38',1,1)

