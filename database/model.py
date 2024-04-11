from sqlalchemy import Column, Integer, String, Numeric, Index
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


class Abstract(Base):
    __tablename__ = "abstract"

    id = Column(Integer, primary_key=True, autoincrement=True)  # 每个task的唯一id
    Table = Column(String(30))  # table name in this schema

    # type = Column(String(30))  # data type correlation data or pvalue data

    def to_dict(self):
        # return {'id': self.id, 'type': self.type, 'table_id': self.table_id}
        return {'id': self.id, 'Table': self.Table}


class Log10p(Base):
    __tablename__ = "log10p"
    id = Column(Integer, primary_key=True, autoincrement=True)
    gene = Column(String(30))
    log10p = Column(Numeric)
    abbr_id = Column(String(30))

    def to_dict(self):
        return {'id': self.id, 'gene': self.gene, 'log10p': self.log10p, 'abbr_id': self.abbr_id}


Index('idx_log10p', Log10p.gene, Log10p.log10p, Log10p.abbr_id)


# def dynamic_corr_class(table_name):
#     class table_name(Base):
#         __tablename__ = table_name
#         gene1 = Column(String(30), primary_key=True)
#         gene2 = Column(String(30), primary_key=True)
#         cor_pearson_binMean = Column("cor_pearson.binMean", Numeric)
#
#         def to_dict(self):
#             return {'gene1': self.gene1, 'gene2': self.gene2,
#                     'cor_pearson.binMean': self.cor_pearson_binMean}
#
#     return table_name

# 动态构建class的方法--->使用type()
def dynamic_corr_class(table_name):
    # Define class attributes (columns) in a dictionary
    attributes = {
        '__tablename__': table_name,
        'id': Column(Integer, primary_key=True, autoincrement=True),
        'gene1': Column(String(30)),
        'gene2': Column(String(30)),
        'cor_pearson_binMean': Column("cor_pearson.binMean", Numeric),
        'to_dict': lambda self: {'id': self.id, 'gene1': self.gene1, 'gene2': self.gene2,
                                 'cor_pearson.binMean': self.cor_pearson_binMean}
    }

    # Create a new class dynamically with a unique class name for each table
    cls = type(table_name, (Base,), attributes)

    # After the class is created, define an Index on the desired columns.
    # Note: Directly use the column attribute of the class.
    Index('idx_corr', cls.gene1, cls.gene2, cls.cor_pearson_binMean)
    return cls
