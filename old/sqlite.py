import sqlite3
import pandas as pd
import utils


class sqlite:
    def __init__(self, db):
        self.db = db

    def __enter__(self):
        self.con = sqlite3.connect(self.db)
        self.cur = self.con.cursor()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.con.close()

    def load_from_df(self, file_path, name, if_exists="replace"):
        delimiter = utils.detect_delimiter(file_path)
        pd.read_csv(file_path, delimiter=delimiter).to_sql(name, self.con, if_exists=if_exists)
        print(f"load table {name}")
        self.con.commit()

    def checkTableExist(self, name):
        sql = f"SELECT name FROM sqlite_master WHERE type='table' AND name=`{name}`"
        return sql

    def createTableSimple(self, name, columns: tuple):
        if not self.checkTableExist(name):
            sql = f"CREATE TABLE {name}{columns}"
            self.cur.execute(sql)
            self.con.commit()
            print(f"create simple table {name}")
        else:
            print(f"Already exist simple table {name}")
