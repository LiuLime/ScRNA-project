"""
Save ScRNA result into MySQL database

@ Liu
"""

import mysql.connector
from mysql.connector import errorcode
from utils.log import logger

import pandas as pd


class sql:
    """
    Example
    -------
    ```
    >>> sql = sql()
    >>> with sql as sql:
    >>>     cid = sql.search("synomys", "synomys","cid","ribitol")[0]
    ```
    """

    def __init__(self, db):
        self.db = db
        self.con = None
        self.cur = None
        self.log = logger()

    def __enter__(self):
        self.con = mysql.connector.connect(host="localhost",
                                           user="root",
                                           charset='utf8mb4')

        self.cur = self.con.cursor()
        try:
            print(f"Databse {self.db} chosen")
            self.cur.execute(f"USE {self.db}")

        except mysql.connector.Error as err:
            if err.errno == errorcode.ER_BAD_DB_ERROR:
                print("Database does not exist...")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cur.close()
        self.con.close()

    def filter_data_by_rank(self):
        with open("./pickTopP.sql", 'r') as file:
            sql_script = file.read()
        try:
            self.cur.execute(sql_script)
            results = self.cur.fetchall()
            return results
        except Exception as e:
            print("filter_data_by_rank error:", e)
