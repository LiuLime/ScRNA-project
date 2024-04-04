"""
Save ScRNA result into MySQL database

@ Liu
"""

import mysql.connector
from mysql.connector import errorcode

import pandas as pd


class sql:
    """
    Example
    -------
    >>> sql = sql()
    >>> with sql as sql:
    >>>     cid = sql.search("synomys", "synomys","cid","ribitol")[0]
    >>>     name = sql.search("cid","cid","InChI",cid)[0]
    >>>
    """

    def __init__(self, db):
        self.db = db
        self.con = None
        self.cur = None

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
                self.createDB()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cur.close()
        self.con.close()

    def execute_query(self, query, params=None, return_object=True):
        """
        Execute a given SQL query with optional parameters and return the results.
        """
        # print("query-----",query)
        try:
            if params:
                self.cur.execute(query, params)
            else:
                self.cur.execute(query)

            if return_object:
                return self.cur.fetchall()
            else:
                self.con.commit()
        except mysql.connector.Error as err:
            print(f"Error: {err}")
        return None

    def createDB(self):
        self.cur.execute(f"CREATE DATABASE {self.db}")
        print(f"Create {self.db} database...")

    def createTable(self, tableName, arguments) -> None:
        """ Create sql table in database, already contain 'id' and 'created' column
        :param:
        * tableName
        * arguments: exp.'tissue varchar(100), tissueid NOT NULL, varchar(100), PRIMARY KEY (tissue,tissueid)'
        """

        sql = f"""CREATE TABLE IF NOT EXISTS {tableName} 
               (id INT NOT NULL AUTO_increment PRIMARY KEY, 
               created timestamp default current_timestamp,
               {arguments})"""

        self.cur.execute(sql)
        self.con.commit()

    def store(self, tableName, column_dict:dict):
        """

        :param kwargs: tissue,celltype,cellid,donar,age,agetype,genename,expression,study,health_type,cell_num
        :return:
        """
        column, value = list(zip(*column_dict.items()))  # tuple
        placeholder = ",".join(['%s'] * len(column))

        sql = f"""INSERT IGNORE INTO {tableName} {column} VALUES ({placeholder})"""
        print(sql)
        self.cur.execute(sql, value)
        self.con.commit()

    def check_exist(self, table, column, query):
        """
        Parameter:
            * table: "tissue",f"table.cell_type"
            * column: tissue, cell_type, donar, age, age_type, cell_id, gene_name, gene_expression
            * query: query content

        return: boolen
        """

        sql = f"SELECT * FROM {table} WHERE {column} = %s"
        self.cur.execute(sql, (query))
        if self.cur.rowcount == 0:
            return False
        else:
            return True

    def search(self, query_column, table, query):
        """
        return all record, Nonetype was returned when not found
        """

        sql = f"SELECT {query_column} FROM {table} WHERE {query}"
        print("search query--", sql)
        self.cur.execute(sql)
        return self.cur.fetchall()

    # def count(self, table, query_column, query):
    #     sql = f"SELECT {query_column} FROM {table} WHERE {define_column} = %s"

    # def load(self, file_path, name="", set_primary_key: tuple = None):
    #     """Create table by loading txt file
    #
    #     """
    #     delimiter = utilss.detect_delimiter(file_path)
    #
    #     with open(file_path, "r") as t:
    #         file = t.readlines()
    #     columns = file[0].strip("\n").split(delimiter)
    #
    #     # create table
    #     column_input = ','.join([f"`{col}`" for col in columns])
    #     columns_sql = ', '.join([f"`{col}` VARCHAR(255)" for col in columns])
    #
    #     sql1 = f"""CREATE TABLE IF NOT EXISTS {name}
    #     (id INT NOT NULL AUTO_increment PRIMARY KEY,
    #     {columns_sql},
    #     created timestamp default current_timestamp
    #     )"""
    #
    #     print(sql1)
    #     self.cur.execute(sql1)
    #     # input txt content
    #     sql2 = f"LOAD DATA INFILE '{file_path}' INTO TABLE {name} IGNORE 1 LINES ({column_input})"
    #     print(sql2)
    #     self.cur.execute(sql2)
    #     if set_primary_key:
    #         pri_sql = f"ALTER TABLE {name} ADD PRIMARY KEY {set_primary_key}"
    #         self.cur.execute(pri_sql)
    #     self.con.commit()
