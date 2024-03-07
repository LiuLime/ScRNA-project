"""
Save ScRNA result into MySQL database

@ Liu
"""

import pymysql


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

    def __init__(self, db="ScRNA"):
        self.db = db

    def __enter__(self):
        self.connection = pymysql.Connection(host="localhost",
                                             user="root",
                                             db=self.db,
                                             charset='utf8mb4')
        self.cur = self.connection.cursor()
        self.cur.execute(f"USE {self.db}")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cur.close()
        self.connection.close()

    def createTable(self, tableName):
        """ Create tissue_cell_type table in database

        tableName = f"tissue_cell_type"
        """

        sql = f"""CREATE TABLE IF NOT EXISTS {tableName} 
               (id INT NOT NULL AUTO_increment PRIMARY KEY, 
               donar varchar(100), 
               tissue varchar(100) NOT NULL, 
               cell_type varchar(100), 
               cell_id varchar(1000) NOT NULL, 
               gene_name varchar(100) NOT NULL,
               gene_expression DOUBLE,
               created timestamp default current_timestamp,
               PRIMARY KEY (cell_id, gene_name))"""

        self.cur.execute(sql)
        self.connection.commit()

    def loadTable(self, file_path, table_name=""):
        """Create table by loading txt file

        """
        with open(file_path, "r") as t:
            file = t.readlines()
        columns = file[0].strip("\n").split("\t")

        # create table
        column_input = ','.join([f"`{col}`" for col in columns])
        columns_sql = ', '.join([f"`{col}` VARCHAR(255)" for col in columns])
        sql1 = f"""CREATE TABLE IF NOT EXISTS {table_name} 
        (id INT NOT NULL AUTO_increment PRIMARY KEY, 
        {columns_sql}, 
        created timestamp default current_timestamp)"""
        self.cur.execute(sql1)
        # input txt content
        sql2 = f"LOAD DATA INFILE '{file_path}' INTO TABLE {table_name} IGNORE 1 LINES ({column_input})"
        self.cur.execute(sql2)
        self.connection.commit()

    def store(self, tableName=None, infoLevel=False, **kwargs):
        """

        :param kwargs: tissue,celltype,cellid,donar,age,agetype,genename,expression,study,health_type,cell_num
        :return:
        """

        tissue = kwargs.get("tissue", "None")
        celltype = kwargs.get("celltype", "None")
        cellid = kwargs.get("cellid", "None")
        donar = kwargs.get("donar", "None")
        age = kwargs.get("age", "None")
        agetype = kwargs.get("agetype", "None")
        genename = kwargs.get("genename", "None")
        expression = kwargs.get("expression", 0)
        study = kwargs.get("study", "None")
        health_type = kwargs.get("health_type", "None")
        cell_num = kwargs.get("cell_num", -1)
        if infoLevel is True:
            sql = f"""INSERT IGNORE INTO infos 
            (tissue, cell_type, donar, age, age_type, study, health_type, cell_num) 
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s)"""
            self.cur.execute(sql, (tissue, celltype, donar, age, agetype, study, health_type, cell_num))
        else:
            sql1 = f"""INSERT IGNORE INTO {tableName} 
            (donar, tissue, cell_type,cell_id, gene_name, gene_expression) 
            VALUES (%s, %s, %s, %s, %s, %s)"""
            self.cur.execute(sql1, (donar, tissue, celltype, cellid, genename, expression))

        self.connection.commit()

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

    def search(self, table, define_column, query_column, query):
        """
        return one record
        """
        sql = f"SELECT {query_column} FROM {table} WHERE {define_column} = %s"
        self.cur.execute(sql, (query))
        return self.cur.fetchone()

    def count(self, table, query_column, query):
        sql = f"SELECT {query_column} FROM {table} WHERE {define_column} = %s"