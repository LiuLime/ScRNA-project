"""
Save ScRNA result into MySQL database

@ Liu
"""

import mysql.connector
from mysql.connector import errorcode
from utils.log import logger
import math
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

    def get_avalible_studies(self) -> list:
        """return the studies list in database exclude log10p"""
        query = f"SELECT * FROM `abstract`"
        self.cur.execute(query)
        results = self.cur.fetchall()
        tasks = [table for idx, table in results if table != "log10p"]
        return tasks

    # def create_index(self, _tasks: list):
    #     """ create index for multiple columns
    #     :param: db: str, database name
    #     :param: _tasks: list, list of table names
    #     """
    #     count = 0
    #     for task in _tasks:
    #         self.log.debug(f"Start create index |||||| {self.db}.{task}")
    #         query1 = f"ALTER TABLE `{task}` MODIFY COLUMN `gene1` varchar(30)"
    #         query2 = f"ALTER TABLE `{task}` MODIFY COLUMN `gene2` varchar(30)"
    #         query3 = f"ALTER TABLE `{task}` MODIFY COLUMN `cor_pearson.binMean` DOUBLE"
    #         query4 = f"DROP INDEX `idx_corr` ON `{task}` "
    #         query5 = f"ALTER TABLE `{task}` ADD INDEX `idx_corr` (`gene1`,`gene2`,`cor_pearson.binMean`)"
    #         self.cur.execute(query1)
    #         self.cur.execute(query2)
    #         self.cur.execute(query3)
    #         self.cur.execute(query4)
    #         self.cur.execute(query5)
    #
    #         self.con.commit()
    #         count += 1
    #         self.log.debug(f"Finish create index----->{count}/{len(_tasks)} | {task}")
    #
    # def create_log10p_table(self) -> None:
    #     """Load log10p file into mysql table, by create table `log10p` and then insert rows"""
    #     drop_table_stmt = "DROP TABLE `log10p`"
    #     self.cur.execute(drop_table_stmt)
    #     self.con.commit()
    #     self.log.debug(f"(create_log10p_table) drop log10p table 😊")
    #
    #     creat_table_stmt = ("CREATE TABLE IF NOT EXISTS `log10p` "
    #                         "(`id` INT PRIMARY KEY AUTO_INCREMENT,"
    #                         "`gene` varchar(30), "
    #                         "`log10p` DOUBLE, "
    #                         "`abbr_id` varchar(15), "
    #                         "INDEX `idx_log10p` (`gene`, `abbr_id`))")
    #     self.cur.execute(creat_table_stmt)
    #     self.con.commit()
    #     self.log.debug(f"(create_log10p_table) create log10p table 🤝")
    #
    # def df_to_sql(self, df: pd.DataFrame, table_name: str) -> None:
    #     """insert dataframe rows into mysql table"""
    #     df_dict = (df.dropna()
    #                .replace(float('inf'), 999)
    #                .replace(float('-inf'), -999))
    #     # self.log.debug(f"length of dataframe after dropna{len(df_dict)}")
    #
    #     df_dict = df_dict.to_dict(orient="list")
    #     df_list_record = list(zip(*df_dict.values()))  # mysql table record content, each row as tuple
    #     self.log.debug(f"(df_to_sql) {len(df_list_record)} NOT NA records waiting import to mysql····")
    #
    #     columns = f"{tuple(df_dict.keys())}".replace('\'', "`")  # mysql table column name
    #     trunksize = 5000
    #     truncktime = math.ceil(len(df_list_record) / trunksize)
    #
    #     for times in range(truncktime):
    #         batch = df_list_record[times * trunksize: (times + 1) * trunksize]
    #         value = ','.join(str(i) for i in batch)
    #         # stmt = f"INSERT INTO `{table_name}` {columns} VALUES {value}"
    #         stmt = "INSERT INTO `{}` {} VALUES {}".format(table_name, columns, value)  # batch-> list[tuples]
    #         self.cur.execute(stmt)
    #     self.con.commit()
    #     self.log.debug(f"(df_to_sql) finish load file ✌️✌️✌️")
    #
    # def join_corr_log10p(self, tasks: list):
    #
    #     # query = f"""
    #     #         CREATE TABLE new_table AS
    #     #         SELECT t1.*, lg1.log10p AS gene1_log10p, lg2.log10p AS gene2_log10p
    #     #         FROM `{task}` AS t1
    #     #         JOIN `log10p` AS lg1 ON t1.gene1 = lg1.gene AND lg1.abbr_id = '{task}'
    #     #         JOIN `log10p` AS lg2 ON t1.gene2 = lg2.gene AND lg2.abbr_id = '{task}'
    #     #         WHERE t1.`cor_pearson.binMean` > {corr_threshold} AND lg1.log10p > {log10p_threshold} AND lg2.log10p > {log10p_threshold};
    #     #         """
    #
    #     for task in tasks:
    #         table_name = f"joined_{task}"
    #         query1 = f"""CREATE TABLE IF NOT EXISTS `{table_name}`
    #                             (`id` INT PRIMARY KEY AUTO_INCREMENT,
    #                             `gene1` varchar(30),
    #                             `gene2` varchar(30),
    #                              `cor_pearson.binMean` DOUBLE,
    #                             `log10p1` DOUBLE,
    #                             `log10p2` DOUBLE
    #
    #                             )
    #
    #         """
    #         self.cur.execute(query1)
    #         self.con.commit()
    #         self.log.debug(f"(join_by_corr_log10p) {self.db} create `{table_name}` table 🤝")
    #
    #         query2 = f"""INSERT INTO `{table_name}` (`gene1`,`gene2`,`cor_pearson.binMean`,`log10p1`,`log10p2`)
    #                    SELECT  a.gene1,
    #                            a.gene2,
    #                            a.`cor_pearson.binMean`,
    #                            a.`log10p`      AS p1,
    #                            log10p.`log10p` AS p2
    #
    #                     FROM (SELECT {task}.`gene1`,
    #                                  {task}.`gene2`,
    #                                  {task}.`cor_pearson.binMean`,
    #                                  log10p.`log10p`
    #
    #                   FROM {task}
    #                            INNER JOIN log10p ON ({task}.`gene1` = log10p.`gene` and log10p.`abbr_id` = '{task}')) AS a
    #                      INNER JOIN log10p ON (a.`gene2` = log10p.`gene` and log10p.`abbr_id` = '{task}')
    #                 ;
    #                 """
    #
    #         self.cur.execute(query2)
    #         self.con.commit()
    #         self.log.debug(f"(join_by_corr_log10p) {self.db} `{table_name}` finish ✌️✌️✌️")
    #     # return self.cur.fetchall()
    #
    # def createTable(self, tableName, arguments) -> None:
    #     """ Create sql table in database, already contain 'id' and 'created' column
    #     :param:
    #     * tableName
    #     * arguments: exp.'tissue varchar(100), tissueid NOT NULL, varchar(100), PRIMARY KEY (tissue,tissueid)'
    #     """
    #
    #     sql = f"""CREATE TABLE IF NOT EXISTS {tableName}
    #            (id INT NOT NULL AUTO_increment PRIMARY KEY,
    #            created timestamp default current_timestamp,
    #            {arguments})"""
    #
    #     self.cur.execute(sql)
    #     self.con.commit()
    #
    # def store(self, tableName, column_dict: dict):
    #     """
    #
    #     :param kwargs: tissue,celltype,cellid,donar,age,agetype,genename,expression,study,health_type,cell_num
    #     :return:
    #     """
    #     column, value = list(zip(*column_dict.items()))  # tuple
    #     placeholder = ",".join(['%s'] * len(column))
    #
    #     sql = f"""INSERT IGNORE INTO {tableName} {column} VALUES ({placeholder})"""
    #     print(sql)
    #     self.cur.execute(sql, value)
    #     self.con.commit()
    #
    # def check_exist(self, table, column, query):
    #     """
    #     Parameter:
    #         * table: "tissue",f"table.cell_type"
    #         * column: tissue, cell_type, donar, age, age_type, cell_id, gene_name, gene_expression
    #         * query: query content
    #
    #     return: boolen
    #     """
    #
    #     sql = f"SELECT * FROM {table} WHERE {column} = %s"
    #     self.cur.execute(sql, (query))
    #     if self.cur.rowcount == 0:
    #         return False
    #     else:
    #         return True
    #
    # def search(self, query_column, table, query):
    #     """
    #     return all record, Nonetype was returned when not found
    #     """
    #
    #     sql = f"SELECT {query_column} FROM {table} WHERE {query}"
    #     print("search query--", sql)
    #     self.cur.execute(sql)
    #     return self.cur.fetchall()
