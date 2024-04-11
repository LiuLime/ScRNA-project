from sqlalchemy import create_engine, MetaData, select, func
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import sessionmaker
from .model import Base
from .model import Abstract, Log10p, dynamic_corr_class
from utils.log import logger
import threading


class databaseHelper:
    def __init__(self, db_path):
        self.db_path = db_path
        self.engine = create_engine('mysql+mysqlconnector://root@localhost:3306/' + self.db_path)
        self.lock = threading.Lock()
        self.log = logger()
        self.log.debug("init start")

    def connect(func):
        def wrapper(self, *args, **kwargs):
            self.lock.acquire()
            conn = self.engine.connect()
            success, information = func(self, *args, **kwargs)
            conn.close()
            self.lock.release()
            return success, information

        return wrapper

    @connect
    def init(self) -> (bool, str):
        Base.metadata.create_all(self.engine, checkfirst=True)
        return True, "初始化了...而已"

    @connect
    def get_available_tasks(self) -> (bool, list):
        success = True
        available_tasks = []
        try:
            Session = sessionmaker(bind=self.engine)
            session = Session()
            all_tasks = session.query(Abstract).all()
            tasks = [task.to_dict() for task in all_tasks]
            available_tasks = []
            for task in tasks:
                if task["Table"] == 'log10p':
                    pass
                else:
                    available_tasks.append(task)
        except Exception:
            success = False
        session.close()
        return success, available_tasks


    def get_log10p_records_above_threshold(self, threshold) -> (bool, list):
        success = True
        try:
            Session = sessionmaker(bind=self.engine)
            session = Session()
            compliant_table = session.query(Log10p).filter(func.abs(Log10p.log10p) >= threshold).all()

            tasks = [task.to_dict() for task in compliant_table]
            # print("tasks first record--------", tasks[0])
        except Exception as e:
            self.log.error(f"ERROR----->{e}")
            success = False
        session.close()
        return success, tasks


    def get_correlation_records_above_threshold(self, table_name, threshold) -> (bool, list):
        success = True
        try:
            Session = sessionmaker(bind=self.engine)
            session = Session()
            corr_table = dynamic_corr_class(table_name)
            compliant_table = session.query(corr_table).filter(
                func.abs(corr_table.cor_pearson_binMean) >= threshold).all()

            tasks = [task.to_dict() for task in compliant_table]
        except Exception as e:
            self.log.error(f"ERROR----->{e}")
            success = False
        session.close()
        return success, tasks


