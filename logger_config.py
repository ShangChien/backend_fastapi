import logging
from logging import Logger

def get_logger(name) -> Logger:
    logger: Logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # 创建一个文件处理器，用于将日志输出到文件
    file_handler = logging.FileHandler('app.log')
    file_handler.setLevel(logging.DEBUG)

    # 定义日志格式
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)

    # 将处理器添加到日志记录器
    logger.addHandler(file_handler)

    return logger
