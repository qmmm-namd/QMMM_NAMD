#!/usr/bin/env python3 
from time import time

def timer(function):
    def cal_time(*args, **kwargs):
        start = time()
        function(*args, **kwargs)
        return time() - start
    return cal_time