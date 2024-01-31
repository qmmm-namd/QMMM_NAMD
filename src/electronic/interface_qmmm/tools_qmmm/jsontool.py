#!/usr/bin/env python3 
import json
from os import path

def load_json(json_file):
    if not path.exists(json_file):
        print('  File (%s) not exists!'%json_file)
        exit()

    with open(json_file) as jf:
        return json.load(jf)


def dump_json(json_file, obj, format=None):
    json.encoder.FLOAT_REPR = lambda f: format("%.18g" % f)
    with open(json_file, 'w') as jf:
        jf.write(json.dumps(obj, indent=2))

    return 
    

if __name__ == '__main__':
    # test 
    print(load_json('../additive/config.json'))
