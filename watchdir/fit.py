from unimol_tools import MolTrain, MolPredict
import numpy as np
from pathlib import Path as P
import pickle
import argparse

parser = argparse.ArgumentParser(description='unimol fine-tuning')
parser.add_argument('--target', '-T', help='拟合的目标值',)
args = parser.parse_args()


list_P = list(P.cwd().glob('*.pkl'))
if len(list_P) == 0:
    raise Exception('No *.pkl found')
elif len(list_P) > 1:
    raise Exception('Too many *.pkl found')

with open(list_P[0],'rb') as f:
    raw = pickle.load(f)
    
data={
    'target':raw[args.target],
    'atoms':raw['atoms'],
    'coordinates':raw['coordinates'],
}

clf = MolTrain(task='regression', data_type='molecule', epochs=10, batch_size=16, metrics='r2')
pred = clf.fit(data = data)

# clf = MolPredict(load_model='../exp')
# res = clf.predict(data = data)
