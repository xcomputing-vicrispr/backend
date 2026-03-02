import matplotlib, warnings, sys, json, os
warnings.filterwarnings("ignore")
matplotlib.use('Agg')
import numpy as np

base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
azimuth_path = os.path.abspath(os.path.join(base_dir, 'Azimuth'))
sys.path.insert(0, azimuth_path)
from azimuth import model_comparison as mc

seq_json = sys.argv[1]
seq_list = json.loads(seq_json)
seqs = np.array(seq_list, dtype=object)
preds = mc.predict(seqs)

print(json.dumps(list(preds)))


