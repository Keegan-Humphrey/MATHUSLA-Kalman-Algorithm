import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import joblib
import pandas as pd
    
flows = joblib.load('flows.pkl')
cutflow = pd.DataFrame(flows)
print(cutflow)