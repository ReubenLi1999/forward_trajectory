import numpy as np
import dask.dataframe as dd
data = np.array(
    dd.read_csv('..//output//coordinates_earth_fixed_system.txt', sep='\s+',
    dtype=np.float128)
)
print(data)
