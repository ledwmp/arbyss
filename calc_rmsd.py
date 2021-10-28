from scipy.spatial.transform import Rotation as R
import numpy as np

def align_chunks(chunk_a,chunk_b):
    array_a = np.array([[i._coord._X,i._coord._Y,i._coord._Z] for i in chunk_a])
    array_b = np.array([[i._coord._X,i._coord._Y,i._coord._Z] for i in chunk_b])

    blah = R.align_vectors(array_a,array_b)
    print(blah[1])
