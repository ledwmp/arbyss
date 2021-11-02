from scipy.spatial.transform import Rotation as R
import numpy as np

def rmsd(points_a,points_b):
    d = np.linalg.norm(points_a - points_b,axis=1)
    rmsd = np.sqrt((1./len(d))*np.sum(d**2))
    return rmsd

def align_chunks(chunk_a,chunk_b):
    """Method to superimpose large chunks, then determine the RMSD of all iterations of small chunks
    """
    array_a = np.array([[i._coord._X,i._coord._Y,i._coord._Z] for i in chunk_a[0]])
    array_b = np.array([[i._coord._X,i._coord._Y,i._coord._Z] for i in chunk_b[0]])

    centroid_a = np.average(array_a,axis=0)
    centroid_b = np.average(array_b,axis=0)

    array_a -= centroid_a
    array_b -= centroid_b

    h = array_a.T @ array_b #calculate covariance matrix
    u,s,vt = np.linalg.svd(h) #calculate SVD of covariance matrix
    v = vt.T

    d = np.linalg.det(v @ u.T)
    e = np.array([[1,0,0],[0,1,0],[0,0,d]])

    r = v @ e @ u.T
    tt = centroid_b - (r @ centroid_a)
    tmp_list = []
    for i in array_a:
        point = (r @ i)
        #point = (r @ i) + tt
        tmp_list.append(np.reshape(point,(1,3)))
    array_a_remap = np.vstack(tmp_list)
    tmp_score = {}
    for k in range(1,len(chunk_a[1])-1):
        for l in range(1,len(chunk_b[1])-1):
            #fix thing here, need to evaulate at 0 and -1 too
            align_a = np.array([[i._coord._X,i._coord._Y,i._coord._Z] for i in chunk_a[1][k-1:k+1]])
            align_b = np.array([[i._coord._X,i._coord._Y,i._coord._Z] for i in chunk_b[1][l-1:l+1]])

            align_a -= centroid_a
            align_b -= centroid_b

            tmp_list = []
            for i in align_a:
                point = (r @ i)
                #point = (r @ i) + tt
                tmp_list.append(np.reshape(point,(1,3)))
            align_a_remap = np.vstack(tmp_list)
            tmp_score[(chunk_a[1][k]._resraw,chunk_b[1][l]._resraw)] = rmsd(align_b,align_a_remap)
    return tmp_score
