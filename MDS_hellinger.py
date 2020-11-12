import pandas as pd
import numpy as np
from sklearn import manifold
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


def perform_mds(di):

    # Create an incidence matrix where the rows are samples and columns
    # are OTU's.
    di = di.fillna(0)

    # Calculate the total per sample and normalize.
    tot = di.sum(1)
    di = di.divide(tot, axis=0)

    # Sort the columns in descending abundance order
    ii = np.argsort(di.sum(0))
    di = di.iloc[:, ii[::-1]]

    # Shannon diversity per sample
    # shdiver = -(di * np.log(di)).sum(1)

    # Simpson diversity per sample
    # sidiver = 1 - (di**2).sum(1)

    # Inverse Simpson diversity per sample
    # isdiver = 1 / (di**2).sum(1)

    # Pairwise Hellinger distances
    n = di.shape[0]
    hdist = np.zeros((n, n))
    dq = np.sqrt(np.asarray(di))
    for i in range(n):
        for j in range(i):
            d = np.sum((dq[i, :] - dq[j, :])**2)
            d = np.sqrt(d) / np.sqrt(2)
            hdist[i, j] = d
            hdist[j, i] = d
    na = di.index.tolist()
    hdist = pd.DataFrame(hdist, index=na, columns=na)

    # Use multidimensional scaling to visualize the data
    mds = manifold.MDS(n_components=2, dissimilarity='precomputed')
    pos = mds.fit(hdist).embedding_

    return pos, mds