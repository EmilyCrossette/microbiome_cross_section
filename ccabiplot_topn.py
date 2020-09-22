import matplotlib.pyplot as plt
import numpy as np

def cca(A, B, qa, qb):
    """
    Canonical Correlation analysis with PC pre-processing.

    Parameters:
    -----------
    A : two dimensional array
        Data for the first set of variables.  The variables are the
        rows and the cases are the columns.
    B : two dimensional array
        Data for the second set of variables.  The variables are the
        rows and the cases are the columns.
    qa : integer
        The number of PCs to retain for the first set of variables.
    qb : integer
        The number of PCs to retain for the second set of variables.

    Note that A and B should have the same number of columns, and the
    columns should be aligned so that the k^th column of A and the
    k^th column of B hold data for the same unit.

    Returns:
    --------
    wa : two-dimensional array
        The coefficients for the first set of variables.
    wb : two-dimensional array
        The coefficients for the second set of variables.
    s : The canonical correlations
    """

    assert(A.shape[1] == B.shape[1])
    n = A.shape[1]

    # Center the data
    Ac = A - A.mean(1)[:, None]
    Bc = B - B.mean(1)[:, None]

    # Use PCA to linearly transform A and B to
    # have identity covariance.
    ua, sa, vta = np.linalg.svd(Ac, 0)
    ub, sb, vtb = np.linalg.svd(Bc, 0)
    ua = ua[:, 0:qa]
    ub = ub[:, 0:qb]
    ua = np.sqrt(n-1) * ua / sa[0:qa]
    ub = np.sqrt(n-1) * ub / sb[0:qb]
    Ax = np.dot(ua.T, Ac)
    Bx = np.dot(ub.T, Bc)

    # Get the canonical coefficients in the transformed coordinates
    C = np.dot(Ax, Bx.T) / n
    u, s, vt = np.linalg.svd(C, 0)

    # Map back to the original coordinates
    wa = np.dot(ua, u)
    wb = np.dot(ub, vt.T)

    return wa, wb, s


def biplot(A, wa, n_arrows, arr_labels=None, obj_labels=None):
    """
    Generate a biplot.

    Parameters
    ----------
    A : 2-dimensional array
        The raw data, rows are objects and columns are variables
    wa : w-dimensional array
        The loadings, rows are variables and columns are components
    n_arrow : boolean array
        If provided, number of arrows to plot
    obj_labels : array of strings
        If provided, object labels to be drawn on the plot.
    arr_labels : array of strings
        If provided, arrow labels to be drawn on plot

    Returns a matplotlib axis object that can be used to manipulate
    the plot.
    """

    # Center the data
    Ac = A - A.mean(1)[:, None]

    # Standardized coordinates of the points
    qp = np.dot(Ac.T, wa)
    qp /= np.sqrt(np.sum(qp**2, 0))

    # Standardized coordinates of the arrows
    wx = wa / np.sqrt(np.sum(wa**2, 0))

    plt.clf()

    # Plot points for the objects
    plt.plot(qp[:, 0], qp[:, 1], 'o', color="orange")

    # Optionally plot object labels
    if obj_labels is not None:
        for i in range(qp.shape[0]):
            plt.text(qp[i, 0], qp[i, 1], obj_labels[i], va="top", ha="center", color="orange")

    # Draw arrows
    # Get top n arrow magnitudes
    arrow_mag = [wa[i, 0]**2 + wa[i, 1]**2 for i in range(len(wa))]
    idx = np.argpartition(arrow_mag, -n_arrows)[-n_arrows:]
    arrows = [True if i in idx else False for i in range(len(wa))]
    
    for i in range(len(arrows)):
        if arrows is not None and not arrows[i]:
            continue
        plt.arrow(0, 0, 0.6*wx[i, 0], 0.6*wx[i, 1], color="purple")

        # Space the arrow label outward from the arrow head
        (x, y) = tuple(wx[i, 0:2])
        a = np.arctan2(y, x)
        r = np.sqrt(x**2 + y**2)
        r *= 0.6
        x = r * np.cos(a)
        y = r * np.sin(a)
        
        if x > 0:
            plt.text(x, y, arr_labels[i], va="center", ha="left", color="purple")
        else:
            plt.text(x, y, arr_labels[i], va="center", ha="right", color="purple")

    return plt.gca()


def test():

    # Sample size
    n = 80

    # Number of features for variable set 1
    pa = 400

    # Number of features for variable set 2
    pb = 200

    # Number of PCs to retain
    qa = 5
    qb = 4

	# This will be the dominant factor for both A and B.
    v1 = np.random.normal(size=n)

    # The coefficients for A are identical.
    A = np.outer(np.ones(pa), v1) + 2*np.random.normal(size=(pa, n))

    # The coefficients for B are alternating.
    c = (-1)**np.mod(np.arange(pb), 2)
    B = np.outer(c, v1) + 2*np.random.normal(size=(pb, n))

    wa, wb, s = cca(A, B, 10, 10)

    # Check the canonical correlation
    aa = np.dot(wa[:, 0], A - A.mean(1)[:, None])
    bb = np.dot(wb[:, 0], B - B.mean(1)[:, None])
    s0 = np.corrcoef(aa, bb)[0, 1]
    assert(np.abs((s[0] - s0) / s[0]) < 0.05)

    # Check the first set of coefficients
    rb = np.corrcoef(wb[:, 0], c)[0, 1]
    assert(np.abs(rb) > 0.8)

    # Check the second set of coefficients
    cv = np.std(wa[:, 0]) / np.abs(np.mean(wa[:,0]))
    assert(np.abs(cv) < 0.5)

    # Check the left biplot
    ax = biplot(A, wa)
    plt.savefig("biplot_left.pdf")

    # Check the right biplot
    arrows = [wb[i, 0]**2 + wb[i, 1]**2 > 0.0002 for i in range(pb)]
    ax = biplot(B, wb, arrows=arrows, obj_labels=[str(i) for i in range(n)])
    plt.savefig("biplot_right.pdf")