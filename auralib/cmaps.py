"""
AuraQI module providing various custom colormaps.

Author:   Wes Hamlyn
Created:   19-Dec-2016
Last Mod:  19-Dec-2016

"""

import numpy as np
from matplotlib.colors import ListedColormap


def petrel(direction='normal'):
    """
    Returns the equivalent of the colormap named "Petrel" in RokDoc.
    """
    
    import numpy as np
    import matplotlib as mpl
    
    n1 = 80
    r1 = np.linspace(255, 191, n1)
    g1 = np.linspace(255, 0, n1)
    b1 = np.linspace(0, 0, n1)
    
    n2 = 20
    r2 = np.linspace(191, 97, n2)
    g2 = np.linspace(0, 69, n2)
    b2 = np.linspace(0, 0, n2)
    
    n3 = 20
    r3 = np.linspace(97, 202, n3)
    g3 = np.linspace(69, 202, n3)
    b3 = np.linspace(0, 202, n3)
    
    r4 = np.linspace(202, 77, n3)
    g4 = np.linspace(202, 77, n3)
    b4 = np.linspace(202, 77, n3)
    
    r5 = np.linspace(77, 0, n2)
    g5 = np.linspace(77, 0, n2)
    b5 = np.linspace(77, 191, n2)
    
    r6 = np.linspace(0, 161, n1)
    g6 = np.linspace(0, 255, n1)
    b6 = np.linspace(191, 255, n1)
    
    r = np.hstack([r1, r2, r3, r4, r5, r6])
    g = np.hstack([g1, g2, g3, g4, g5, g6])
    b = np.hstack([b1, b2, b3, b4, b5, b6])
    cmap = np.vstack([r, g, b]).T / 255.0
    
    if direction == 'reverse':
        cmap = np.flipud(cmap)
        
    return ListedColormap(cmap, name='petrel')
    