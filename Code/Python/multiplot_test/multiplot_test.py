 import matplotlib.pyplot as plt
 import numpy as np




 aa = np.linspace(0, 10)

 a1 = aa*5
 a2 = aa*-5
 a3 = [10 for a in aa]
 a4 = [0 for a in aa]
 a5 = -aa*-5
 a5 = aa*-5
 a6 = aa*5

 fig, axs = plt.subplots(2, 3)
 [[ax1, ax3, ax5], [ax2, ax4, ax6]] = axs
 ax1.plot(aa, a1)
 ax2.plot(aa, a2)
 ax3.plot(-aa, a3)
 ax4.plot(aa, a4)
 ax5.plot(aa, a5)
 ax6.plot(-aa, a6)
