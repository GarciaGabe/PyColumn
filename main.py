import time
from funcs import read_file, clear_console, clear_console
from concrete_section import concrete_section
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')

# COMPRESSION = POSITIVE

clear_console()
print('Starting...')
start_time = time.time()
As = 0

fck, vertices, reinf_bars = read_file('examples//fig5-4-3-JM.txt')
section1 = concrete_section(vertices, reinf_bars, fck)

#fig1, ax1 = section1.plot_concrete_section()
#ax1.set_title("Concrete Section")
#fig2, ax2 = section1.plot_verify_diagram(1.4*1000)

#section1.verify_section(1.4*4000, angle = 0, iprint = True)

As = section1.design_section(1.4*2000, 1.4*8000, 1.4*8000, iprint = True)


print('\n\n--------------------------------------------------------')
print(f'Final reinforcement area: {As:.2f} cm2')

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time:.2f} seconds")

plt.show()
