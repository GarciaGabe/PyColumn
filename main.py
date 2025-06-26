from funcs import read_file, clear_console
from concrete_section import concrete_section
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')

# COMPRESSION = POSITIVE

clear_console()
print('Starting...')

fck, vertices, reinf_bars = read_file('datpil.txt')
section1 = concrete_section(vertices, reinf_bars, fck)

section1.plot_concrete_section()
section1.plot_verify_diagram(100)

#fck, vertices, reinf_bars = read_file('fig4-6-1-JM.txt')
#ex461JM = concrete_section(vertices, reinf_bars, fck)
#ex461JM.plot_concrete_section()
#ex461JM.plot_verify_diagram(1500)

#fck, vertices, reinf_bars = read_file('section2.txt')
#section2 = concrete_section(vertices, reinf_bars, fck)
#section2.plot_concrete_section()

axial = np.linspace(0, 2400, 10)
diagrams = []

for NSd in axial:
    diagrams.append(section1.verify_section(NSd))

diagrams = np.array(diagrams)


colors = ['red', 'blue', 'black', 'orange', 'green', 'gray', 'purple', 'pink', 'brown', 'yellow']

for diagram, color in zip(diagrams, colors):
    plt.plot(diagram[:, 2], diagram[: , 1], color = color, linewidth = 1.5, label = round(diagram[0, 0],0))

plt.legend()
plt.show()
