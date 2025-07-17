import time
from funcs import clear_console
from column import RectangularColumn
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')

# COMPRESSION = POSITIVE

clear_console()
print('Starting...')
start_time = time.time()
As = 0

test_column = RectangularColumn(40, 19, 30, 300)
test_column.set_forces([100, 100], [4000, 4000], [1500, -1500])
test_column.calculate_second_order_effects(method='k_approx')
As = test_column.designColumn()

print('\n\n--------------------------------------------------------')
print(f'Final reinforcement area: {As:.2f} cm2')


end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time:.2f} seconds")
print('--------------------------------------------------------\n\n\n')
