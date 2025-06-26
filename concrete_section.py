import matplotlib
import matplotlib.pyplot as plt
from funcs import *
import scipy.optimize as opt
matplotlib.use('TkAgg')


class concrete_section:

    def __init__(self, vertices, reinf_bars, fck ):
        
        self.vertices, self.reinf_bars = translate_to_geometric_center(vertices, reinf_bars)
        self.concrete = concrete_characteristics(fck)

    def find_x0(self, effective_heights, rotated_reinf_bars, rotated_vertices, angle, NSd, height):
        
        def func(x0):
            tensions = calculate_reinf_tensions(x0, effective_heights, rotated_reinf_bars, self.concrete['peak_compressive_strain'], self.concrete['ultimate_compressive_strain'], height)
            compressed_vertices = find_compression_vertices(rotated_vertices, x0, self.concrete)
            compressed_vertices = rotate_axes(compressed_vertices, rotated_reinf_bars, -angle)[0]
            NRd = calculate_axial_and_moments(compressed_vertices, self.reinf_bars, tensions, self.concrete)[0]
            return NSd - NRd
        
        print("Starting optimize function:")
        
        #x0 = opt.root_scalar(func, x0 = height/2).root
        limits = [0 , height]
        while (func(limits[0]) * func(limits[1]) > 0):
            print('bissec limits increased: ', limits[1])
            limits[1] *= 1000

        x0 = (limits[0] + limits[1]) / 2

        while abs(func(x0)) > 1e-4:
            #print(limits)
            results = [func(limits[0]), func(limits[1])]
            res_bissec = func(x0)

            if (results[0] * res_bissec < 0):
                limits[1] = x0
            elif (results[1] * res_bissec < 0):
                limits[0] = x0
            x0 = (limits[0] + limits[1]) / 2
        
        print("Found x0: " + str(x0) + "\n")    
        return x0

    def verify_section(self, NSd):

        diagram = np.zeros((360, 3))

        for angle in range(0,360):

            print("Starting calculating with angle: " + str(angle) + " degrees")

            rotated_vertices, rotated_reinf_bars = rotate_axes(self.vertices, self.reinf_bars, angle)
            y_max, y_min, height = calculate_heights(rotated_vertices)
            effective_heights = calculate_effective_heights(y_max, rotated_reinf_bars)
            
            x0 = self.find_x0(effective_heights, rotated_reinf_bars, rotated_vertices, angle, NSd, height)

            tensions = calculate_reinf_tensions(x0, effective_heights, rotated_reinf_bars, self.concrete['peak_compressive_strain'], self.concrete['ultimate_compressive_strain'], height)
            compressed_vertices = find_compression_vertices(rotated_vertices, x0, self.concrete)
            compressed_vertices = rotate_axes(compressed_vertices, rotated_reinf_bars, -angle)[0]
            NRd, Mx, My = calculate_axial_and_moments(compressed_vertices, self.reinf_bars, tensions, self.concrete)
            diagram[angle] = [NRd, Mx, My]
            print(f'NRd: {NRd:.2f} kN, Mx: {Mx/100:.2f} kN.m, My: {My/100:.2f} kN.m')
        return diagram
   
    def plot_concrete_section(self):
        plot_section(self.vertices, self.reinf_bars)

    def plot_verify_diagram(self, NSd, x0 = None):
        diagram = self.verify_section(NSd) / 100
        fig, ax = plt.subplots(figsize=(8, 8)) # Create a new square figure and axes
        ax.plot(diagram[:, 2], diagram[: , 1], color = 'red', linewidth = 2)

        ax.spines['right'].set_color('black')
        ax.spines['top'].set_color('black')

        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.axvline(x=0, color='k', linewidth=1) # Plot Y-axis at X=0
        ax.axhline(y=0, color='k', linewidth=1) # Plot X-axis at Y=0
        
        max_abs_val = np.max(np.abs(diagram[:, 1:3])) # Get max absolute value from Mx and My (columns 1 and 2)
        # Ensure the limit is at least a small positive value to handle cases where diagram values are all zero
        limit = max(max_abs_val * 1.1, 1.0) # Add a small margin, with a minimum of 1.0
        ax.set_xlim([-limit, limit])
        ax.set_ylim([-limit, limit])
        ax.set_aspect('equal', adjustable='box') # Ensure equal aspect ratio and adjust the plot box
        ax.set_xlabel('MRd,y (kN.m)')
        ax.set_ylabel('MRd,x (kN.m)')
        ax.grid(alpha = 0.4)
        
        plt.show()
    
    def calculate_resisting_moments(self, NSd = 0, angle = 0, x0 = None, plot = False):
        print("Starting calculating with angle: " + str(angle) + " degrees")

        rotated_vertices, rotated_reinf_bars = rotate_axes(self.vertices, self.reinf_bars, angle)
        y_max, y_min, height = calculate_heights(rotated_vertices)
        effective_heights = calculate_effective_heights(y_max, rotated_reinf_bars)
            
        # searching for x0

        if (x0 is None):
            x0 = self.find_x0(effective_heights, rotated_reinf_bars, rotated_vertices, angle, NSd, height)

        tensions = calculate_reinf_tensions(x0, effective_heights, rotated_reinf_bars, self.concrete['peak_compressive_strain'], self.concrete['ultimate_compressive_strain'], height)
        compressed_vertices = find_compression_vertices(rotated_vertices, x0, self.concrete)
        compressed_vertices = rotate_axes(compressed_vertices, rotated_reinf_bars, -angle)[0]
        NRd, Mx, My = calculate_axial_and_moments(compressed_vertices, self.reinf_bars, tensions, self.concrete)

        if plot:
            plot_section(self.vertices, self.reinf_bars, compressed_vertices, tensions[1])

        return np.array([NRd, Mx, My])
    

