import matplotlib
import matplotlib.pyplot as plt
from funcs import *
import scipy.optimize as opt
matplotlib.use('TkAgg')


class concrete_section:

    def __init__(self, vertices, reinf_bars, fck ):
        
        self.vertices, self.reinf_bars = translate_to_geometric_center(vertices, reinf_bars)
        self.concrete = concrete_characteristics(fck)

    def verify_section(self, NSd):

        diagram = np.zeros((360, 3))

        for angle in range(0,360):

            print("Starting calculating with angle: " + str(angle) + " degrees")

            rotated_vertices, rotated_reinf_bars = rotate_axes(self.vertices, self.reinf_bars, angle)
            y_max, y_min, height = calculate_heights(rotated_vertices)
            effective_heights = calculate_effective_heights(y_max, rotated_reinf_bars)
            

            # searching for x0

            def func(x0):
                tensions = calculate_reinf_tensions(x0, effective_heights, rotated_reinf_bars, self.concrete['peak_compressive_strain'], self.concrete['ultimate_compressive_strain'], height)
                compressed_vertices = find_compression_vertices(rotated_vertices, x0, self.concrete)
                compressed_vertices = rotate_axes(compressed_vertices, rotated_reinf_bars, -angle)[0]
                NRd = calculate_axial_and_moments(compressed_vertices, self.reinf_bars, tensions, self.concrete)[0]
                return NSd - NRd
            
            funcx = lambda x0: func(x0)
            print("Starting optimize function:")
            x0 = opt.root_scalar(funcx, x0 = height/2).root
            print("Found x0: " + str(x0) + "\n")

            tensions = calculate_reinf_tensions(x0, effective_heights, rotated_reinf_bars, self.concrete['peak_compressive_strain'], self.concrete['ultimate_compressive_strain'], height)
            compressed_vertices = find_compression_vertices(rotated_vertices, x0, self.concrete)
            compressed_vertices = rotate_axes(compressed_vertices, rotated_reinf_bars, -angle)[0]
            NRd, Mx, My = calculate_axial_and_moments(compressed_vertices, self.reinf_bars, tensions, self.concrete)
            diagram[angle] = [NRd, Mx, My]

        return diagram
   
    def plot_concrete_section(self):
        plot_section(self.vertices, self.reinf_bars)

    def plot_verify_diagram(self, NSd, x0 = None):
        diagram = self.verify_section(NSd)
        ax = plt.gca() # Get the current axes
        ax.plot(diagram[:, 1], diagram[: , 2], color = 'red', linewidth = 2) # Corrected typo: linewidht -> linewidth

        max_abs_val = np.max(np.abs(diagram[:, 1:3])) # Get max absolute value from Mx and My
        limit = max_abs_val * 1.1 # Add a small margin
        ax.set_xlim([-limit, limit])
        ax.set_ylim([-limit, limit])
        ax.set_aspect('equal', adjustable='box') # Ensure equal aspect ratio and adjust the plot box
        ax.set_xlabel('Resisting Moment X (kN.m)')
        ax.set_ylabel('Resisting Moment Y (kN.m)')
        ax.grid()
        

        plt.show()

    def calculate_resisting_moments(self, NSd = 0, angle = 0, x0 = None):
        print("Starting calculating with angle: " + str(angle) + " degrees")

        rotated_vertices, rotated_reinf_bars = rotate_axes(self.vertices, self.reinf_bars, angle)
        y_max, y_min, height = calculate_heights(rotated_vertices)
        effective_heights = calculate_effective_heights(y_max, rotated_reinf_bars)
            
        # searching for x0

        if (x0 is None):

            def func(x0):
                tensions = calculate_reinf_tensions(x0, effective_heights, rotated_reinf_bars, self.concrete['peak_compressive_strain'], self.concrete['ultimate_compressive_strain'], height)
                compressed_vertices = find_compression_vertices(rotated_vertices, x0, self.concrete)
                compressed_vertices = rotate_axes(compressed_vertices, rotated_reinf_bars, -angle)[0]
                NRd = calculate_axial_and_moments(compressed_vertices, self.reinf_bars, tensions, self.concrete)[0]
                return NSd - NRd
            
            funcx = lambda x0: func(x0)
            print("Starting optimize function:")
            x0 = opt.root_scalar(funcx, x0 = height/2).root
            print("Found x0: " + str(x0) + "\n")

        tensions = calculate_reinf_tensions(x0, effective_heights, rotated_reinf_bars, self.concrete['peak_compressive_strain'], self.concrete['ultimate_compressive_strain'], height)
        compressed_vertices = find_compression_vertices(rotated_vertices, x0, self.concrete)
        compressed_vertices = rotate_axes(compressed_vertices, rotated_reinf_bars, -angle)[0]
        NRd, Mx, My = calculate_axial_and_moments(compressed_vertices, self.reinf_bars, tensions, self.concrete)
        return np.array([NRd, Mx, My])