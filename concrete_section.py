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
        
        #x0 = opt.root_scalar(func, x0 = height/2).root
        limits = [0 , height]
        while (func(limits[0]) * func(limits[1]) > 0):
            limits[1] *= 1000
            if limits[1] > 10e300:
                raise Exception('Overflow in scalar limit')

        x0 = (limits[0] + limits[1]) / 2

        while abs(func(x0)) > 1e-2:
            results = [func(limits[0]), func(limits[1])]
            res_bissec = func(x0)

            if (results[0] * res_bissec < 0):
                limits[1] = x0
            elif (results[1] * res_bissec < 0):
                limits[0] = x0
            x0 = (limits[0] + limits[1]) / 2
        
        return x0

    def verify_section(self, NSd, angle = None , iprint = False):

        if angle is not None:
            rotated_vertices, rotated_reinf_bars = rotate_axes(self.vertices, self.reinf_bars, angle)
            y_max, y_min, height = calculate_heights(rotated_vertices)
            effective_heights = calculate_effective_heights(y_max, rotated_reinf_bars)
            
            x0 = self.find_x0(effective_heights, rotated_reinf_bars, rotated_vertices, angle, NSd, height)

            tensions = calculate_reinf_tensions(x0, effective_heights, rotated_reinf_bars, self.concrete['peak_compressive_strain'], self.concrete['ultimate_compressive_strain'], height)
            compressed_vertices = find_compression_vertices(rotated_vertices, x0, self.concrete)
            compressed_vertices = rotate_axes(compressed_vertices, rotated_reinf_bars, -angle)[0]
            NRd, Mx, My = calculate_axial_and_moments(compressed_vertices, self.reinf_bars, tensions, self.concrete)
            if iprint: print(f'Angle: {angle}, NRd: {NRd:.2f}, Mx: {Mx:.2f}, My: {My:.2f}')
            return [NRd, Mx, My]
        else:
            diagram = np.zeros((360, 3))

            for angle in range(0, 360):

                rotated_vertices, rotated_reinf_bars = rotate_axes(self.vertices, self.reinf_bars, angle)
                y_max, y_min, height = calculate_heights(rotated_vertices)
                effective_heights = calculate_effective_heights(y_max, rotated_reinf_bars)
            
                x0 = self.find_x0(effective_heights, rotated_reinf_bars, rotated_vertices, angle, NSd, height)

                tensions = calculate_reinf_tensions(x0, effective_heights, rotated_reinf_bars, self.concrete['peak_compressive_strain'], self.concrete['ultimate_compressive_strain'], height)
                compressed_vertices = find_compression_vertices(rotated_vertices, x0, self.concrete)
                compressed_vertices = rotate_axes(compressed_vertices, rotated_reinf_bars, -angle)[0]
                NRd, Mx, My = calculate_axial_and_moments(compressed_vertices, self.reinf_bars, tensions, self.concrete)
                diagram[angle] = [NRd, Mx, My]
                if iprint: print(f'Angle: {angle}, NRd: {NRd:.2f}, Mx: {Mx:.2f}, My: {My:.2f}')
            return diagram
   
    def plot_concrete_section(self):
        return plot_section(self.vertices, self.reinf_bars)

    def plot_verify_diagram(self, NSd, x0 = None):
        diagram = self.verify_section(NSd) / 100.0
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
        
        return fig, ax
    
    def calculate_resisting_moments(self, NSd = 0, angle = 0, x0 = None, plot = False):
 
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
    
    def set_steel_area(self, total_rebar_area):
        steel_proportion = self.reinf_bars[:, 2] / np.sum(self.reinf_bars[:, 2])
        self.reinf_bars[:, 2] = steel_proportion * total_rebar_area
    
    def design_section(self, NSd, MSxd, MSyd, iprint = False):

        if (MSxd == 0 and MSyd == 0):
            A_c = calculate_polygon_properties(self.vertices)[0]
            As = (1.01*(NSd -  A_c * self.concrete['fc']) / steel_stress_strain_curve(self.concrete['peak_compressive_strain'], STEEL_YOUNG_MODULUS, STEEL_YIELD_STRESS))
            return np.max([As, 0.00])


        if iprint: print('Starting Design Process...')

        theta_d = np.degrees( np.arctan2( MSyd,  MSxd ))
        final_theta_r = 1000    
        iter = 0
        
        while( (abs(theta_d - final_theta_r) > 0.001 if abs(theta_d) < 1e-6 else abs((theta_d - final_theta_r)/theta_d ) > 0.001) and iter < 11):
            iter += 1
            if iprint: print('------------------------------\n\nDesingning section, iter: {0} - Current Theta Difference = {1:.2f}%'.format(iter, abs((theta_d - final_theta_r)/theta_d )*100 if theta_d != 0 else 9999 ) )
            
            flag = True
            while flag:
                try:
                    diagrama = self.verify_section(NSd)
                    flag = False
                except Exception as e:
                    print(f'Exception raised = {e}')
                    as_atual = np.sum( self.reinf_bars[:,2] )
                    print(f'Setting steel area to = {as_atual * 10:.2f} cm2')
                    self.set_steel_area(as_atual * 10)
                finally:
                    continue

            
            theta_r = np.degrees(np.arctan2(diagrama[:, 2], diagrama[:, 1]))
            
            closest_alpha = np.argmin(np.abs(theta_r - theta_d))
            closest_alpha_plus = closest_alpha + 1
            closest_alpha_minus = closest_alpha - 1
        
            def func(alpha):
                N, Mx, My = self.calculate_resisting_moments(NSd, angle = alpha, plot = False)
                theta_r = np.degrees( np.arctan2( My , Mx ))
                return theta_d - theta_r

            if( func(closest_alpha) * func(closest_alpha_plus) > 0):
                limits = [closest_alpha_minus, closest_alpha]
            else:
                limits = [closest_alpha, closest_alpha_plus]

            if(limits[0] > limits[1]):
                aux = limits[0]
                limits[0] = limits[1]
                limits[1] = aux
              

            if (func(closest_alpha) <= 1e-2):
                alpha_0 = closest_alpha
            else:
                alpha_0 = (limits[0] + limits[1]) / 2
        
            while abs(func(alpha_0)) > 1e-2:
                results = [func(limits[0]), func(limits[1])]
                res_bissec = func(alpha_0)

                if (results[0] * res_bissec < 0):
                    limits[1] = alpha_0
                elif (results[1] * res_bissec < 0):
                    limits[0] = alpha_0
                alpha_0 = (limits[0] + limits[1]) / 2
        
            if iprint: print(f'Found neutral line angle = {alpha_0:.2f} degrees')
            if iprint: print(f'Theta_d = {theta_d:.5f} degrees')
            N, Mx, My = self.calculate_resisting_moments(NSd, angle = alpha_0, plot = False)
            if iprint: print(f'Theta_r = {np.degrees( np.atan( My / Mx )):.5f} degrees')

            A_c = calculate_polygon_properties(self.vertices)[0]
            A_s0 = 1.01*(NSd -  A_c * self.concrete['fc']) / steel_stress_strain_curve(self.concrete['peak_compressive_strain'], STEEL_YOUNG_MODULUS, STEEL_YIELD_STRESS)
        
        
            if A_s0 < 0:
                A_s0 = 0.01


            if iprint: print(f'A_s0 = {A_s0:.2f} cm2')
        
            self.set_steel_area(A_s0)
            NRd, MRxd, MRyd = self.calculate_resisting_moments(NSd, angle = alpha_0, plot = False)
            MRd0 = (MRxd**2 + MRyd**2)**0.5
            MSd0 = (MSxd**2 + MSyd**2)**0.5

            if iprint: print('NRd, MRxd, MRyd for A_s0 = {0:.2f}, {1:.2f}, {2:.2f}'.format(NRd, MRxd, MRyd))

            psi0 = (MRd0 - MSd0) / MSd0

            if psi0 > 0: return 0

            if iprint: print('psi0 = {:.2f}'.format(psi0))

            psiu = -1

            A_su = 0.02 * A_c

            while(psiu < 0):
                A_su *= 2
                self.set_steel_area(A_su)
                NRd, MRxd, MRyd = self.calculate_resisting_moments(NSd, angle = alpha_0, plot = False)
                MRdu = (MRxd**2 + MRyd**2)**0.5
                psiu = (MRdu - MSd0) / MSd0

            area_limits = [A_s0, A_su]
            if iprint: print('psiu = {:.2f}'.format(psiu))
            if iprint: print('Area limits: {0:.2f} cm2, {1:.2f} cm2'.format(area_limits[0], area_limits[1]))
            
            def func(As):
                self.set_steel_area(As)
                NRd, MRxd, MRyd = self.calculate_resisting_moments(NSd, angle = alpha_0, plot = False)
                MRd = (MRxd**2 + MRyd**2)**0.5
                psi = (MRd - MSd0) / MSd0
                return psi
        
            As = (area_limits[0] + area_limits[1]) / 2

            while abs(func(As)) > 1e-4:
                results = [func(area_limits[0]), func(area_limits[1])]
                res_bissec = func(As)

                if (results[0] * res_bissec < 0):
                    area_limits[1] = As
                elif (results[1] * res_bissec < 0):
                    area_limits[0] = As

                As = (area_limits[0] + area_limits[1]) / 2
        
            As *= 1.01
            self.set_steel_area(As)
            if iprint: print(f'Found As = {As:.2f} cm2')
            NRd, MRxd, MRyd = self.calculate_resisting_moments(NSd, angle = alpha_0, plot = False)
            if iprint: print(f'NRd = {NRd:.2f} kN, MRx = {MRxd/100:.2f} kN.m, MRy = {MRyd/100:.2f} kN.m')
            final_theta_r = np.degrees( np.arctan2( MRyd , MRxd ))
            if iprint: print(f'Final theta_r = {final_theta_r:.2f} degrees')
        return As
        
