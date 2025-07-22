import numpy as np
from consts import *
from os import system
from concrete_section import ConcreteSection
import matplotlib.pyplot as plt


def alpha_b(M1, M2):
    """
    Calculate the alpha_b factor based on the moments M1 and M2.
    :param M1: Moment at one end of the column.
    :param M2: Moment at the other end of the column.
    :return: Calculated alpha_b value.
    """
    if M1 == 0 and M2 == 0:
        return 0.40
    aux1 = np.sign( M1 * M2 )
    aux2 = np.min(np.abs( [M1 , M2] )) / np.max(np.abs( [M1 , M2] ))
    alpha_b = 0.60 + 0.40 * aux1 * aux2
    return np.max([alpha_b, 0.40])

class RectangularColumn:

    def __init__(self, width, height, fck, l_e ,cover = 2.5, stirrups = 0.63, consider_eta_c = True):
        """
        width: Lenght of the column in direction X [cm].
        height: Lenght of the column in direction Y [cm].
        """
        self.width = width
        self.height = height
        self.fck = fck
        self.l_e = l_e
        self.cover = cover
        self.stirrups = stirrups
        self.consider_eta_c = consider_eta_c

    
    def set_forces(self, N, Mx, My):
        self.N = {'bottom': N[0], 'top': N[1]}
        self.Mx = {'bottom': Mx[0], 'top': Mx[1]}
        self.My = {'bottom': My[0], 'top': My[1]}

    def calculate_second_order_effects(self, method='k_approx'):
        """
        Calculate second-order effects based on the method specified.
        :param method: Method to use for second-order effects calculation.Options are 'k_approx' or 'exact'.
        
        """
        if self.N == None or self.Mx is None or self.My is None:
            raise ValueError("Forces must be set before calculating second-order effects.")

        i_x = self.width * self.height**3 / 12
        i_y = self.height * self.width**3 / 12
        area = self.width * self.height
        lambda_x = self.l_e / ( np.sqrt(i_x / area) )
        lambda_y = self.l_e / ( np.sqrt(i_y / area) )

        e1_x = self.Mx['top'] / (self.N['top'] )
        e1_y = self.My['top'] / (self.N['top'] )

        alpha_b_x = alpha_b(self.Mx['bottom'], self.Mx['top'])
        alpha_b_y = alpha_b(self.My['bottom'], self.My['top'])        
        
        lambda_1_x = np.min( [np.max( [( 25 + 12.5 * e1_x/ self.height ) / alpha_b_x, 35]  ), 90] )
        lambda_1_y = np.min( [np.max( [( 25 + 12.5 * e1_y/ self.height ) / alpha_b_y, 35]  ), 90] )

        if method == 'k_approx':
            print("Calculating second-order effects using k_approx method...\n")
            print(f"lambda_x: {lambda_x:.1f}, lambda_y: {lambda_y:.1f}")
            print(f"lambda_1_x: {lambda_1_x:.1f}, lambda_1_y: {lambda_1_y:.1f}")
            print(f"alpha_b_x: {alpha_b_x:.2f}, alpha_b_y: {alpha_b_y:.2f}\n")
            self.N['mid'] = float(self.N['top'] + self.N['bottom'])/2
            
            if lambda_x <= lambda_1_x:
                self.Mx['mid'] = (self.Mx['top'] + self.Mx['bottom']) / 2
            else:
                if lambda_x > 91:
                    raise ValueError("Column is slender in X direction, second-order effects cannot be calculated by k_approx method.")
                else:
                    Nd = np.max( [self.N['top'], self.N['bottom']] ) * 1.0
                    MSd = np.max( np.abs( [self.Mx['top'], self.Mx['bottom']] ) ) * 1.0
                    print(f"Nd: {Nd:.2f} kN, MSd: {MSd:.2f} kNcm")
                    
                    a = 5 * self.height
                    b = (self.height ** 2 * Nd) - (Nd * self.l_e ** 2 / 320) - (5 * self.height * alpha_b_x * MSd)
                    c = - Nd * self.height ** 2 * alpha_b_x * MSd

                    MSd_tot = (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
                    self.Mx['mid'] = float(MSd_tot)
            print(f"Calculated Mx at midspan: {self.Mx['mid']:.2f} kNcm")

            if lambda_y <= lambda_1_y:
                self.My['mid'] = (self.My['top'] + self.My['bottom']) / 2
            else:
                if lambda_y > 91:
                    raise ValueError("Column is slender in Y direction, second-order effects cannot be calculated by k_approx method.")
                else:
                    Nd = np.max( [self.N['top'], self.N['bottom']] ) * 1.0
                    MSd = np.max( np.abs( [self.My['top'], self.My['bottom']] ) ) * 1.0
                    print(f"Nd: {Nd:.2f} kN, MSd: {MSd:.2f} kNcm")
                    
                    a = 5 * self.width
                    b = (self.width ** 2 * Nd) - (Nd * self.l_e ** 2 / 320) - (5 * self.width * alpha_b_x * MSd)
                    c = - Nd * self.width ** 2 * alpha_b_y * MSd

                    MSd_tot = (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
                    self.My['mid'] = float(MSd_tot)
            print(f"Calculated My at midspan: {self.My['mid']:.2f} kNcm")

        #implementar momento mínimo


        elif method == '1/r approx':
            pass


        elif method == 'N-M interaction':    
            pass
        
        else: # method == 'exact'
            pass

    def designColumn(self, only_mid = False):
        print("\n\n|--------------------------------------------------------|")
        print("Designing column...\n")
        print(f"Width: {self.width} cm, Height: {self.height} cm, fck: {self.fck*10} MPa, l_e: {self.l_e} cm")
        print(f"Cover: {self.cover} cm, Stirrups: {self.stirrups} cm\n")
        print(f"Axial Load N: {self.N}")
        print(f"Moments Mx: {self.Mx}")
        print(f"Moments My: {self.My}\n")

        vertices = np.array([[0, 0], [self.width, 0], [self.width, self.height], [0, self.height], [0, 0]])
        
        
        reinf_bars = []
        n_bars = 2

        size = np.max([self.width, self.height])
        bar_each = 2

        while ( (n_bars-1) * bar_each) < (size - 2 * (self.cover + self.stirrups + 1)):
            n_bars += 1

        bar_each = float((size - 2 * (self.cover + self.stirrups + 1)) / (n_bars - 1))
        bar_center = self.cover + self.stirrups + 1

        for i in range(n_bars):
            pos = self.cover + self.stirrups + 1 + i * bar_each
            if self.width > self.height:
                reinf_bars.append( [pos, bar_center, 1.0] )
                reinf_bars.append( [pos, self.height - bar_center, 1.0] )
            else:
                reinf_bars.append( [bar_center, pos, 1.0] )
                reinf_bars.append( [self.width - bar_center, pos, 1.0] )
        reinf_bars = np.array(reinf_bars)

        print(f"Number of bars per face: {n_bars}\n\n")
        section1 = ConcreteSection(vertices, reinf_bars, self.fck, consider_eta_c = self.consider_eta_c)

        As_top = 0
        As_bottom = 0
        if not only_mid:
            print("Designing reinforcement for top section...")
            As_top = section1.design_section(self.N['top'], self.Mx['top'], self.My['top'], iprint=False)
            print(f'Final reinforcement area at top: {As_top:.2f} cm2')
            print("Designing reinforcement for bot section...")
            As_bottom = section1.design_section(self.N['bottom'], self.Mx['bottom'], self.My['bottom'], iprint=False)
            print(f'Final reinforcement area at bottom: {As_bottom:.2f} cm2')
        print("Designing reinforcement for mid section...")
        As_mid = section1.design_section(self.N['mid'], self.Mx['mid'], self.My['mid'], iprint=True)
        print(f'Final reinforcement area at mid: {As_mid:.2f} cm2')
        print('\n\n--------------------------------------------------------')

        return np.max([As_top, As_mid, As_bottom])
