import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.path import Path
from matplotlib.transforms import Bbox
from consts import *
matplotlib.use('TkAgg')

def concrete_characteristics(fck):
    ultimate_compressive_strain = (3.5 / 1000) if fck <= 50 else (2.6 / 1000 + 35 / 1000 * ((90 - fck*10)/100)**4)
    peak_compressive_strain = (2.0 / 1000) if fck <= 50 else (2.0 / 1000 + 0.085 / 1000 * (fck*10 - 50)**0.53)
    lambda_c = 0.80 if fck <= 50 else (0.80 - (fck*10-50)/400)
    alpha_c = 0.85 if fck <= 50 else (0.85 * ( 1 - ( fck*10 - 50 ) / 200 ))
    eta_c = 1.0 if (not CONSIDER_ETA_C or fck <= 40 ) else (fck*10)**(2/3)
    
    dic = {
        'fck': fck,
        'fc': 0.95 * fck / GAMMA_C * eta_c * alpha_c,
        'eta_c': eta_c,
        'peak_compressive_strain': peak_compressive_strain,
        'ultimate_compressive_strain': ultimate_compressive_strain,
        'eta_c': eta_c,
        'lambda_c': lambda_c,
        'alpha_c': alpha_c
    }
    return dic 
    
def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')

def read_file(file_name):

    # Program start: reading vertices
    lines = []

    with open(file_name, 'r') as file:
        for line in file:
            if line.startswith('//'):
                continue
            lines.append(line.strip())

    fck = float(lines[0])
    line_index = 1

    n_vertices = int(lines[line_index])
    vertices = np.zeros((n_vertices, 2))
    line_index += 1

    for i in range(n_vertices):
        vertices[i] = lines[line_index].split()
        line_index += 1

    vertices = np.array(vertices, dtype=float)


    # Reading reinforcing bars

    n_reinf_bars = int(lines[line_index])
    line_index += 1
    reinf_bars = np.zeros((n_reinf_bars, 3))


    for i in range(n_reinf_bars):
        reinf_bars[i] = lines[line_index].split()
        line_index += 1

    return fck, vertices, reinf_bars

def rotate_axes(vertices, reinf_bars, angle):
    rotated_vertices = np.zeros_like(vertices)
    rotated_reinf_bars = np.zeros_like(reinf_bars)

    for i in range(len(vertices)):
        x, y = vertices[i]
        x_line =  x * np.cos(np.radians(angle)) + y * np.sin(np.radians(angle))
        y_line = - x * np.sin(np.radians(angle)) + y * np.cos(np.radians(angle))
        rotated_vertices[i] = [x_line, y_line]

    for i in range(len(reinf_bars)):
        x, y, steel_area = reinf_bars[i]
        x_line =  x * np.cos(np.radians(angle)) + y * np.sin(np.radians(angle))
        y_line = - x * np.sin(np.radians(angle)) + y * np.cos(np.radians(angle))
        rotated_reinf_bars[i] = [x_line, y_line, steel_area]

    return rotated_vertices, rotated_reinf_bars

def calculate_heights(vertices):
    y_max = np.max(vertices[:, 1])
    y_min = np.min(vertices[:, 1])
    height = y_max - y_min
    return y_max, y_min, height

def calculate_effective_heights(y_max, reinf_bars):
    effective_heights = np.zeros(len(reinf_bars))
    for i in range(len(reinf_bars)):
        effective_heights[i] = y_max - reinf_bars[i][1]
    return effective_heights

def steel_stress_strain_curve(steel_strain, young_modulus = 21000, yield_stress = 50):
    design_yield_stress = yield_stress / GAMMA_S
    yield_strain = design_yield_stress / young_modulus
    stress = np.where(np.abs(steel_strain) <= yield_strain,
                      young_modulus * steel_strain,
                      design_yield_stress * np.sign(steel_strain))
    return stress

def calculate_reinf_tensions(x0 ,effective_heights, reinf_bars,peak_compressive_strain , ultimate_compressive_strain, height):
    reinf_strains = np.zeros(len(reinf_bars))
    reinf_stresses = np.zeros(len(reinf_bars))
    effective_height = np.max(effective_heights)
    
    domain_2_limit = (ultimate_compressive_strain / (ultimate_compressive_strain + 10 / 1000)) * effective_height
    domain_34_limit = height

    if (x0 < 0):
        raise Exception("Concrete section is not under bending compression")
    elif (x0 <= domain_2_limit):
        reinf_strains = 10/1000 * ( ( x0 - effective_heights ) / ( effective_height - x0 ) )
    elif (x0 <= domain_34_limit):
        reinf_strains = ultimate_compressive_strain * ( (x0 - effective_heights) / x0 )
    else:
        kappa = 1 - peak_compressive_strain / ultimate_compressive_strain
        reinf_strains = peak_compressive_strain * ( ( x0 - effective_heights ) / ( x0 - kappa * height ) )

    reinf_stresses = steel_stress_strain_curve(reinf_strains, STEEL_YOUNG_MODULUS, STEEL_YIELD_STRESS)
    return reinf_strains, reinf_stresses
    
def calculate_geometric_center(vertices):
    n = len(vertices)
    if n < 3:
        if n == 0:
            center_x_concrete = 0.0
            center_y_concrete = 0.0
        elif n == 1:
            center_x_concrete = vertices[0, 0]
            center_y_concrete = vertices[0, 1]
        else: # n == 2
            center_x_concrete = np.mean(vertices[:, 0])
            center_y_concrete = np.mean(vertices[:, 1])
    else:
        closed_vertices = np.vstack([vertices, vertices[0]])

        area_sum = 0.0
        for i in range(n):
            area_sum += (closed_vertices[i, 0] * closed_vertices[i+1, 1] -
                         closed_vertices[i+1, 0] * closed_vertices[i, 1])
        area = 0.5 * area_sum

        if area == 0:
            center_x_concrete = np.mean(vertices[:, 0])
            center_y_concrete = np.mean(vertices[:, 1])
        else:
            cx_sum = np.sum((closed_vertices[:-1, 0] + closed_vertices[1:, 0]) * (closed_vertices[:-1, 0] * closed_vertices[1:, 1] - closed_vertices[1:, 0] * closed_vertices[:-1, 1]))
            cy_sum = np.sum((closed_vertices[:-1, 1] + closed_vertices[1:, 1]) * (closed_vertices[:-1, 0] * closed_vertices[1:, 1] - closed_vertices[1:, 0] * closed_vertices[:-1, 1]))
            center_x_concrete = (1 / (6 * area)) * cx_sum
            center_y_concrete = (1 / (6 * area)) * cy_sum

    geometric_center = np.array([center_x_concrete, center_y_concrete])
    return geometric_center

def translate_to_geometric_center(vertices, reinf_bars):
    geometric_center = calculate_geometric_center(vertices)
    
    translated_vertices = vertices - geometric_center
    
    translated_reinf_bars = np.copy(reinf_bars)
    translated_reinf_bars[:, 0] = reinf_bars[:, 0] - geometric_center[0]
    translated_reinf_bars[:, 1] = reinf_bars[:, 1] - geometric_center[1]
    
    return translated_vertices, translated_reinf_bars

def find_compression_vertices(vertices, x0, concrete):

    y_max = np.max(vertices[:, 1])
    y_c = y_max - concrete['lambda_c'] * x0

    section_path = Path(vertices, closed=True)
    v_min, v_max = np.min(vertices, axis=0), np.max(vertices, axis=0)
    padding = max(v_max[0] - v_min[0], v_max[1] - v_min[1]) * 2
    clip_box = Bbox.from_extents(v_min[0] - padding, y_c, v_max[0] + padding, v_max[1] + padding)
    clipped_path = section_path.clip_to_bbox(clip_box, inside=True)

    return clipped_path.vertices

def plot_section(vertices, reinf_bars, compressed_vertices = None, bar_tensions = None):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(vertices[:, 0], vertices[:, 1], color='black', linewidth = 2, label='Vertices')
    ax.fill(vertices[:, 0], vertices[:, 1], color='gray', alpha = 0.20)

    if compressed_vertices is not None:
        ax.fill(compressed_vertices[:, 0], compressed_vertices[:, 1], label='Compressed Vertices')

    # The 's' parameter controls marker size (area in points^2).
    # Using it instead of 'linewidths' to correctly represent rebar area.
    # A scaling factor is used for better visualization.
    if bar_tensions is None:
        ax.scatter(reinf_bars[:, 0], reinf_bars[:, 1], s = 50 * reinf_bars[:, 2], color='red', label='Reinforcing Bars')
    else:
        ax.scatter(reinf_bars[:, 0], reinf_bars[:, 1], s = 50 * reinf_bars[:, 2], cmap = bar_tensions, label='Reinforcing Bars')
    

    lim = np.max(np.abs(vertices))*1.05
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal', adjustable='box')
    return fig, ax

def calculate_polygon_properties(compressed_vertices):
    """
    Calculates the area and first moments of inertia (Mx, My) for a given polygon.

    Parameters
    ----------
    compressed_vertices (np.ndarray): A 2D NumPy array of shape (N, 2)
                                      representing the polygon's vertices.

    Returns
    -------
    tuple
           - area (float): The area of the polygon.
           - sx (float): The first moment of area about the x-axis.
           - sy (float): The first moment of area about the y-axis.
    """
    n = len(compressed_vertices)
    if n < 3:
        return 0.0, 0.0, 0.0 # Area, Sx, Sy

    # Close the polygon by appending the first vertex to the end
    closed_vertices = np.vstack([compressed_vertices, compressed_vertices[0]])

    area_sum = 0.0
    sx_sum = 0.0
    sy_sum = 0.0

    for i in range(n):
        xi, yi = closed_vertices[i]
        xi1, yi1 = closed_vertices[i+1]

        cross_product_term = (xi * yi1) - (xi1 * yi)
        area_sum += cross_product_term
        sx_sum += (yi + yi1) * cross_product_term
        sy_sum += (xi + xi1) * cross_product_term

    area = 0.5 * area_sum
    sx = (1/6) * sx_sum
    sy = (1/6) * sy_sum

    if area < 0:
        # The shoelace formula gives a negative area for clockwise ordering.
        # For physical consistency, we expect a counter-clockwise ordering.
        area = abs(area)
    return area, sx, sy

def calculate_axial_and_moments(compressed_vertices, reinf_bars, reinf_tensions, concrete):
    area, Sx, Sy = calculate_polygon_properties(compressed_vertices)

    as_sum = np.zeros(3) # Force, Mx, My
    for bar, tension in zip(reinf_bars, reinf_tensions[1]):
        x, y, steel_area = bar
        as_sum[0] += steel_area * tension
        as_sum[1] += steel_area * y * tension
        as_sum[2] += steel_area * x * tension

    Nd = area*concrete['fc'] + as_sum[0]
    Mx = Sx*concrete['fc'] + as_sum[1]
    My = Sy*concrete['fc'] + as_sum[2]

    return Nd, Mx, My
