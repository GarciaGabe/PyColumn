from funcs import read_file, clear_console
from concrete_section import concrete_section

clear_console()
print('Starting...')
fck, vertices, reinf_bars = read_file('datpil.txt')


section1 = concrete_section(vertices, reinf_bars, fck)
section1.plot_concrete_section()
diagram = section1.plot_verify_diagram(800)

#print(section1.calculate_resisting_moments(NSd = 100, angle = 0))

'''
VALORES DO DIAGRAMA ESTAO INCORRETOS
AINDA PRECISO CHECAR COMO ESTAO OS SINAIS, NAO ESTOU 100% CERTO.

'''