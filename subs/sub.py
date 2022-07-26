import os

energies = [7, 11, 19, 27, 39, 62]

for energy in energies:
    os.system(f'star-submit sub_ampt_python{energy}.xml')

print('donzo')
