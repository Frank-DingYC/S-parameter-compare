import skrf as rf
from pathlib import Path
import os
import re
from itertools import combinations
import shutil
import numpy as np

class S_parameter:
    def __init__(self, path:Path, software:str):
        self.software = software.lower()
        if path.exists():
            temp_list = os.listdir(path)  # Get the list of files in the path
            for file in temp_list:
                suffix = file.split('.')[-1].lower()  # Extract the file extension
                port_number = re.search(r'\d+', suffix)  # Check for any digits in the suffix to determine port number
                if re.search(r'[sp]', suffix) and port_number:  # If the file is an s-parameter file
                    global compare_list  # Declare global variable compare_list
                    key = file.split(suffix)[0].lower().split(self.software)[0]  # Extract the base name without extension and software part
                    if key in compare_list:
                        compare_list[key][1] = path / file  # Update existing entry with new path
                    else:
                        compare_list[key] = [path / file, None, int(port_number[0])]  # Create new entry with port number
        else:
            raise FileNotFoundError("The path is not exist")
    
    def read(self, file):
        self.network = rf.Network(file)

class compare:

    def __init__(self, standard:S_parameter, comparison:S_parameter, compare_list:dict):
        self.standard = standard
        self.comparison = comparison
        self.compare_list = compare_list
    
    def read(self, key):
        if self.compare_list[key][1] is not None:
            self.standard.read(self.compare_list[key][0])
            self.comparison.read(self.compare_list[key][1])
            if not (len(self.standard.network.f) == len(self.comparison.network.f) and \
               (self.standard.network.f == self.comparison.network.f).all()):
                Warning(f"The frequency of {self.standard.software} and {self.comparison.software} are not the same")
                try:
                    self.comparison.network = self.comparison.network.interpolate(self.standard.network.f)
                except ValueError:
                    self.standard.network = self.standard.network.interpolate(self.comparison.network.f)
        else:
            print(f"{self.comparison.software} is not exist: {key}")
    
    def snp_s2p(self, key:str, result_path:Path, cmp_type:str):
        combos = list(combinations(range(self.compare_list[key][2]), 2))
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(len(combos), 4, figsize=(5*4, 3*len(combos)))
        fig1, ax1 = plt.subplots(len(combos), 4, figsize=(5*4, 3*len(combos)))
        for j, combo in enumerate(combos):
            # get the subnetwork for the given port combination when the other port is matched to 50 ohms
            standard_s2p = rf.subnetwork(self.standard.network, combo)
            comparison_s2p = rf.subnetwork(self.comparison.network, combo)
            # plot the S-parameters
            for i, permutation in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
                self.plot_s_parameters(key, ax, i, j, standard_s2p, comparison_s2p, combo, permutation)
                self.compare_s2p(ax1, i, j, standard_s2p, comparison_s2p, combo, permutation, cmp_type)
        (result_path /"Figures").mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig1.tight_layout()
        fig.savefig(result_path /"Figures"/f"{key}_S-parameters.png")
        fig1.savefig(result_path /"Figures"/f"{key}_{cmp_type}_diff.png")
        plt.close(fig)
        plt.close(fig1)

    
    def compare_s2p(self, ax, i, j, standard_s2p, comparison_s2p, combo, permutation, cmp_type):
        frequency = standard_s2p.frequency.f
        if cmp_type == 'abs':
            standard = abs(standard_s2p.s[:, permutation[0], permutation[1]])
            comparison = abs(comparison_s2p.s[:, permutation[0], permutation[1]])
            # calculate the relative difference
            relative_difference  = abs(standard - comparison) / standard * 100
            # plot the difference
            try:
                ax[j][i].plot(frequency, relative_difference, label='Difference')
                ax[j][i].set_title(f"Abs Relative Difference S{combo[permutation[0]]+1}{combo[permutation[1]]+1}")
            except TypeError:
                ax[i].plot(frequency, relative_difference, label='Difference')
                ax[i].set_title(f"Abs Relative Difference S{combo[permutation[0]]+1}{combo[permutation[1]]+1}")

        elif cmp_type == 'db':
            standard = 20*np.log10(np.abs(standard_s2p.s[:, permutation[0], permutation[1]]))
            comparison = 20*np.log10(np.abs(comparison_s2p.s[:, permutation[0], permutation[1]]))
            # calculate the relative difference
            relative_difference  = abs(standard - comparison) / abs(standard) * 100
            # plot the difference
            try:
                ax[j][i].plot(frequency, relative_difference, label='Difference')
                ax[j][i].set_title(f"dB Relative Difference S{combo[permutation[0]]+1}{combo[permutation[1]]+1}")
            except TypeError:
                ax[i].plot(frequency, relative_difference, label='Difference')
                ax[i].set_title(f"dB Relative Difference S{combo[permutation[0]]+1}{combo[permutation[1]]+1}")
        elif cmp_type == 'complex':
            standard = standard_s2p.s[:, permutation[0], permutation[1]]
            comparison = comparison_s2p.s[:, permutation[0], permutation[1]]
            # calculate the relative difference
            relative_difference  = ((standard.real-comparison.real)**2 + (standard.imag-comparison.imag)**2)**0.5 / abs(standard) * 100
            # plot the difference
            try:
                ax[j][i].plot(frequency, relative_difference, label='Difference')
                ax[j][i].set_title(f"Complex Relative Difference S{combo[permutation[0]]+1}{combo[permutation[1]]+1}")
            except TypeError:
                ax[i].plot(frequency, relative_difference, label='Difference')
                ax[i].set_title(f"Complex Relative Difference S{combo[permutation[0]]+1}{combo[permutation[1]]+1}")
        else:
            raise ValueError('cmp_type must be abs, db or complex')
        try:
            ax[j][i].plot(frequency, 5*np.ones_like(frequency), 'r--', label='5% Difference') 
            ax[j][i].set_xlabel("Frequency (GHz)")
            ax[j][i].set_ylabel("Difference (%)")
            ax[j][i].legend()
        except TypeError:
            ax[i].plot(frequency, 5*np.ones_like(frequency), 'r--', label='5% Difference')
            ax[i].set_xlabel("Frequency (GHz)")
            ax[i].set_ylabel("Difference (%)")
            ax[i].legend()

    def plot_s_parameters(self, key, ax, i, j, standard_s2p, comparison_s2p, combo, permutation):
        try:
            ax[j][i].set_title(f"S{combo[permutation[0]]+1}{combo[permutation[1]]+1}")
            ax[j][i].set_xlabel("Frequency (GHz)")
            ax[j][i].set_ylabel("S-parameter")
            ax[j][i].grid(True)
            standard_s2p.plot_s_db(m=permutation[0], n=permutation[1], ax=ax[j][i], label=f"{key}{self.standard.software}")
            comparison_s2p.plot_s_db(m=permutation[0], n=permutation[1], ax=ax[j][i], label=f"{key}{self.comparison.software}")
            ax[j][i].legend()
        except TypeError:
            ax[i].set_title(f"S{combo[permutation[0]]+1}{combo[permutation[1]]+1}")
            ax[i].set_xlabel("Frequency (GHz)")
            ax[i].set_ylabel("S-parameter")
            ax[i].grid(True)
            standard_s2p.plot_s_db(m=permutation[0], n=permutation[1], ax=ax[i], label=f"{key}{self.standard.software}")
            comparison_s2p.plot_s_db(m=permutation[0], n=permutation[1], ax=ax[i], label=f"{key}{self.comparison.software}")
            ax[i].legend()

    def compare_all(self, result_path = Path("./Results"), cmp_type = 'db'):
        if not result_path.exists():
            os.makedirs(result_path)
            print(f"Directory '{result_path}' created.")
        else:
            shutil.rmtree(result_path)
            print(f"Directory '{result_path}' already exists. Deleting and recreating.")
            os.makedirs(result_path)
            print(f"Directory '{result_path}' recreated.")
        for key in self.compare_list:
            self.read(key)
            self.snp_s2p(key, result_path, cmp_type)

if __name__ == "__main__":
    compare_list = {}
    cst = S_parameter(Path("./Standard"), "_cst")
    tidy3d = S_parameter(Path("./Comparison"), "_tidy3d")
    cmp = compare(cst, tidy3d, compare_list)
    cmp.compare_all(cmp_type='abs')