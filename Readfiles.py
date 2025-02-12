import skrf as rf
from pathlib import Path
import os
import re
from itertools import combinations
import shutil

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
                    self.comparison.network.interpolate(self.standard.network.f)
                except:
                    self.standard.network.interpolate(self.comparison.network.f)
        else:
            print(f"{self.comparison.software} is not exist: {key}")
    
    def snp_s2p(self, key:str, result_path:Path):
        if not result_path.exists():
            os.makedirs(result_path)
            print(f"Directory '{result_path}' created.")
        else:
            shutil.rmtree(result_path)
        
        combos = list(combinations(range(self.compare_list[key][2]), 2))
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(len(combos), 4, figsize=(5*4, 3*len(combos)))
        for j, combo in enumerate(combos):
            # get the subnetwork for the given port combination when the other port is matched to 50 ohms
            standard_s2p = rf.subnetwork(self.standard.network, combo)
            comparison_s2p = rf.subnetwork(self.comparison.network, combo)
            # plot the S-parameters
            for i, permutation in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
                ax[j][i].set_title(f"S{combo[permutation[0]]+1}{combo[permutation[1]]+1}")
                ax[j][i].set_xlabel("Frequency (GHz)")
                ax[j][i].set_ylabel("S-parameter")
                ax[j][i].grid(True)
                standard_s2p.plot_s_db(m=permutation[0], n=permutation[1], ax=ax[j][i], label=f"{key}{self.standard.software}")
                comparison_s2p.plot_s_db(m=permutation[0], n=permutation[1], ax=ax[j][i], label=f"{key}{self.comparison.software}")
                ax[j][i].legend()
        (result_path /"Figures").mkdir(parents=True, exist_ok=True)
        plt.tight_layout()
        fig.savefig(result_path /"Figures"/f"{key}.png")
        plt.close(fig)

    def compare_all(self, result_path = Path("./Results")):
        for key in self.compare_list:
            self.read(key)
            self.snp_s2p(key, result_path)

if __name__ == "__main__":
    compare_list = {}
    cst = S_parameter(Path("./Standard"), "_cst")
    tidy3d = S_parameter(Path("./Comparison"), "_tidy3d")
    cmp = compare(cst, tidy3d, compare_list)
    cmp.compare_all()