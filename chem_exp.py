import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
from scipy.signal import find_peaks

def plot_maldi_data(file_paths, legends, x_section, svg_path):
    # Setting up the number of subplots based on the number of files
    fig, axs = plt.subplots(len(file_paths), 1, sharex=True, figsize=(10,8))
    if len(file_paths) == 1:
        axs = [axs]  # Ensure axs is a list for consistency
    
    # Eliminating margins and adjusting space between subplots
    plt.subplots_adjust(hspace=0, top=1, bottom=0, right=1, left=0)
    bsa_mass_list = [2435.2427, 2003.1000, 1955.9596, 1888.9268, 1850.8993, 1823.8996, 1667.8131, 1633.6621, 1578.5981, 1567.7427, 1519.7461, 1511.8427, 1497.6314, 1479.7954, 1399.6926, 1388.5708, 1386.6206, 1364.4803, 1362.6722, 1349.5460, 1305.7161, 1283.7106, 1177.5591, 1163.6306, 1052.4499, 1050.4924, 1024.4550, 1015.4877, 1014.6193, 1011.4200, 1002.5830]
    albumin_mass_list = [3863.0108, 3293.6279, 2976.4932, 2460.3166, 2374.1455, 2284.1464, 2281.1822, 2008.9457, 1858.9657, 1839.8114, 1773.8990, 1687.8398, 1581.7213, 1465.7759, 1345.7375, 1209.5204, 1190.6026]
    hemoglobin_mass_list = [2969.6094, 2367.1938, 1833.8918, 1529.7342, 1279.7256, 1071.5543]
    
    # Processing each file and plotting
    for i, file in enumerate(file_paths):
        try:
            # Read the data
            data = pd.read_csv(file, sep=' ', header=None, names=["x", "y"])

            # Skip empty data files
            if data.empty:
                print(f"Warning: {file} is empty.")
                continue
            
            data_len = data.shape[0]
            # data = data.rolling(window=30).mean()
            x_ = data["x"]
            y_ = data["y"]
            
            x_range = (data["x"] > x_section[0]) & (data["x"] < x_section[1])
            ## Softmax ##
            # Maximize signal-to-noise ratio 
            # y_max = np.max(np.array(y_[x_range]))
            # y_ = np.exp(y_)
            # y_ = np.divide(y_, np.exp(y_max))
            
            # Find all y_ local maximum
            peak_indices, _ = find_peaks(y_)
            x_val = []
            for j in range(len(x_)):
                if j in peak_indices:
                    x_val.append(x_[j])

            # Plot the data in a subplot
            axs[i].plot(x_, y_, label=legends[i], color='black', linewidth=1)
            
            # Select useful y_ local maximum and x_val in that point
            real_peak_x_val = []
            for j in range(len(x_val)):
                if x_[peak_indices[j]] > x_section[0] and x_[peak_indices[j]] < x_section[1]:
                    if y_[peak_indices[j]] > 0.1:
                        real_peak_x_val.append(x_[peak_indices[j]])
            
            y_range = data[x_range]["y"]
            
            # Visualize useful y_ local maximum points
            for j in range(len(real_peak_x_val)):
                axs[i].vlines(real_peak_x_val[j], 0, max(y_range), color="r", ls=":", label="{0}".format(real_peak_x_val[j]))

            axs[i].legend()
            axs[i].set_xlim(x_section)
            
            if not y_range.empty:
                ## Data Overview ##
                axs[i].set_ylim([0, max(y_range)]) 

                ## Range Analysis (Softmax version) ##
                # axs[i].set_ylim([0, 1])
            else:
                axs[i].set_ylim([0, 0.1])  # Or a default range if y_range is empty

            # Remove x-axis ticks and labels for all but the bottom plot
            if i < len(file_paths) - 1:
                axs[i].tick_params(axis='x', which='both', bottom=False, labelbottom=False)

            # Remove all box lines except the bottom (x-axis)
            axs[i].spines['top'].set_visible(False)
            axs[i].spines['right'].set_visible(False)
            axs[i].spines['left'].set_visible(False)

            # For the bottom subplot, only keep the x-axis
            if i == len(file_paths) - 1:
                axs[i].spines['bottom'].set_linewidth(1)
                axs[i].set_xlabel('M/Z')
            else:
                axs[i].spines['bottom'].set_visible(True)

        except Exception as e:
            print(f"Error processing {file}: {e}")

    # Adjust layout and save the figure
    plt.tight_layout()

    plt.savefig(svg_path, format='svg')
    plt.show()
    
    # To avoid multiple checking
    observed_data_done_bsa = []
    observed_data_done_albumin = []
    observed_data_done_hemoglobin = []
    for i in range(len(real_peak_x_val)):
        observed_data_done_bsa.append(0)
        observed_data_done_albumin.append(0)
        observed_data_done_hemoglobin.append(0)
    real_data_done_bsa = []
    real_data_done_albumin = []
    real_data_done_hemoglobin = []
    for i in range(len(bsa_mass_list)):
        real_data_done_bsa.append(0)
    for i in range(len(albumin_mass_list)):
        real_data_done_albumin.append(0)
    for i in range(len(hemoglobin_mass_list)):
        real_data_done_hemoglobin.append(0)
    
    bsa_point = 0
    albumin_point = 0
    hemoglobin_point = 0
    for one_peptide_mass_idx in range(len(bsa_mass_list)):
        for observed_mass_idx in range(len(real_peak_x_val)):
            if real_peak_x_val[observed_mass_idx] < bsa_mass_list[one_peptide_mass_idx] + 0.5 and real_peak_x_val[observed_mass_idx] > bsa_mass_list[one_peptide_mass_idx] - 0.5 and observed_data_done_bsa[observed_mass_idx] == 0 and real_data_done_bsa[one_peptide_mass_idx] == 0:
                bsa_point += 1
                observed_data_done_bsa[observed_mass_idx] = 1
                real_data_done_bsa[one_peptide_mass_idx] = 1
    for one_peptide_mass_idx in range(len(albumin_mass_list)):
        for observed_mass_idx in range(len(real_peak_x_val)):
            if real_peak_x_val[observed_mass_idx] < albumin_mass_list[one_peptide_mass_idx] + 0.5 and real_peak_x_val[observed_mass_idx] > albumin_mass_list[one_peptide_mass_idx] - 0.5 and observed_data_done_albumin[observed_mass_idx] == 0 and real_data_done_albumin[one_peptide_mass_idx] == 0:
                albumin_point += 1
                observed_data_done_albumin[observed_mass_idx] = 1
                real_data_done_albumin[one_peptide_mass_idx] = 1
    for one_peptide_mass_idx in range(len(hemoglobin_mass_list)):
        for observed_mass_idx in range(len(real_peak_x_val)):
            if real_peak_x_val[observed_mass_idx] < hemoglobin_mass_list[one_peptide_mass_idx] + 0.5 and real_peak_x_val[observed_mass_idx] > hemoglobin_mass_list[one_peptide_mass_idx] - 0.5 and observed_data_done_hemoglobin[observed_mass_idx] == 0 and real_data_done_hemoglobin[one_peptide_mass_idx] == 0:
                hemoglobin_point += 1
                observed_data_done_hemoglobin[observed_mass_idx] = 1
                real_data_done_hemoglobin[one_peptide_mass_idx] = 1

    print("============ Result ============")
    print("==== Verified Peptide ratio ====")
    print("BSA : ", bsa_point, "/", len(bsa_mass_list), ", ", round(bsa_point/len(bsa_mass_list),4) * 100, "%")
    print("Albumin : ", albumin_point, "/", len(albumin_mass_list), ", ", round(albumin_point/len(albumin_mass_list), 5) * 100, "%")
    print("Hemoglobin : ", hemoglobin_point, "/", len(hemoglobin_mass_list), ", ", round(hemoglobin_point/len(hemoglobin_mass_list), 4) * 100, "%")
    print("================================")

input=["/Path/To/Input/Data.txt"]
plot_maldi_data(input,"experimental data",[1050,1120], "/Path/To/Output/File.svg")