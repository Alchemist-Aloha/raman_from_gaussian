import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path



def plot_raman_spectrum(min, max, nstep, freqs, intensities, lorentzian_width,
                             gaussian_width, lorentzian_percentage):
    frequencies = np.linspace(min, max, nstep)
    central_frequencies = freqs
    raman_intensities = intensities

    lorentzian_lineshapes = np.zeros((len(central_frequencies), len(frequencies)))
    gaussian_lineshapes = np.zeros((len(central_frequencies), len(frequencies)))

    for i, (central_frequency, intensity) in enumerate(zip(central_frequencies, raman_intensities)):
        lorentz = (1 / (1 + ((frequencies - central_frequency) / lorentzian_width)**2))
        lorentz /= np.trapezoid(lorentz, frequencies)  # Normalize area
        lorentzian_lineshapes[i, :] = lorentz * intensity

        gauss = (np.exp(-0.5 * ((frequencies - central_frequency) / gaussian_width)**2) /
                 (gaussian_width * np.sqrt(2 * np.pi)))
        gauss /= np.trapezoid(gauss, frequencies)  # Normalize area
        gaussian_lineshapes[i, :] = gauss * intensity

    gaussian_spectrum = np.sum(gaussian_lineshapes, axis=0)
    lorentzian_spectrum = np.sum(lorentzian_lineshapes, axis=0)
    output = (1 - lorentzian_percentage) * gaussian_spectrum + lorentzian_percentage * lorentzian_spectrum

    plt.figure(figsize=(5, 4))
    plt.plot(frequencies, lorentzian_spectrum, "r-", label="Lorentzian")
    plt.plot(frequencies, gaussian_spectrum, "b-", label="Gaussian")
    plt.plot(frequencies, output, "g-", label="Combined")
    plt.xlabel("Frequency")
    plt.ylabel("Raman Intensity")
    plt.title("Raman Spectrum: Lorentzian vs Gaussian")
    plt.legend()
    plt.grid(True)
    plt.show()
    return frequencies, gaussian_spectrum, lorentzian_spectrum, output


def generate_rr_spectrum(min, max, nstep, input_filename, lorentzian_width,
                            gaussian_width, lorentzian_percentage):
    """Generates a Raman spectrum from gaussian log file and saves it to a CSV file."""
    input_filepath = Path(input_filename)
    filestem = input_filepath.stem
    filedir = input_filepath.parent
    print(filestem)
    frequencies = np.linspace(min, max, nstep)
    data_array = extract_rr(input_filename)
    central_frequencies = data_array[1:, 0]
    raman_intensities = data_array[1:, 1]
    np.savetxt(
    filedir / (filestem + "_rrcross.csv"),
    data_array[1:, :], 
    delimiter=",", 
    header="Raman Shift (cm^-1),Sigma", 
    comments='', 
    fmt='%.8e' # Specifies scientific notation output for high precision
)

    lorentzian_lineshapes = np.zeros((len(central_frequencies), len(frequencies)))
    gaussian_lineshapes = np.zeros((len(central_frequencies), len(frequencies)))

    for i, (central_frequency, intensity) in enumerate(zip(central_frequencies, raman_intensities)):
        lorentz = (1 / (1 + ((frequencies - central_frequency) / lorentzian_width)**2))
        lorentz /= np.trapezoid(lorentz, frequencies)
        lorentzian_lineshapes[i, :] = lorentz * intensity

        gauss = (np.exp(-0.5 * ((frequencies - central_frequency) / gaussian_width)**2) /
                 (gaussian_width * np.sqrt(2 * np.pi)))
        gauss /= np.trapezoid(gauss, frequencies)
        gaussian_lineshapes[i, :] = gauss * intensity

    gaussian_spectrum = np.sum(gaussian_lineshapes, axis=0)
    lorentzian_spectrum = np.sum(lorentzian_lineshapes, axis=0)
    output = (1 - lorentzian_percentage) * gaussian_spectrum + lorentzian_percentage * lorentzian_spectrum
    plt.figure(figsize=(5, 4))
    plt.plot(frequencies, lorentzian_spectrum, "r-", label="Lorentzian")
    plt.plot(frequencies, gaussian_spectrum, "b-", label="Gaussian")
    plt.plot(frequencies, output, "g-", label="Combined")
    plt.xlabel("Frequency")
    plt.ylabel("Raman Intensity")
    plt.title("Raman Spectrum: Lorentzian vs Gaussian")
    plt.legend()
    plt.grid(True)
    plt.show()
    np.savetxt(filedir / (filestem + f"_rrspec_lorentzian{lorentzian_width}.csv"), np.column_stack((frequencies, lorentzian_spectrum)), delimiter=',', header='Raman Shift (cm^-1),Raman Intensity', comments='')
    np.savetxt(filedir / (filestem + f"_rrspec_gaussian{gaussian_width}.csv"), np.column_stack((frequencies, gaussian_spectrum)), delimiter=',', header='Raman Shift (cm^-1),Raman Intensity', comments='')
    # optionally save combined spectrum
    # np.savetxt(filedir / (filestem + f"_rrspec_lorentzian{lorentzian_width}_gaussian{gaussian_width}_lorentzianweight{lorentzian_percentage}.csv"), np.column_stack((frequencies, output)), delimiter=',', header='Raman Shift (cm^-1),Raman Intensity', comments='')
    return frequencies, gaussian_spectrum, lorentzian_spectrum, output



# --- 1. Data Setup (omitted for brevity, assumes raman_data.txt is created) ---

def extract_rr(file_path):
    """
    Extracts Omega and Sigma values using re and converts them directly to a
    NumPy array, avoiding the use of Pandas.
    """
    # Read the entire content of the file
    with open(file_path, 'r') as file:
        text = file.read()
    # Initialize lists to hold the extracted energies and sigmas
    energies = []
    sigmas = []

    # A flag to indicate whether we are in the correct section to extract data
    rr = False
    for line in text.splitlines():
        # Check for the start marker
        if "Information on Transitions" in line:
            rr = True
            continue  # Move to the next line without processing this one

        # Check for the end marker
        if "Final Spectrum" in line:
            rr = False
            break # Exit the loop as we are done

        # If we are in the correct section, extract the data
        if rr:
            if "Energy =" in line:
                match = re.search(r"Energy =\s*([\d\.]+)", line)
                if match:
                    energies.append(float(match.group(1)))
            elif "Sigma =" in line:
                match = re.search(r"Sigma =\s*([-\d\.E\+]+)", line)
                if match:
                    sigmas.append(float(match.group(1)))

    print("Energies:")
    print(energies)
    print("Sigmas:")
    print(sigmas)
    return np.column_stack((energies, sigmas))


def extract_nrr(file_path):
    """
    Extracts frequencies and Raman activities for all incident light frequencies from Gaussian output.
    Returns a NumPy array: columns are [frequency, RamAct Fr=1, RamAct Fr=2, ...].
    Dynamically handles any number of incident light frequencies.
    """
    with open(file_path, 'r') as file:
        text = file.read()

    # Extract incident light frequencies - handle variable number
    incident_light_pattern = r"Incident light \(cm\*\*-1\):\s+((?:[\d]+\.[\d]+\s*)+)"
    incident_light_match = re.search(incident_light_pattern, text)
    incident_light = []
    if incident_light_match:
        # Split the matched string and convert to floats
        incident_light = [float(x) for x in incident_light_match.group(1).split()]
        print(f"Incident light frequencies: {incident_light} cm⁻¹")
    else:
        # Fallback pattern for cases where numbers might not have decimal points
        incident_light_pattern_fallback = r"Incident light \(cm\*\*-1\):\s+((?:[\d]+\.[\d]{2}\s*)+)"
        incident_light_match = re.search(incident_light_pattern_fallback, text)
        if incident_light_match:
            incident_light = [float(x) for x in incident_light_match.group(1).split()]
            print(f"Incident light frequencies: {incident_light} cm⁻¹")
        else:
            print("Warning: Could not find incident light frequencies")
            return np.array([])

    num_incident_lights = len(incident_light)
    print(f"Number of incident light frequencies: {num_incident_lights}")

    # Find all frequency blocks - each block contains up to 3 modes
    freq_pattern = r"Frequencies --\s+([-\d\.]+)(?:\s+([-\d\.]+))?(?:\s+([-\d\.]+))?"
    freq_matches = re.findall(freq_pattern, text)

    # Dynamically create patterns for all RamAct Fr= entries
    ramact_patterns = []
    ramact_matches_list = []
    
    for fr_num in range(1, num_incident_lights + 1):
        # Pattern that matches floating point numbers (including scientific notation)
        # number_pattern = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
        number_pattern = r"\S+"
        pattern = rf"RamAct Fr=\s*{fr_num}--\s*({number_pattern})(?:\s+({number_pattern}))?(?:\s+({number_pattern}))?"
        ramact_patterns.append(pattern)
        matches = re.findall(pattern, text)
        ramact_matches_list.append(matches)
        print(f"Found {len(matches)} blocks for RamAct Fr={fr_num}")

    # Debug: Check if we have consistent block counts
    expected_blocks = len(freq_matches)
    print(f"Expected {expected_blocks} blocks based on frequency data")
    
    for fr_idx in range(num_incident_lights):
        if len(ramact_matches_list[fr_idx]) != expected_blocks:
            print(f"WARNING: RamAct Fr={fr_idx+1} has {len(ramact_matches_list[fr_idx])} blocks, expected {expected_blocks}")

    # Collect all data
    frequencies = []
    ramact_data = [[] for _ in range(num_incident_lights)]  # List of lists for each Fr

    # Process each frequency block (each block has up to 3 modes)
    for i, freq_match in enumerate(freq_matches):
        for j, freq_str in enumerate(freq_match):
            if freq_str:  # Skip empty matches
                frequencies.append(float(freq_str))
                
                # Extract corresponding Raman activities for each incident light
                for fr_idx in range(num_incident_lights):
                    if (i < len(ramact_matches_list[fr_idx]) and 
                        j < len(ramact_matches_list[fr_idx][i]) and 
                        ramact_matches_list[fr_idx][i][j]):
                        try:
                            ramact_data[fr_idx].append(float(ramact_matches_list[fr_idx][i][j]))
                        except Exception as e:
                            print(f"Error converting RamAct Fr={fr_idx+1} value '{ramact_matches_list[fr_idx][i][j]}' to float: {e}")
                            ramact_data[fr_idx].append(0.0)
                    else:
                        # Debug: Print when we're using default values
                        if i >= len(ramact_matches_list[fr_idx]):
                            print(f"Missing block {i} for RamAct Fr={fr_idx+1}, using 0.0")
                        else:
                            print(f"Missing value at block {i}, position {j} for RamAct Fr={fr_idx+1}, using 0.0")
                        ramact_data[fr_idx].append(0.0)

    print(f"Extracted {len(frequencies)} frequencies")
    print(f"Frequencies: {frequencies[:5] if len(frequencies) >= 5 else frequencies}...")
    
    for fr_idx in range(num_incident_lights):
        print(f"RamAct Fr={fr_idx+1}: {ramact_data[fr_idx][:5] if len(ramact_data[fr_idx]) >= 5 else ramact_data[fr_idx]}...")
        print(f"Total RamAct Fr={fr_idx+1} values: {len(ramact_data[fr_idx])}")

    # Ensure all lists have the same length
    min_length = min(len(frequencies), min(len(data) for data in ramact_data))
    if min_length < len(frequencies):
        print(f"Warning: Truncating data to {min_length} entries due to mismatched lengths")
        frequencies = frequencies[:min_length]
        ramact_data = [data[:min_length] for data in ramact_data]

    # Return as NumPy array with columns: [frequency, RamAct Fr=1, Fr=2, ...]
    all_data = [frequencies] + ramact_data
    return np.column_stack(all_data), incident_light

def generate_nrr_spectrum(min, max, nstep, input_filename, lorentzian_width,
                            gaussian_width, lorentzian_percentage):
    """Generates a Raman spectrum from gaussian log file and saves it to a CSV file."""
    input_filepath = Path(input_filename)
    filestem = input_filepath.stem
    filedir = input_filepath.parent
    print(filestem)
    frequencies = np.linspace(min, max, nstep)
    data_array,incident_freqs = extract_nrr(input_filename)
    central_frequencies = data_array[:, 0]
    raman_intensities_array = data_array[:, 1:]
    header = "Raman Shift (cm^-1)"
    for i, freq in enumerate(incident_freqs):
        header += f",RamAct ({freq} cm^-1)"
    np.savetxt(
    filedir / (filestem + "_nrrcross.csv"),
    data_array[:, :], 
    delimiter=",", 
    header=header, 
    comments='', 
    fmt='%.8e' # Specifies scientific notation output for high precision
)
    lorentzian_list=[]
    gaussian_list=[]
    output_list=[]
    for i,freq in enumerate(incident_freqs):
        raman_intensities = raman_intensities_array[:, i]
        lorentzian_lineshapes = np.zeros((len(central_frequencies), len(frequencies)))
        gaussian_lineshapes = np.zeros((len(central_frequencies), len(frequencies)))

        for i, (central_frequency, intensity) in enumerate(zip(central_frequencies, raman_intensities)):
            lorentz = (1 / (1 + ((frequencies - central_frequency) / lorentzian_width)**2))
            lorentz /= np.trapezoid(lorentz, frequencies)
            lorentzian_lineshapes[i, :] = lorentz * intensity

            gauss = (np.exp(-0.5 * ((frequencies - central_frequency) / gaussian_width)**2) /
                    (gaussian_width * np.sqrt(2 * np.pi)))
            gauss /= np.trapezoid(gauss, frequencies)
            gaussian_lineshapes[i, :] = gauss * intensity

        gaussian_spectrum = np.sum(gaussian_lineshapes, axis=0)
        lorentzian_spectrum = np.sum(lorentzian_lineshapes, axis=0)
        output = (1 - lorentzian_percentage) * gaussian_spectrum + lorentzian_percentage * lorentzian_spectrum
        lorentzian_list.append(lorentzian_spectrum)
        gaussian_list.append(gaussian_spectrum)
        output_list.append(output)
        plt.figure(figsize=(5, 4))
        plt.plot(frequencies, lorentzian_spectrum, "r-", label="Lorentzian")
        plt.plot(frequencies, gaussian_spectrum, "b-", label="Gaussian")
        plt.plot(frequencies, output, "g-", label="Combined")
        plt.xlabel("Frequency")
        plt.ylabel("Raman Intensity")
        plt.title(f"Raman Spectrum: {freq} cm^-1 Incident Light")
        plt.legend()
        plt.grid(True)
        plt.show()
    np.savetxt(filedir / (filestem + f"_nrrspec_lorentzian{lorentzian_width}.csv"), np.column_stack([frequencies]+lorentzian_list), delimiter=',', header=header, comments='')
    np.savetxt(filedir / (filestem + f"_nrrspec_gaussian{gaussian_width}.csv"), np.column_stack([frequencies]+ gaussian_list), delimiter=',', header=header, comments='')
    return frequencies, gaussian_list, lorentzian_list, output_list
# Example usage:
# nrr,inc_freqs = extract_nrr(r"C:\Users\kentc\OneDrive - University of Rochester\Documents\Lab\gaussian\bluehive\ys7\ZnP_Acc\ZnP_Acc_rr_s1_wb97xd_deftzvp.log")
# print(nrr,inc_freqs)  # columns: frequency, RamAct Fr=1, Fr=2, Fr=3

