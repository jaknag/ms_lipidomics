import requests
import json
import re
from urllib import request
import gzip
import shutil
import requests
import pandas as pd
import time
import os
import numpy as np
import multiprocessing
import re

def downloadandprocess_exposure(ftp_path, reported_trait):
    # Download process
    print("Downloading: ", ftp_path)

    file_name = ftp_path.split('/')[-1]  # Get the file name (e.g., GCST90449363.tsv.gz)
    file_base_name = file_name.replace(".tsv.gz", "")  # Remove the .tsv.gz extension
    # Format the reported_trait for the file name (replace spaces with underscores)
    reported_trait_formatted = reported_trait.replace(" ", "_")
    output_file = f"filtered_{file_base_name}_{reported_trait_formatted}.tsv"

    if os.path.exists(output_file):
        print("This file already exists: ", output_file)

    else:
        start_time = time.time()
        request.urlretrieve(url=ftp_path, filename=file_name)
        print('Processing: ', output_file)

        # Define the chunk size
        chunk_size = 1000000
        # Initialize an empty list to collect filtered data
        filtered_chunks = []

        # Process the file in chunks
        for chunk in pd.read_csv(file_name, sep='\t', compression='gzip', chunksize=chunk_size):
            # Filter rows where p_val <= 5e-8
            filtered_chunk = chunk[chunk['neg_log_10_p_value'] > threshold]

            # Append the filtered chunk to the list
            filtered_chunks.append(filtered_chunk)

        # Concatenate all filtered chunks into a single DataFrame
        filtered_data = pd.concat(filtered_chunks, ignore_index=True)
        # Add the 'p_val' column after filtering
        filtered_data['p_val'] = 10 ** (-filtered_data['neg_log_10_p_value'])

        # Save the filtered data to a new file
        filtered_data.to_csv(output_file, sep='\t', index=False)
        print(f"Filtered data saved to: {output_file}")

        # Delete the original file
        os.remove(file_name)
        print(f"Deleted the original file: {file_name}")

        print("--- %s seconds ---" % (time.time() - start_time))


# Function to handle multiprocessing for downloading and processing all URLs
def downloadandprocess_all_exposures(ftp_paths, reported_traits):
    # Create a pool of worker processes
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        # Distribute the download tasks to worker processes
        pool.starmap(downloadandprocess_exposure, zip(ftp_paths, reported_traits))



ms_path = '/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended'
exposures_path = f"{ms_path}/exposures/VLDL-TG/"

os.chdir(ms_path)

ftp_list = pd.read_csv('ftp_list.csv')

# Filter rows where 'study_tag' starts with 'meta_full'
filtered_ftp_list = ftp_list[ftp_list['study_tag'].str.startswith('meta_EUR')]
ftp_paths = filtered_ftp_list['ftp_path'].tolist()
reported_traits = filtered_ftp_list['reported_trait'].tolist()

os.chdir(exposures_path)

# Set the threshold
# If threshold required is p<5e-8, then insert this value below
# Filtering in this data is through -log10 of the p value
threshold = -np.log10(1e-6)

downloadandprocess_all_exposures(ftp_paths, reported_traits)
