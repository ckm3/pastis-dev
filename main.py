#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 17:28:29 2021

@author: rodrigo, modified by Kaiming
"""
import numpy as np
import pandas as pd
import multiprocessing as mp
import argparse
import os

argparser = argparse.ArgumentParser()
argparser.add_argument("--batch_id", type=int, default=0)
args = argparser.parse_args()

# Import relevant modules from PASTIS
from pastis import limbdarkening

# Initialise if needed
if not hasattr(limbdarkening, "LDCs"):
    limbdarkening.initialize_limbdarkening(["TESS"])

import draw as d
import simulation as s

# Read parameters
from parameters import SCENARIO, NSIMU_PER_TIC_STAR, THETAMIN_DEG, RANDOM_SEED


def get_flat_attributes(obj, parent_key='', sep='.'):
    attributes = {}
    for key, value in obj.__dict__.items():
        if isinstance(value, list):
            for i, item in enumerate(value):
                new_key = f"{parent_key}{sep}{key}[{i}]" if parent_key else f"{key}[{i}]"
                if hasattr(item, '__dict__'):
                    attributes.update(get_flat_attributes(item, new_key, sep=sep))
                else:
                    attributes[new_key] = item
        elif not key.startswith('_') and not callable(value):
            new_key = f"{parent_key}{sep}{key}" if parent_key else key
            if hasattr(value, '__dict__'):
                attributes.update(get_flat_attributes(value, new_key, sep=sep))
            else:
                attributes[new_key] = value
    return attributes


def gen_files(params, part_num, pd_tess, **kwargs):
    # Draw parameters for scenario

    input_dict, flag = d.draw_parameters(
        params, SCENARIO, nsimu=NSIMU_PER_TIC_STAR, thetamin_deg=THETAMIN_DEG, **kwargs
    )

    # Create objects
    object_list, rej = s.build_objects(input_dict, np.sum(flag), True, verbose=False)

    # Compute model light curves
    lc = s.lightcurves(object_list, scenario=SCENARIO, lc_cadence_min=2.0)

    if not os.path.exists(f"./simulations/{SCENARIO}/" ):
        os.makedirs(f"./simulations/{SCENARIO}/", exist_ok=True)

    output_df = pd.DataFrame()
    for simu_number in range(len(lc)):
        attributes_dict = get_flat_attributes(object_list[simu_number][0])
        for key in list(attributes_dict.keys()):
            if 'drift' in key.lower():
                del attributes_dict[key]
            if "ticid" in key.lower():
                attributes_dict["TIC"] = attributes_dict[key]
                del attributes_dict[key]

        df = pd.DataFrame([attributes_dict])
        output_df = pd.concat([output_df, df])

        # save simulations
        simu_name = (
            f"./simulations/{SCENARIO}/lcs/"
            + SCENARIO
            + "-simu-"
            + str(part_num)
            + "-"
            + str(simu_number)
            + ".csv"
        )
        if not os.path.exists(simu_name):
            os.makedirs(os.path.dirname(simu_name), exist_ok=True)
        print("Saving slice:", part_num, "simulation:", simu_number)
        np.savetxt(simu_name, lc[simu_number][0], delimiter=",")  # as np array
    
    if output_df.empty:
        return
    output_df = pd.merge(output_df, pd_tess, on="TIC", how="inner")
    output_df.to_csv(f"./simulations/{SCENARIO}/" + SCENARIO + "-parameters-" + str(part_num) + ".csv", index=False)



print("Reading input files")

filenames = [
    "filled_spoc_gaia.csv_3-10k.csv"
]

# filenames = filenames * 5  # quick and dirty way to repeat stars, I love it

full_data = pd.DataFrame([])
full_data_PD = pd.DataFrame([])

for file in filenames:
    print("Reading:", file)

    # read files
    data_pd = pd.read_csv(file).sample(frac=1, random_state=RANDOM_SEED)

    params_pd = data_pd[["Rad", "Tmag", "Av", "mass", "Teff", "logg", "MH", "ebv", "B","TIC","distance"]].copy()

    full_data = pd.concat([full_data, params_pd])
    full_data_PD = pd.concat([full_data_PD, data_pd])

def process_batch(start, end, part, full_data, full_data_PD):
    print(start, end, "Part:", part)
    params = full_data.iloc[start:end].values.T
    gen_files(params, part, full_data_PD, method="uniform")

if __name__ == "__main__":
    # Split into batches and process
    batch_size = 1000 # send this number of objects to each process
    start = args.batch_id * batch_size
    num_batches = 1

    with mp.Pool(num_batches) as pool:
        results = []
        for part in range(num_batches):
            end = min(start + batch_size, len(full_data_PD))
            results.append(pool.apply_async(process_batch, (start, end, args.batch_id + part, full_data, full_data_PD)))
            start = end

        for result in results:
            result.get()  # Ensure all processes complete

    print("All batches processed.")
