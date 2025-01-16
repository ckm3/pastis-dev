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
argparser.add_argument("--total_batches", type=int, default=1)
argparser.add_argument("--scenario", type=str, default="BEB") #PLA, EB, BEB, BTP, PIB, TRIPLE
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
SCENARIO = args.scenario


def get_flat_attributes(obj, parent_key='', sep='.'):
    attributes = {}

    if isinstance(obj, list):
        for i, item in enumerate(obj):
            new_key = parent_key if len(obj) == 1 else f"{parent_key}[{i}]"
            attributes.update(get_flat_attributes(item, new_key, sep=sep))
        return attributes

    items = obj.items() if  isinstance(obj, dict) else obj.__dict__.items()

    for key, value in items:
        if isinstance(value, list):
            for i, item in enumerate(value):
                new_key = f"{parent_key}{sep}{key}[{i}]" if parent_key else f"{key}[{i}]"
                if hasattr(item, '__dict__'):
                    attributes.update(get_flat_attributes(item, new_key, sep=sep))
                else:
                    attributes[new_key] = item
        elif not key.startswith('_') and not callable(value):
            new_key = f"{parent_key}{sep}{key}" if parent_key else key
            if hasattr(value, '__dict__') or isinstance(value, dict):
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
    object_list, rejection_list = s.build_objects(input_dict, flag, True, verbose=False)

    # Compute model light curves
    lc = s.lightcurves(object_list, scenario=SCENARIO, lc_cadence_min=2.0)

    if not os.path.exists(f"./simulations/{SCENARIO}/" ):
        os.makedirs(f"./simulations/{SCENARIO}/", exist_ok=True)

    output_df = pd.DataFrame()
    for simu_number in range(len(lc)):
        attributes_dict = get_flat_attributes(object_list[simu_number])
        for key in list(attributes_dict.keys()):
            if 'drift' in key.lower():
                del attributes_dict[key]
            if "ticid" in key.lower():
                attributes_dict["TIC"] = attributes_dict[key]
                del attributes_dict[key]

        df = pd.DataFrame([attributes_dict])
        output_df = pd.concat([output_df, df])

        # save simulations
        simu_name = f"./simulations/{SCENARIO}/lcs/{SCENARIO}-simu-{part_num}-{simu_number}.csv"
        if not os.path.exists(simu_name):
            os.makedirs(os.path.dirname(simu_name), exist_ok=True)
        print("Saving slice:", part_num, "simulation:", simu_number)
        np.savetxt(simu_name, lc[simu_number][0], delimiter=",")  # as np array
    
    if output_df.empty:
        return
    output_df = pd.merge(output_df, pd_tess, on="TIC", how="inner")
    output_df.dropna(axis=1, how='all', inplace=True)
    output_df.to_csv(f"./simulations/{SCENARIO}/{SCENARIO}-parameters-{part_num}.csv", index=False)

    rej_df = pd.DataFrame()
    for rej in rejection_list:
        content_dict = get_flat_attributes(rej)
        for key in list(content_dict.keys()):
            if "ticid" in key.lower():
                content_dict["TIC"] = int(content_dict[key])
                del content_dict[key]
        df = pd.DataFrame([content_dict])
        rej_df = pd.concat([rej_df, df])

    if rej_df.empty:
        return
    rej_df = pd.merge(rej_df, pd_tess, on="TIC", how="inner")
    rej_df.to_csv(f"./simulations/{SCENARIO}/{SCENARIO}-rejections-{part_num}.csv", index=False)



print("Reading input files")

filenames = [
    "filled_spoc_gaia_3-10k_ruwe105.csv"
]

# filenames = filenames * 5  # quick and dirty way to repeat stars, I love it

full_data = pd.DataFrame([])
full_data_PD = pd.DataFrame([])

for file in filenames:
    print("Reading:", file)

    # read files
    data_pd = pd.read_csv(file).dropna().sample(frac=1, random_state=RANDOM_SEED)
    # remove ticid in konwn tfop
    data_pd = data_pd[~data_pd.TIC.isin(np.genfromtxt("known_tfop.txt"))]

    params_pd = data_pd[["Rad", "Tmag", "Av", "mass", "Teff", "logg", "MH", "Gmag", "BP-RP", "B","TIC","distance"]].copy()

    full_data = pd.concat([full_data, params_pd])
    full_data_PD = pd.concat([full_data_PD, data_pd])

def process_batch(start, end, part, full_data, full_data_PD):
    print(start, end, "Part:", part)
    params = full_data.iloc[np.arange(start, end)%len(full_data)].values.T
    gen_files(params, part, full_data_PD, method="hsu")

if __name__ == "__main__":
    # Split into batches and process
    batch_size = 100000 # send this number of objects to each process
    start = args.batch_id * batch_size
    num_batches = 48  # number of simulations
    num_batches = min(num_batches, args.total_batches - args.batch_id)
    print(num_batches, args.batch_id, args.total_batches)
    # num_batches = (len(full_data_PD) + batch_size - 1) // batch_size  # Ceiling division

    with mp.Pool(num_batches) as pool:
        results = []
        for part in range(num_batches):
            end = start + batch_size
            results.append(pool.apply_async(process_batch, (start, end, args.batch_id + part, full_data, full_data_PD)))
            start = end

        for result in results:
            result.get()  # Ensure all processes complete

    print("All batches processed.")
