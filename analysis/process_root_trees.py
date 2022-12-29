"""
A script to process the ROOT tree outputs of the Geant4 simulations of muons in ice.
The output is a csv file containing track lengths of the secondaries below a user-defined energy threshold
for each of the processes in which the secondaries were created.

Author: Tania Kozynets
tetiana.kozynets@nbi.ku.dk
December 29, 2022

"""
import sys

sys.path.append("../")

from argparse import ArgumentParser
from g4mu.utils.sorting import *
from particle import Particle
from numpy_indexed import group_by
from scipy.interpolate import interp1d
import pandas as pd
import numpy as np
from ROOT import TFile
import ROOT
from tqdm import tqdm
import glob
import os


# GLOBALS
# Refractive index (assume wavelength-independent)
n_refr = 1.33


def get_event_data(ttree, event_num):

    """
    Extract Geant4 step-level information from the ROOT trees for event with number "event_num".
    Input:
    ttree: input ROOT TTree object
    event_num: number of the event to process (int)

    Returns:
    secondary_data_dict: dictionary containing the secondary data for each propagation step
    """

    # Get the tree data for the current event
    ttree.GetEntry(event_num)
    # PDG IDs of the track-generating particles (e.g. muon, electron...)
    sec_pids = np.array(ttree.pid).astype(int)
    # Unique Geant4 track IDs
    track_ids = np.array(ttree.tid).astype(int)
    # Track IDs of the track parents
    parent_track_ids = np.array(ttree.parid).astype(int)

    # List of all processes (str)
    process_names_list = list(ttree.proc)
    process_names = np.array(
        [str(process_names_list[j]) for j in range(len(process_names_list))]
    )

    # Track lengths of the secondary particles
    secondary_track_lengths = np.array(ttree.track_length).astype(float)
    # Kinetic energies of the secondary particles
    secondary_ekin = np.array(ttree.ekin).astype(float)

    secondary_data_dict = {
        "secondary_pid": sec_pids,
        "track_id": track_ids,
        "parent_track_id": parent_track_ids,
        "parent_process": process_names,
        "secondary_ekin": secondary_ekin,
        "secondary_track_length": secondary_track_lengths,
    }

    return secondary_data_dict


def group_data_by_track_id(ttree, secondary_data_dict):

    """
    Extract Geant4 track-level information from the step-level information recorded as dict (see get_event_data above).
    Input:
    secondary_data_dict: the output of get_event_data for the event of interest (dict)
    ttree: input ROOT TTree object

    Returns:
    grouped_secondary_data: dictionary containing the secondary data grouped into tracks
    muon_data: dictionary containing the primary muon data (energy in GeV and track length in m)
    """

    # Each track contains multiple propagation steps. We need to group the tree data by tracks using track IDs:
    grouping = group_by(secondary_data_dict["track_id"])
    # The PID of the track-generating particle
    grouped_pids = grouping.first(secondary_data_dict["secondary_pid"])[1]
    # The track ID of the current track's parent
    grouped_parent_track_ids = grouping.first(secondary_data_dict["parent_track_id"])[1]
    # The track ID of the current track
    grouped_track_ids = grouping.first(secondary_data_dict["track_id"])[1]
    # The process in which the track-generating particle was originally created
    grouped_processes = grouping.first(secondary_data_dict["parent_process"])[1]
    # The initial kinetic energy of the secondary is the maximum energy of the track
    # (as the secondary loses energy over the propagation steps stored in the track)
    grouped_ekin = grouping.max(secondary_data_dict["secondary_ekin"])[1]
    # The total track length is the max of the track lengths stored in propagation steps
    # (what is stored in each step is the cumulative track length value, so we just need the final one)
    grouped_track_lengths = grouping.max(secondary_data_dict["secondary_track_length"])[
        1
    ]

    # Total muon track length in meters
    muon_track_length_m = float(grouped_track_lengths[grouped_pids == 13][0]) / 1e3
    # Total muon kinetic energy in GeV
    muon_energy_gev = np.array(ttree.ekini)[0]

    grouped_secondary_data = {
        "secondary_pid": grouped_pids,
        "track_id": grouped_track_ids,
        "parent_track_id": grouped_parent_track_ids,
        "parent_process": grouped_processes,
        "secondary_ekin": grouped_ekin,
        "secondary_track_length": grouped_track_lengths,
    }

    muon_data = {"energy": muon_energy_gev, "track_length": muon_track_length_m}

    return grouped_secondary_data, muon_data


def get_grandparent_processes(track_ids, parent_track_ids, parent_processes):
    """
    Since it is possible that some electrons were created not as the result of direct muon ionization,
    but as the result of the secondary electron ionization ('eIoni' process), we need to search over the "grantparent"
    processes in addition to the direct parent processes.

    Inputs:
    - track_ids: Geant4 track IDs of the secondary tracks (np.array, int)
    - parent_ids: Geant4 track IDs of the parent tracks for each track in track_ids (np.array, int)
    - parent_processes: names of the processes that created each track in track_ids (np.array, str)

    Returns:
    np.array of grandparent process names (str)
    """

    # To get the grandparent's process, create a lookup table over the track IDs (this is MUCH faster than np.where search!)
    track_dict = {track_ids[x]: x for x in range(len(track_ids))}
    # The first parent is a muon and was not created in any processes, so we just append a 'None' at the beginning
    grandparent_processes = np.append(
        "None",
        parent_processes[
            np.array([track_dict[parent] for parent in parent_track_ids[1:]])
        ],
    )

    return grandparent_processes


def get_secondary_track_lengths_for_process(
    process, grouped_secondary_data, energy_cut
):

    """
    Extract sums of the secondary track lengths with and without the Frank-Tamm correction
    (see Eq. 6 in https://arxiv.org/pdf/1206.5530.pdf).

    Inputs:
    - process: name of the process in which the secondaries where created; one of "muIoni", "muPairProd", "muBrems", or "muonNuclear".
    This includes both direct parent and "grandparent" processes.
    - grouped_secondary_data: secondary data extracted from the ROOT trees and grouped into tracks as dict (first output of group_data_by_track_id)
    - energy_cut: upper energy threshold (in GeV) of the secondary particles to process (float)

    Returns:
    - total_sec_track_length_m: sum of the secondary track lengths without the Frank-Tamm correction in meters (float)
    - total_sec_track_length_with_tamm_m: sum of the secondary track lengths with the Frank-Tamm correction in meters (float)
    - sec_ekin_sum: sum of the secondary energies below energy_cut in GeV (float)

    """

    grandparent_processes = get_grandparent_processes(
        grouped_secondary_data["track_id"],
        grouped_secondary_data["parent_track_id"],
        grouped_secondary_data["parent_process"],
    )
    # We search for an electron/positron in all processes but the muon nuclear interaction
    if process != "muonNuclear":

        mask = (
            (np.abs(grouped_secondary_data["secondary_pid"]) == 11)
            & (grouped_secondary_data["secondary_ekin"] < energy_cut)
            & (
                (grandparent_processes == process)
                | (grouped_secondary_data["parent_process"] == process)
            )
        )
    # Otherwise, we search for any charged particles (hadrons)
    # TODO: check the actual charge instead of just removing the photons
    else:
        mask = (
            (grouped_secondary_data["secondary_ekin"] < energy_cut)
            & (
                (grandparent_processes == process)
                | (grouped_secondary_data["parent_process"] == process)
            )
            & (grouped_secondary_data["secondary_pid"] != 22)
        )

    if np.sum(mask) != 0:

        # The PID of the secondary
        sec_pid = grouped_secondary_data["secondary_pid"][mask][0]
        # The mass of the secondary
        sec_mass = Particle.from_pdgid(sec_pid).mass / 1e3
        # The kinetic energy of the secondary
        sec_ekin = grouped_secondary_data["secondary_ekin"][mask]
        # The total energy of the secondary
        sec_etot = sec_ekin + sec_mass
        # Lorentz gamma
        sec_gamma = sec_etot / sec_mass
        # Lorentz beta
        sec_beta = np.sqrt(1 - (1 / sec_gamma) ** 2)
        # Frank-Tamm correction factors as per https://arxiv.org/pdf/1206.5530.pdf
        frank_tamm_factors = (1 - (1 / (n_refr * sec_beta) ** 2)) / (
            1 - (1 / n_refr**2)
        )
        frank_tamm_factors = np.where(frank_tamm_factors < 0.0, 0.0, frank_tamm_factors)

        # Sum of the secondary track lengths without the Frank-Tamm correction
        total_sec_track_length_m = (
            np.sum(grouped_secondary_data["secondary_track_length"][mask]) / 1e3
        )
        # Sum of the secondary track lengths with the Frank-Tamm correction
        total_sec_track_length_with_tamm_m = (
            np.sum(
                grouped_secondary_data["secondary_track_length"][mask]
                * frank_tamm_factors
            )
            / 1e3
        )

        sec_ekin_sum = np.sum(sec_ekin)

    else:
        return None

    return total_sec_track_length_m, total_sec_track_length_with_tamm_m, sec_ekin_sum


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument(
        "-id", "--input-dir", type=str, help="Input root tree directory", required=True
    )
    parser.add_argument(
        "-od", "--output-dir", type=str, help="Output csv file directory", required=True
    )
    parser.add_argument(
        "-ecut",
        "--energy-cut",
        type=float,
        help="Upper energy limit for the secondaries to process (GeV)",
        default=0.5,
    )
    args = parser.parse_args()

    # Process charged secondaries below this energy
    energy_cut = args.energy_cut

    # Processes in which the secondaries were generated
    extract_processes = ["muIoni", "muPairProd", "muonNuclear", "muBrems"]

    # Placeholder for the secondary particle records
    secondary_data = {}
    for process in extract_processes:
        secondary_data[process] = []

    # Get the file names of the Geant4 sim output trees.
    input_files = glob.glob(os.path.join(args.input_dir, "*.root"))
    input_files.sort(key=natural_keys)

    print("Found %s input ROOT file(s)" % len(input_files))

    for file_path in input_files:
        print("Processing file %s" % file_path)
        # Get the ROOT tree and its branches
        inp = TFile(file_path)
        ttree = inp.Get("muons_in_ice")

        # Get the number of events in the current file
        n_entries = ttree.GetEntries()
        print("Found %s events" % n_entries)

        # Loop through the events in the current file

        for i in tqdm(range(n_entries)):

            secondary_data_dict = get_event_data(ttree, event_num=i)

            grouped_secondary_data, muon_data = group_data_by_track_id(
                ttree, secondary_data_dict
            )
            for process in extract_processes:

                try:
                    (
                        total_sec_track_length_m,
                        total_sec_track_length_with_tamm_m,
                        sec_ekin_sum,
                    ) = get_secondary_track_lengths_for_process(
                        process, grouped_secondary_data, args.energy_cut
                    )

                    secondary_data[process].append(
                        [
                            total_sec_track_length_m,
                            total_sec_track_length_with_tamm_m,
                            sec_ekin_sum,
                            np.round(muon_data["energy"], 6),
                            muon_data["track_length"],
                        ]
                    )
                # If the get_secondary_track_lengths_for_process returned None (no secondaries passed the selection criteria),
                # skip this event
                except TypeError:
                    continue

    for process in extract_processes:

        print("Saving the csv output for process %s" % process)

        secondary_data_as_df = pd.DataFrame.from_records(
            np.array(secondary_data[process]),
            columns=[
                "sum_of_sec_track_lengths_without_FT_correction [m]",
                "sum_of_sec_track_lengths_with_FT_correction [m]",
                "sum_of_sec_kin_energies [GeV]",
                "muon_kinetic_energy [GeV]",
                "muon_track_length [m]",
            ],
        )

        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)

        secondary_data_as_df.to_csv(
            os.path.join(args.output_dir, "secondaries_%s.csv" % process),
            index=False,
        )
