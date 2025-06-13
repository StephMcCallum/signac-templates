#!/usr/bin/env python
"""Initialize the project's data space.

Iterates over all defined state points and initializes
the associated job workspace directories.
The result of running this file is the creation of a signac workspace:
    - signac.rc file containing the project name
    - signac_statepoints.json summary for the entire workspace
    - workspace/ directory that contains a sub-directory of every individual statepoint
    - signac_statepoints.json within each individual statepoint sub-directory.

"""

import signac
import flow
import logging
from collections import OrderedDict
from itertools import product


def get_parameters():
    ''''''
    parameters = OrderedDict()
    #system
    parameters["N"] = [128]
    parameters["length"] = [1] #chain length for a polymer chain
    parameters["harmonic_bonds"] = [True]
    parameters["bead_mass"] = [1.0]
    #pack parameters
    parameters["density"] = [0.6,0.8,1.0]
    parameters["packing_expand_factor"] = [11]
    parameters["edge"] = [2]
    parameters["overlap"] = [1]
    parameters["fix_orientation"] = [True]
    #forcefield
    parameters["epsilon"] = [1.0]
    parameters["lpar"] = [1.0]
    parameters["lper"] = [0.5]
    parameters["r_cut"] = [2.5]
    parameters["bond_k"] = [100]
    parameters["bond_r0"] = [0]
    #simulation
    parameters["kT"] = [2.0]
    parameters["n_equil_steps"] = [1e6]
    parameters["n_prod_steps"] = [1e8]
    parameters["shrink_kT"] = [6.0]
    parameters["n_shrink_steps"] = [5e8]
    parameters["shrink_period"] = [10000]
    parameters["dt"] = [0.001]
    parameters["tau_kT"] = [10]
    parameters["gsd_write_freq"] = [1e6]
    parameters["log_write_freq"] = [1e4]
    parameters["sim_seed"] = [42]
    parameters["system_seed"] = [34]
    return list(parameters.keys()), list(product(*parameters.values()))


def main():
    project = signac.init_project()
    param_names, param_combinations = get_parameters()
    # Create workspace of jobs
    for params in param_combinations:
        statepoint = dict(zip(param_names, params))
        job = project.open_job(statepoint)
        job.init()
        job.doc.setdefault("node","p100")
        job.doc.setdefault("equilibrated", False)
        job.doc.setdefault("sampled", False)
        job.doc.setdefault("runs", 0)
        job.doc.setdefault("production_runs",0)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
