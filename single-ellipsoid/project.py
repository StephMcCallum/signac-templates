"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help
"""
import signac
import pickle
from flow import FlowProject, directives
from flow.environment import DefaultSlurmEnvironment
import os
from unyt import Unit

class Ellipsoids(FlowProject):
    pass


class Borah(DefaultSlurmEnvironment):
    hostname_pattern = "borah"
    template = "borah.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="gpu-v100",
            help="Specify the partition to submit to."
        )

@Ellipsoids.label
def system_built(job):
    return job.isfile("shrink_restart.gsd")


@Ellipsoids.label
def initial_run_done(job):
    return job.doc.runs > 0


@Ellipsoids.label
def equilibrated(job):
    return job.doc.equilibrated

@Ellipsoids.label
def production_done(job):
    return job.isfile("production-restart.gsd")

@Ellipsoids.label
def restart_rigid_ellipsoid(): #this function needs to updated to be more extensible
    local_coords = [(0.0, 0.0, 0.0), (1.049999999999999, 0.0, 0.0), (1.0, 0.0, 0.0), (-1.0000000000000009, 0.0, 0.0)]
    rigid_constrain = hoomd.md.constrain.Rigid()
    rigid_constrain.body["R"] = {
        "constituent_types": ["X", "A", "T", "T"],
        "positions": local_coords,
        "orientations": [[1, 0, 0, 0]] * len(local_coords),
    }
    return rigid_constrain

@Ellipsoids.post(system_built)
@Ellipsoids.operation(
    directives={"ngpu": 1, "ncpu": 1, "executable": "python -u"}, name="build"
)
def build(job):
    """Build ellipsoid system and run shrink simulation."""
    import unyt
    from unyt import Unit
    import flowermd
    from flowermd.base import Simulation,Pack
    from flowermd.library import EllipsoidForcefield, EllipsoidChain
    from flowermd.utils import get_target_box_number_density
    from flowermd.utils.constraints import create_rigid_ellipsoid_chain
    import hoomd
    with job:
        print("------------------------------------")
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        print("Building initial frame.")
        job.doc.n_particles = int(job.sp.N * job.sp.length)
        ellipsoid_chain = EllipsoidChain(num_mols=job.sp.N,
                                         lengths=job.sp.length,
                                         lpar=job.sp.lpar,
                                         bead_mass=job.sp.bead_mass
                                        )
        system = Pack(molecules=ellipsoid_chain,
                      density=job.sp.density*Unit("nm**-3"), 
                      packing_expand_factor=job.sp.packing_expand_factor,
                      edge=job.sp.edge,
                      overlap=job.sp.overlap,
                      fix_orientation=job.sp.fix_orientation,
                     )
        system.to_gsd(job.fn("init_frame.gsd"))
        print("Finished.")

        # Set up Simulation obj
        gsd_path = job.fn(f"trajectory{job.doc.runs}.gsd")
        log_path = job.fn(f"log{job.doc.runs}.txt")

        ff = EllipsoidForcefield(epsilon=job.sp.epsilon,
                                 lpar=job.sp.lpar,
                                 lperp=job.sp.lper,
                                 r_cut=job.sp.r_cut,
                                 bond_k=job.sp.bond_k,
                                 bond_r0=job.sp.bond_r0
                                )
        rigid_frame, rigid = create_rigid_ellipsoid_chain(system.hoomd_snapshot)
        
        sim = Simulation(
            initial_state=job.fn("init_frame.gsd"),
            forcefield=ff.hoomd_forces,
            constraint=rigid,
            dt=job.sp.dt,
            gsd_write_freq=job.sp.gsd_write_freq,
            gsd_file_name=gsd_path,
            log_write_freq=job.sp.log_write_freq,
            log_file_name=log_path,
            seed=job.sp.sim_seed,
        )
        sim.pickle_forcefield(job.fn("forcefield.pickle"))
        # Store more unit information in job doc
        tau_kT = job.sp.dt * job.sp.tau_kT
        job.doc.tau_kT = tau_kT
        job.doc.real_time_step = sim.real_timestep.to("fs").value
        job.doc.real_time_units = "fs"
        target_box = get_target_box_number_density(density=job.sp.density*Unit("nm**-3"),n_beads=job.sp.N * 
                                                   job.sp.length)
        job.doc.target_box = target_box.value
        shrink_kT_ramp = sim.temperature_ramp(
                n_steps=job.sp.n_shrink_steps,
                kT_start=job.sp.shrink_kT,
                kT_final=job.sp.kT
        )
        sim.run_update_volume(
                final_box_lengths=target_box,
                n_steps=job.sp.n_shrink_steps,
                period=job.sp.shrink_period,
                tau_kt=tau_kT,
                kT=shrink_kT_ramp
        )
        sim.save_restart_gsd(job.fn("shrink_restart.gsd"))
        print("Shrinking simulation finished...")

@Ellipsoids.pre(system_built)
@Ellipsoids.post(initial_run_done)
@Ellipsoids.operation(
    directives={"ngpu": 1, "ncpu": 1, "executable": "python -u"},
    name="run"
)
def run(job):
    import unyt
    from unyt import Unit
    import flowermd
    from flowermd.base import Simulation
    import hoomd
    import gsd.hoomd
    from flowermd.utils.constraints import create_rigid_ellipsoid_chain
    with job:
        print("------------------------------------")
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        print("Restarting and continuing simulation...")
        system = gsd.hoomd.open(job.fn("restart.gsd"),'r')
        with open(job.fn("forcefield.pickle"), "rb") as f:
            hoomd_ff = pickle.load(f)
        gsd_path = job.fn(f"trajectory{job.doc.runs}.gsd")
        log_path = job.fn(f"log{job.doc.runs}.txt")

        sim = Simulation(
            initial_state=job.fn("shrink_restart.gsd"),
            forcefield=hoomd_ff,
            constraint=restart_rigid_ellipsoid(),
            dt=job.sp.dt,
            gsd_write_freq=job.sp.gsd_write_freq,
            gsd_file_name=gsd_path,
            log_write_freq=job.sp.log_write_freq,
            log_file_name=log_path,
            seed=job.sp.sim_seed,
        )
        print("Running simulation.")
        sim.run_NVT(
            n_steps=job.sp.n_equil_steps,
            kT=job.sp.kT,
            tau_kt=job.doc.tau_kT,
        )
        sim.save_restart_gsd(job.fn("restart.gsd"))
        job.doc.runs += 1
        print("Simulation finished.")

@Ellipsoids.pre(equilibrated)
@Ellipsoids.post(production_done)
@Ellipsoids.operation(
    directives={"ngpu": 1, "ncpu": 1, "executable": "python -u"},
    name="production"
)
def production_run(job):
    import unyt
    from unyt import Unit
    import flowermd
    from flowermd.base import Simulation
    import hoomd
    import gsd.hoomd
    from flowermd.utils.constraints import create_rigid_ellipsoid_chain
    with job:
        print("------------------------------------")
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        print("Restarting and continuing simulation...")
        print("Running the production run...")
        system = gsd.hoomd.open(job.fn("restart.gsd"),'r')
        with open(job.fn("forcefield.pickle"), "rb") as f:
            hoomd_ff = pickle.load(f)

        gsd_path = job.fn(f"production.gsd")
        log_path = job.fn(f"production.txt")

        sim = Simulation(
            initial_state=job.fn("restart.gsd"),
            forcefield=hoomd_ff,
            constraint=restart_rigid_ellipsoid(),
            dt=job.sp.dt,
            gsd_write_freq=int(5e5),
            gsd_file_name=gsd_path,
            log_write_freq=job.sp.log_write_freq,
            log_file_name=log_path,
            seed=job.sp.sim_seed,
        )
        print("Running simulation.")
        sim.run_NVT(
            n_steps=job.sp.n_prod_steps*2,
            kT=job.sp.kT,
            tau_kt=job.doc.tau_kT,
        )
        sim.save_restart_gsd(job.fn("production-restart.gsd"))
        job.doc.production_runs += 1
        print("Simulation finished.")
