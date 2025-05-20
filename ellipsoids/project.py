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
            default="shortgpu-v100",
            help="Specify the partition to submit to."
        )


class Fry(DefaultSlurmEnvironment):
    hostname_pattern = "fry"
    template = "fry.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="v100," "batch",
            help="Specify the partition to submit to."
        )


@Ellipsoids.label
def system_built(job):
    return job.isfile("init_frame.gsd")


@Ellipsoids.label
def initial_run_done(job):
    return job.doc.runs > 0


@Ellipsoids.label
def equilibrated(job):
    return job.doc.equilibrated

@Ellipsoids.label
def production_done(job):
    return job.isfile("production-restart.gsd")

@Ellipsoids.post(initial_run_done)
@Ellipsoids.operation(
    directives={"ngpu": 1, "ncpu": 1, "executable": "python -u"}, name="run"
)
def run(job):
    """Run initial single-chain simulation."""
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
        job.doc.n_particles = int(job.doc.num_mols * job.doc.lengths)
        ellipsoid_chain = EllipsoidChain(num_mols=job.doc.num_mols, lengths=job.doc.lengths,lpar=1.0,bead_mass=1.0)
        system = Pack(molecules=ellipsoid_chain, density=job.sp.density*Unit("nm**-3"), 
                      packing_expand_factor=11,edge=2,overlap=1,fix_orientation=True)
        system.to_gsd(job.fn("init_frame.gsd"))
        print("Finished.")

        # Set up Simulation obj
        gsd_path = job.fn(f"trajectory{job.doc.runs}.gsd")
        log_path = job.fn(f"log{job.doc.runs}.txt")

        ff = EllipsoidForcefield(epsilon=1.0,lpar=1.0,lperp=0.5,r_cut=2.0,bond_k=100,bond_r0=0,angle_k=30,angle_theta0=1.9)
        ff.hoomd_forces
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
        target_box = get_target_box_number_density(density=job.sp.density*Unit("nm**-3"),n_beads=job.doc.num_mols * 
                                                   job.doc.lengths)
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
        sim.run_NVT(n_steps=job.sp.n_equil_steps, kT=job.sp.kT, tau_kt=tau_kT)
        sim.save_restart_gsd(job.fn("restart.gsd"))
        job.doc.runs = 1
        print("Simulation finished.")

@Ellipsoids.pre(initial_run_done)
@Ellipsoids.post(equilibrated)
@Ellipsoids.operation(
    directives={"ngpu": 1, "ncpu": 1, "executable": "python -u"},
    name="run-longer"
)
def run_longer(job):
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
        rigid_frame, rigid = create_rigid_ellipsoid_chain(system.hoomd_snapshot)
        with open(job.fn("forcefield.pickle"), "rb") as f:
            hoomd_ff = pickle.load(f)
        gsd_path = job.fn(f"trajectory{job.doc.runs}.gsd")
        log_path = job.fn(f"log{job.doc.runs}.txt")

        sim = Simulation(
            initial_state=job.fn("restart.gsd"),
            forcefield=hoomd_ff,
            constraint=rigid,
            dt=job.sp.dt,
            gsd_write_freq=job.sp.gsd_write_freq,
            gsd_file_name=gsd_path,
            log_write_freq=job.sp.log_write_freq,
            log_file_name=log_path,
            seed=job.sp.sim_seed,
        )
        print("Running simulation.")
        sim.run_NVT(
            n_steps=1e7,
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
        rigid_frame, rigid = create_rigid_ellipsoid_chain(system.hoomd_snapshot)
        with open(job.fn("forcefield.pickle"), "rb") as f:
            hoomd_ff = pickle.load(f)

        gsd_path = job.fn(f"production.gsd")
        log_path = job.fn(f"production.txt")

        sim = Simulation(
            initial_state=job.fn("restart.gsd"),
            forcefield=hoomd_ff,
            constraint=rigid,
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
   

@Ellipsoids.pre(production_done)
@Ellipsoids.operation(
    directives={"ngpu": 1, "ncpu": 1, "executable": "python -u"},
    name="production_run_longer"
)
def production_run_longer(job):
    import unyt
    from unyt import Unit
    import flowermd
    from flowermd.base import Simulation
    import hoomd
    with job:
        print("------------------------------------")
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        print("Restarting and continuing simulation...")
        print("Continuing the production run...")
        with open(job.fn("forcefield.pickle"), "rb") as f:
            hoomd_ff = pickle.load(f)

        gsd_path = job.fn(f"production{job.doc.production_runs+1}.gsd")
        log_path = job.fn(f"production{job.doc.production_runs+1}.txt")
        ref_values = get_ref_values(job)

        sim = Simulation(
            initial_state=job.fn("production-restart.gsd"),
            forcefield=hoomd_ff,
            reference_values=ref_values,
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
        print("Simulation finished.")
        job.doc.production_runs += 1
        
if __name__ == "__main__":
    Ellipsoids(environment=Fry).main()
