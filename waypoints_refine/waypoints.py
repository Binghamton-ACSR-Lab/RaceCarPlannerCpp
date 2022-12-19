import numpy as np
import matplotlib.pyplot as plt

wheel_base = 0.35


def model(init_state, control):
    return np.array([control[1]*np.cos(control[0]+init_state[2]),
                     control[1]*np.sin(control[0]+init_state[2]),
                     control[1]*np.tan(control[0])/wheel_base])


def generate_path():

    t = np.array([2.5, 0.5, 1.2, 0.7, 0.8])
    init_state = np.array([5.0, 0.0, np.pi/2])
    states = [init_state]
    steps = np.array(t/0.1, dtype=int)

    controls = np.array([[-0.4, 0.5], [-0.1, 1.0], [-0.2, 1.0], [0.3, 1.0], [-0.2, 1.0]])

    for i in range(0, 5):
        for j in range(0, steps[i]):
            init_state = init_state+0.1*model(init_state, controls[i])
            states.append(np.copy(init_state))

    return states


state = generate_path()
np_state = np.array(state)
plt.plot(np_state[:, 0], np_state[:, 1])
plt.show()


from elastica import *

class WaypointsRefineSimulator(
    BaseSystemCollection, Constraints, Forcing, CallBacks
):
    pass

waypoints_sim = WaypointsRefineSimulator()
final_time = 50.0
dt = 0.01


# setting up test params
n_elem = len(state)-1
start = np.copy(state[0])
# start[2]=0

print("start:")
print(start)

phi0 = start[2]
direction = np.array([np.cos(phi0), np.sin(phi0), 0.0])
normal = np.array([0.0, 0.0, 1.0])

start[2]=0


print("direction:")
print(direction)
print("normal:")
print(normal)

base_length = 0.1

base_radius = 0.025
base_area = np.pi * base_radius ** 2
density = 1000
youngs_modulus = 1e3
# For shear modulus of 1e4, nu is 99!
poisson_ratio = 0.5
shear_modulus = youngs_modulus / (poisson_ratio + 1.0)

positions = np_state.transpose()
directions = np.copy(positions)
positions[2, :] = 0

print(positions)

waypoints_rod = CosseratRod.straight_rod(
    n_elem,
    start,
    direction,
    normal,
    base_length,
    base_radius,
    density,
    0.0,  # internal damping constant, deprecated in v0.3.0
    youngs_modulus,
    shear_modulus=shear_modulus,
    position=positions,
)

waypoints_sim.append(waypoints_rod)
waypoints_sim.constrain(waypoints_rod).using(
    OneEndFixedBC, constrained_position_idx=(0,), constrained_director_idx=(0,))
waypoints_sim.constrain(waypoints_rod).using(
    GeneralConstraint,constrained_position_idx=(-1,),translational_constraint_selector=np.array([True, True, True]),
)




# end_force_x = 1.0
# end_force = np.array([end_force_x, 0.0, 0.0])
# stretch_sim.add_forcing_to(stretchable_rod).using(
#     EndpointForces, 0.0 * end_force, end_force, ramp_up_time=1e-2
# )

# # add damping
# dl = base_length / n_elem
# # old damping model (deprecated in v0.3.0) values
# # dt = 0.01 * dl
# # damping_constant = 1.0
# dt = 0.1 * dl
# damping_constant = 0.1
# stretch_sim.dampen(stretchable_rod).using(
#     AnalyticalLinearDamper,
#     damping_constant=damping_constant,
#     time_step=dt,
# )

# Add call backs
class WaypointsRefineCallBack(CallBackBaseClass):
    """
    Tracks the velocity norms of the rod
    """

    def __init__(self, step_skip: int, callback_params: dict):
        CallBackBaseClass.__init__(self)
        self.every = step_skip
        self.callback_params = callback_params

    def make_callback(self, system, time, current_step: int):

        if current_step % self.every == 0:

            self.callback_params["time"].append(time)
            # Collect only x
            self.callback_params["position"].append(
                system.position_collection.copy()
            )
            return


recorded_history = defaultdict(list)
recorded_history["position"].append(waypoints_rod.position_collection.copy())

total_steps = int(final_time / dt)
waypoints_sim.collect_diagnostics(waypoints_rod).using(
    WaypointsRefineCallBack, step_skip=int(total_steps/7), callback_params=recorded_history
)

waypoints_sim.finalize()
timestepper = PositionVerlet()
# timestepper = PEFRL()


print("Total steps", total_steps)
integrate(timestepper, waypoints_sim, final_time, total_steps)

print(recorded_history["position"])





for index,position in enumerate(recorded_history["position"]):
    plt.plot(position[0,:],position[1,:],label="{}".format(index))
plt.legend()
plt.show()


if PLOT_FIGURE:
    # First-order theory with base-length
    expected_tip_disp = end_force_x * base_length / base_area / youngs_modulus
    # First-order theory with modified-length, gives better estimates
    expected_tip_disp_improved = (
            end_force_x * base_length / (base_area * youngs_modulus - end_force_x)
    )

    fig = plt.figure(figsize=(10, 8), frameon=True, dpi=150)
    ax = fig.add_subplot(111)
    ax.plot(recorded_history["time"], recorded_history["position"], lw=2.0)
    ax.hlines(base_length + expected_tip_disp, 0.0, final_time, "k", "dashdot", lw=1.0)
    ax.hlines(
        base_length + expected_tip_disp_improved, 0.0, final_time, "k", "dashed", lw=2.0
    )
    if SAVE_FIGURE:
        fig.savefig("axial_stretching.pdf")
    plt.show()

if SAVE_RESULTS:
    import pickle

    filename = "axial_stretching_data.dat"
    file = open(filename, "wb")
    pickle.dump(stretchable_rod, file)
    file.close()

    tv = (
        np.asarray(recorded_history["time"]),
        np.asarray(recorded_history["velocity_norms"]),
    )

    def as_time_series(v):
        return v.T

    np.savetxt(
        "velocity_norms.csv",
        as_time_series(np.stack(tv)),
        delimiter=",",
    )
