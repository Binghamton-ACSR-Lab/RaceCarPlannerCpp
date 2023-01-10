import math

import numpy as np
import matplotlib.pyplot as plt
from acsr_geometry import ObstacleForce
import shapely

# with open('single.csv') as f:
#     lines=f.readlines()
#     for line in lines:
#         myarray = np.fromstring(line, dtype=float, sep=',')
#         print(myarray)
#
# def model(init_state, control):
#     wheel_base = 0.35
#     return np.array([control[1]*np.cos(control[0]+init_state[2]),
#                      control[1]*np.sin(control[0]+init_state[2]),
#                      control[1]*np.tan(control[0])/wheel_base])


# def generate_path():
#
#     t = np.array([2.5, 0.5, 1.2, 0.7, 0.8])
#     init_state = np.array([5.0, 0.0, np.pi/2])
#     states = [init_state]
#     steps = np.array(t/0.1, dtype=int)
#
#     controls = np.array([[-0.4, 0.5], [-0.1, 1.0], [-0.2, 1.0], [0.3, 1.0], [-0.2, 1.0]])
#
#     for i in range(0, 5):
#         for j in range(0, steps[i]):
#             init_state = init_state+0.1*model(init_state, controls[i])
#             states.append(np.copy(init_state))
#
#     return states
#
#
# state = generate_path()
# np_state = np.array(state)
# plt.plot(np_state[:, 0], np_state[:, 1])
# plt.show()


from elastica import *

class WaypointsRefineSimulator(
    BaseSystemCollection, Constraints, Forcing, Damping,CallBacks
):
    pass

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


def process(state,obstacls,plot = False):
    n_elem = len(state)-1
    positions = np.zeros((3,n_elem+1))
    positions[0:2,:] = state.T[0:2,:]


    if state.shape[1]>2:
        phi0 = state[0,2]
    else:
        phi0 =math.atan2(state[1,1]-state[0,1],state[1,0]-state[0,0])

    start=np.zeros((3,))
    start[0:2]=state[0]

    direction = np.array([np.cos(phi0), np.sin(phi0), 0.0])
    normal = np.array([0.0, 0.0, 1.0])
    final_time = 1
    dt = 0.01

    base_length = 0.1
    base_radius = 0.2
    base_area = np.pi * base_radius ** 2
    density = 1e5
    youngs_modulus = 1e6
    # For shear modulus of 1e4, nu is 99!
    poisson_ratio = 0.03
    shear_modulus = youngs_modulus / (poisson_ratio + 1.0)


    plot_history=[np.copy(positions)]

    total_its = 5
    for i in range(0,total_its):
        lengths = np.linalg.norm(positions[0:2,1:]-positions[0:2,0:-1],axis=0)
        rest_lengths = lengths*0.2
        waypoints_sim = WaypointsRefineSimulator()
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
            # rest_lengths = rest_lengths,
        )
        waypoints_rod.rest_lengths = rest_lengths

        waypoints_sim.append(waypoints_rod)
        waypoints_sim.constrain(waypoints_rod).using(
            FixedConstraint,
            constrained_position_idx=(0,-1),
            constrained_director_idx=(0,),

        )
        waypoints_sim.dampen(waypoints_rod).using(
            AnalyticalLinearDamper,
            damping_constant=1.0,
            time_step=dt,
        )
        waypoints_sim.constrain(waypoints_rod).using(
            GeneralConstraint,
            constrained_position_idx=(1,),
            translational_constraint_selector=np.array([False, True, True]),
        )
        waypoints_sim.add_forcing_to(waypoints_rod).using(
            ObstacleForce,
            data=obstacls
        )
        # waypoints_sim.constrain(waypoints_rod).using(
        #     GeneralConstraint,
        #     constrained_position_idx=(-1,),
        #     translational_constraint_selector=np.array([True, True, True]),
        # )
        # for i in range(1,n_elem-2):
        #     waypoints_sim.constrain(waypoints_rod).using(
        #         GeneralConstraint,
        #         constrained_position_idx=(i,),
        #         # constrained_director_idx=(i,),
        #         translational_constraint_selector=np.array([False, False, True]),
        #         # rotational_constraint_selector=np.array([True, True, False]),
        #     )

        recorded_history = defaultdict(list)
        recorded_history["position"].append(waypoints_rod.position_collection.copy())

        total_steps = 1000
        waypoints_sim.collect_diagnostics(waypoints_rod).using(
            WaypointsRefineCallBack, step_skip=1, callback_params=recorded_history
        )

        waypoints_sim.finalize()
        timestepper = PositionVerlet()
        # timestepper = PEFRL()
        integrate(timestepper, waypoints_sim, 0.5, total_steps,progress_bar=False)
        positions = recorded_history["position"][-1]

        if plot:
            plot_history.append(positions)

    if plot:
        shapes=[]
        for d in obstacls:
            d = np.reshape(d,(-1,2))
            shapes.append(shapely.Polygon(d))
        for poly in shapes:
            x,y = poly.exterior.xy
            plt.plot(x,y,'-g')

        for index,position in enumerate(plot_history):
            plt.plot(position[0,:],position[1,:],label="{}".format(index))
        plt.legend()
        plt.show()

    return positions

if __name__=='__main__':

    obstacles=[]
    with open('../data/map/refined_obstacle.txt') as f:
        lines=f.readlines()
    for line in lines:
        myarray = np.fromstring(line, dtype=float, sep=',')
        obstacles.append(myarray)

    states = np.loadtxt("../data/map/sst_data.txt")
    process(states,obstacles,True)


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










