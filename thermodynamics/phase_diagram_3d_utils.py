# importing required libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.figure_factory as ff

from plotly.offline import iplot
from typing import List, Callable
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull


# run command:
# python .\thermodynamics\phase_diagram_3d_utils.py

class Simplex:
    """Modified from https://github.com/materialsproject/pymatgen/blob/v2022.7.19/pymatgen/util/coord.py#L367
    Thank you, Pymatgen :)
    """

    def __init__(self, vertices: List[list]):
        self.vertices = np.array(vertices)
        self.dim = self.vertices.shape[0]
        self.augmented_matrix = np.concatenate([vertices, np.ones((self.dim, 1))], axis=-1)

    def bary_to_cart(self, bary_coords):
        return np.dot(bary_coords, self.augmented_matrix[:, :-1])

    @staticmethod
    def sample_between_bary_coords(b1, b2, n=31):
        """Sample n barycentric points (including endpoints) in between b1 and b2"""
        dc_vec = b2 - b1
        vec_increment = dc_vec / (n - 1)
        sample_points_initial = np.repeat([b1], n, axis=0)
        sample_points_addition = np.multiply(np.repeat([vec_increment], n, axis=0),
                                             np.arange(0, n)[:, np.newaxis])
        return sample_points_initial + sample_points_addition

    @staticmethod
    def get_point_between(b1, b2, frac):
        """Get a point in between b1 and b2 with its distance to b1 as frac * (b2-b1)"""
        dc_vec = b2 - b1
        return b1 + frac * dc_vec


class TernarySurface:
    def __init__(
            self,
            vertices,
            z_func: Callable,
            resolution=None,
            u_range=(0, 1),
            v_range=(0, 1),
    ):
        self.simplex = Simplex(vertices)
        self.z_func = z_func
        self.u_range = u_range
        self.v_range = v_range
        self.u_vals = None
        self.v_vals = None

        res_value = (100, 100)

        self.resolution = resolution if resolution is not None else res_value
        # self.delauney_tri = self.run_delauney()
        self.faces = None
        # self.get_surface()

    def _get_u_values_and_v_values(self):
        res = self.resolution
        if len(res) == 1:
            u_res = v_res = res[0]
        else:
            u_res, v_res = res

        u_values = np.linspace(*self.u_range, u_res + 1)
        v_values = np.linspace(*self.v_range, v_res + 1)

        # create a mesh based on u, v values
        uu, vv = np.meshgrid(u_values, v_values)
        uv_sum_values = uu + vv
        # the sum should not be greater than 1
        mask = (uv_sum_values <= 1)
        u_values = uu[mask]
        v_values = vv[mask]
        return u_values, v_values

    def get_cart_coord(self, u, v):
        w = 1 - u - v
        x, y, *_ = self.simplex.bary_to_cart([u, v, w])
        z = self.z_func(u, v)
        return [x, y, z]

    # def run_delauney(self):
    #     u_values, v_values = self._get_u_values_and_v_values()
    #     coords = [
    #         self.get_cart_coord(u, v)[:2]
    #         for u, v in zip(u_values, v_values)
    #     ]
    #     tri = Delaunay(
    #         points=np.array(coords)
    #     )
    #     return tri

    def get_surface(self):
        u_values, v_values = self._get_u_values_and_v_values()
        v_func = np.vectorize(self.z_func)
        z_values = v_func(u_values, v_values)
        # convert the u and v values from the barycentric coordinate to the cartesian coordinate
        uv_bary = np.stack((u_values, v_values, 1 - u_values - v_values), axis=1)
        uv_cart = self.simplex.bary_to_cart(uv_bary)
        u_cart = uv_cart[:, 0]
        v_cart = uv_cart[:, 1]
        return u_cart, v_cart, z_values


# define the axes attributes

# define the thermodynamic variables
temp = 400
a_ref_fcc = -1500
b_ref_fcc = -1000
c_ref_fcc = -2000
omega_ab_fcc = 5000
omega_ac_fcc = 8000
omega_bc_fcc = 7000
omega_abc_fcc = 10000

a_ref_liquid = 0
b_ref_liquid = 0
c_ref_liquid = 0
omega_ab_liquid = 5000
omega_ac_liquid = 3000
omega_bc_liquid = 2000
omega_abc_liquid = 10000


def get_entropy(x):
    if x <= 0:
        return 0
    return x * np.log(x)


def get_binary_contrib(x1, x2, l0, l1, l2):
    return x1 * x2 * (l0 + (x1 - x2) * l1 + (x1 - x2) ** 2 * l2)


def get_ternary_contrib(x1, x2, x3, l0, l1, l2, l3):
    return x1 * x2 * x3 * (
            l0 +
            (1 + 2 * x1 - x2 - x3) * l1 / 3 +
            (1 + 2 * x2 - x3 - x1) * l2 / 3 +
            (1 + 2 * x3 - x1 - x2) * l3 / 3
    )


def get_gibbs_fcc(x_a, x_b):
    x_c = 1 - x_a - x_b

    gibbs_ref = x_a * a_ref_fcc + x_b * b_ref_fcc + x_c * c_ref_fcc
    gibbs_ideal = 8.31 * temp * (get_entropy(x_a) + get_entropy(x_b) + get_entropy(x_c))
    gibbs_excess = (
            get_binary_contrib(x_a, x_b, omega_ab_fcc, -3000, -5000) +
            get_binary_contrib(x_a, x_c, omega_ac_fcc, -2000, -4000) +
            get_binary_contrib(x_b, x_c, omega_bc_fcc, -7000, -8000) +
            get_ternary_contrib(x_a, x_b, x_c,
                                omega_abc_fcc, 40000, -5000, -2000)
    )
    return gibbs_ref + gibbs_ideal + gibbs_excess


def get_gibbs_liquid(x_a, x_b):
    x_c = 1 - x_a - x_b

    gibbs_ref = x_a * a_ref_liquid + x_b * b_ref_liquid + x_c * c_ref_liquid
    gibbs_ideal = 8.31 * temp * (get_entropy(x_a) + get_entropy(x_b) + get_entropy(x_c))
    gibbs_excess = (
            x_a * x_b * omega_ab_liquid +
            x_a * x_c * omega_ac_liquid +
            x_b * x_c * omega_bc_liquid +
            x_a * x_b * x_c * omega_abc_liquid
    )
    return gibbs_ref + gibbs_ideal + gibbs_excess


default_vertices = [[0, 0, 0],
                    [1, 0, 0],
                    [1 / 2, np.sqrt(3) / 2, 0]]

gibbs_surface_fcc = TernarySurface(vertices=default_vertices,
                                   z_func=get_gibbs_fcc)
gibbs_surface_liquid = TernarySurface(vertices=default_vertices,
                                      z_func=get_gibbs_liquid)

# creating the dataset for plotting
xs_fcc, ys_fcc, zs_fcc = gibbs_surface_fcc.get_surface()
xs_liquid, ys_liquid, zs_liquid = gibbs_surface_liquid.get_surface()

# merge all the data into a single array
gibbs_fcc = np.vstack([xs_fcc, ys_fcc, zs_fcc]).T
gibbs_liquid = np.vstack([xs_liquid, ys_liquid, zs_liquid]).T
gibbs_combined = np.append(gibbs_fcc, gibbs_liquid, axis=0)

# get the minimum at each x,y coordinate
gibbs_min_df = pd.DataFrame(gibbs_combined, columns=["x", "y", "z"])
min_idx = gibbs_min_df.groupby(["x", "y"]).idxmin().to_numpy().flatten()
min_idx = np.sort(min_idx)

# obtain the convex hull
convex_hull = ConvexHull(points=gibbs_combined)
hull_vertices_idx = convex_hull.vertices
# get only the vertices that belong to the min bound
min_mask = np.in1d(hull_vertices_idx, min_idx)
hull_vertices_idx = hull_vertices_idx[min_mask]
hull_points = convex_hull.points[hull_vertices_idx]

hull_points_projected = np.copy(hull_points)
hull_points_projected[:, -1] = np.max(gibbs_combined) + 100

# # creating figure
# fig1 = plt.figure()
# ax = fig1.add_subplot(projection='3d', computed_zorder=True)
#
# # creating the plot
# plot_fcc = ax.plot_trisurf(xs_fcc, ys_fcc, zs_fcc, color='orange', alpha=0.5)
# plot_liquid = ax.plot_trisurf(xs_liquid, ys_liquid, zs_liquid, color='blue', alpha=0.5)
#
# # setting title and labels
# ax.set_xlabel('x-axis')
# ax.set_ylabel('y-axis')
# ax.set_zlabel('z-axis')

# displaying the plot
# plt.show()


# displaying the Plotly 3D plot
u = xs_fcc
v = ys_fcc
z = zs_fcc

points2D = np.vstack([u, v]).T
tri = Delaunay(points2D)
simplices = tri.simplices

fig1 = ff.create_trisurf(x=u, y=v, z=z,
                         colormap=['#FFBA33', '#FFBA33', '#FFBA33'],
                         show_colorbar=True,
                         simplices=simplices,
                         title="Boy's Surface")

u2 = xs_liquid
v2 = ys_liquid
z2 = zs_liquid

points2D_2 = np.vstack([u2, v2]).T
tri2 = Delaunay(points2D_2)
simplices2 = tri2.simplices

fig2 = ff.create_trisurf(x=u2, y=v2, z=z2,
                         simplices=simplices2,
                         colormap=['#33DEFF', '#33DEFF', '#33DEFF'],
                         title="Boy's Surface")

convex_hull_points = hull_points_projected
fig3 = go.Figure(
    data=[
        go.Scatter3d(
            x=convex_hull_points[:, 0],
            y=convex_hull_points[:, 1],
            z=convex_hull_points[:, 2],
            text=[f"idx: {idx}" for idx in hull_vertices_idx],
            mode="markers",
            marker=dict(size=1)
        )
    ]
)

u3 = convex_hull.points[:, 0]
v3 = convex_hull.points[:, 1]
z3 = convex_hull.points[:, 2]

hull_simplices = convex_hull.simplices


def is_on_lower_hull(input_arr, test_arr=hull_vertices_idx) -> bool:
    """See if all the elements in the input array exist in the test arr"""
    bool_arr = np.isin(input_arr, test_arr)
    return bool_arr.sum() == len(input_arr)


# create the index array for the ternary boundary points
res = 100
bound_2_idx = 0
bound_3_idx = bound_2_idx + res
bound_2_indices = [bound_2_idx]
bound_3_indices = [bound_3_idx]

for _ in range(1, res + 1):
    num_points = res - _
    bound_2_idx = bound_3_idx + 1
    bound_3_idx = bound_2_idx + num_points
    bound_2_indices.append(bound_2_idx)
    bound_3_indices.append(bound_3_idx)

bound_1_indices = np.array(
    list(range(res + 1)), dtype=np.int32
)
bound_2_indices = np.array(
    bound_2_indices, dtype=np.int32
)
bound_3_indices = np.array(
    bound_3_indices, dtype=np.int32
)

# add addition indices when there are more than 1 phase
num_phases = 2
num_points_per_phase = sum(range(1, res + 2))
ternary_boundary_indices = []
for indices in [bound_1_indices, bound_2_indices, bound_3_indices]:
    additional_indices = np.copy(indices)
    all_indices = [additional_indices]
    for _ in range(num_phases - 1):
        new_indices = additional_indices + num_points_per_phase
        all_indices.append(new_indices)
        additional_indices = new_indices
    new_boundary_indices = np.concatenate(all_indices)
    ternary_boundary_indices.append(new_boundary_indices)


def is_on_ternary_boundary(input_arr) -> bool:
    """See if all simplex vertex indices in the input array fall on the same ternary boundary
    i.e., excluding vertical simplices
    """
    for test_arr in ternary_boundary_indices:
        if is_on_lower_hull(input_arr, test_arr=test_arr):
            return True
    return False


hull_simplices_filtered = np.array(
    [
        indices for indices in hull_simplices if (is_on_lower_hull(indices) and not is_on_ternary_boundary(indices))
    ]
)

fig4 = ff.create_trisurf(x=u3, y=v3, z=z3,
                         simplices=hull_simplices_filtered,
                         colormap=['#33DEFF', '#33DEFF', '#33DEFF'])

fig4["data"][0].update(opacity=0.6)

# create the mapping from the convexhull.points to the selected vertices
keys = hull_vertices_idx
vals = np.arange(len(keys), dtype=np.int32)

mapping_arr = np.zeros(keys.max() + 1, dtype=vals.dtype)
mapping_arr[keys] = vals
hull_simplices_filtered_new_indexing = mapping_arr[hull_simplices_filtered]

fig5 = ff.create_trisurf(x=hull_points_projected[:, 0], y=hull_points_projected[:, 1], z=hull_points_projected[:, 2],
                         simplices=hull_simplices_filtered_new_indexing,
                         color_func=["#33DEFF"] * len(hull_simplices_filtered_new_indexing))

data = [fig1.data[0],  # fig1.data[1],
        fig2.data[0],  # fig2.data[1],
        fig3.data[0],
        fig4.data[0], fig4.data[1],
        fig5.data[1]]

iplot(
    dict(
        data=data,
        layout=dict(
            title=f"T={temp}K",
            scene_camera=dict(
                eye=dict(x=0, y=-1.25, z=1.25)
            )
        )
    )
)
