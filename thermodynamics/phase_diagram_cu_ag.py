import numpy as np
import pandas as pd

from manim import *
from itertools import combinations
from scipy.optimize import fsolve
from scipy.interpolate import CubicSpline
from scipy.spatial import ConvexHull

config.frame_size = (1080, 1920)
config.flush_cache = True
config.disable_caching = True


# cd .\thermodynamics\
# manim -p .\phase_diagram_cu_ag.py CuAgPhaseDiagram
class CuAgPhaseDiagram(Scene):
    def __init__(self):
        super().__init__()
        self.gibbs_df = pd.read_csv("../data/gibbs_energy/binary_Cu_Ag_mass_normalized.csv")

    def extract_gibbs_data(self, input_df, phase_name: str) -> tuple:
        """Extract the x, y data for a gibbs curve"""
        # get the x, y values for the given phase
        data_df = input_df[input_df["phase"] == phase_name]
        # sort and drop duplicates based on the x column
        data_df = data_df.sort_values(by="x", ignore_index=True).drop_duplicates(subset=["x"], ignore_index=True)
        return data_df["x"].to_numpy(), data_df["y"].to_numpy()

    def interpolate_gibbs_data(self, temp, phase_name: str) -> tuple:
        """Extract the x, y data for a gibbs curve"""
        # sorted unique temperature list
        temps_sorted = np.sort(self.gibbs_df["temp"].unique())
        # get the temperature ranges
        idx_upper_temp = np.searchsorted(temps_sorted, temp)
        idx_lower_temp = idx_upper_temp - 1
        # make sure the idx stays in the range of temps_sorted
        idx_lower_temp, idx_upper_temp = np.clip([idx_lower_temp, idx_upper_temp], 0, len(temps_sorted) - 1)
        lower_temp, upper_temp = temps_sorted[idx_lower_temp], temps_sorted[idx_upper_temp]

        # get the x, y data for the temperatures above and below
        lower_temp_df = self.gibbs_df[self.gibbs_df["temp"] == lower_temp]
        x_lower, y_lower = self.extract_gibbs_data(lower_temp_df, phase_name)
        upper_temp_df = self.gibbs_df[self.gibbs_df["temp"] == upper_temp]
        x_upper, y_upper = self.extract_gibbs_data(upper_temp_df, phase_name)

        # linearly interpolate a new line between the lower and upper temperature curves
        if upper_temp != lower_temp:
            y_new = y_lower + (temp - lower_temp) / (upper_temp - lower_temp) * (y_upper - y_lower)
        else:
            y_new = y_lower

        return x_lower, y_new

    def get_relevant_mobjects(self, current_temp, ax):
        """Helper function to get all relevant mobjects at a given temperature"""
        current_temp_df = self.gibbs_df[self.gibbs_df["temp"] == current_temp]
        fcc_x_arr, fcc_y_arr = self.extract_gibbs_data(current_temp_df, "FCC")
        liquid_x_arr, liquid_y_arr = self.extract_gibbs_data(current_temp_df, "LIQUID")
        # get the cubic spline fits
        fcc_spline_fit = CubicSpline(fcc_x_arr, fcc_y_arr)
        liquid_spline_fit = CubicSpline(liquid_x_arr, liquid_y_arr)

        fcc_curve = ax.plot(
            lambda x: fcc_spline_fit(x),
            color=ORANGE, x_range=(ax.x_range[0], ax.x_range[1], 1e-3), use_smoothing=False
        )
        liquid_curve = ax.plot(
            lambda x: liquid_spline_fit(x),
            color=BLUE, x_range=(ax.x_range[0], ax.x_range[1], 1e-3), use_smoothing=False
        )

        # determine the number of intersections between the fcc and liquid curves
        detect_intersect_range = np.arange(0, 1 + 1e-6, 1e-6)
        fcc_y_vals = fcc_spline_fit(detect_intersect_range)
        liquid_y_vals = liquid_spline_fit(detect_intersect_range)
        diff_y_vals = fcc_y_vals - liquid_y_vals

        # if the fcc curve is above the liquid curve in the entire x range,
        # there is no intersection and only liquid phase
        if sum(diff_y_vals > 0) == len(diff_y_vals):
            print(f"Only liquid @{current_temp}")

        return fcc_curve, liquid_curve

    @staticmethod
    def get_common_tangent_points(fcc_fit, liquid_fit):
        x_vals = np.arange(0, 1 + 1e-3, 1e-3)
        fcc_y_vals = fcc_fit(x_vals)
        liquid_y_vals = liquid_fit(x_vals)

        # make a convex hull with all the curves
        fcc_arr = np.array(list(zip(x_vals, fcc_y_vals)))
        liquid_arr = np.array(list(zip(x_vals, liquid_y_vals)))
        all_coords = np.append(fcc_arr, liquid_arr, axis=0)
        convex_hull = ConvexHull(all_coords)
        hull_vertex_indices = convex_hull.vertices

        # find the common tangent points
        # get all the possible combinations of the endpoint indices
        endpoint_combs = combinations([0, 1000, 1001, 2001], 2)
        endpoint_comb_sets = [set(comb) for comb in endpoint_combs]

        # take the difference between each pair of consecutive vertex indices
        consecutive_diff = np.diff(hull_vertex_indices)
        # check if the difference is not 1
        non_consecutive_indices = np.nonzero(consecutive_diff != 1)
        # initialize an empty list to store the indices for the common tangent points
        common_tangent_indices = []
        # iterate through all the non-zero differences
        for _ in np.nditer(non_consecutive_indices):
            # get the vertex indices corresponding to each gap
            facet_indices = hull_vertex_indices[_:_ + 2]
            # remove indices that are just endpoints
            # TODO: what if the common tangents are the endpoints?
            facet_indices_set = {*facet_indices}
            if facet_indices_set not in endpoint_comb_sets:
                common_tangent_indices.append(facet_indices)

        # flatten the common tangent indices
        common_tangent_indices_flatten = np.array(common_tangent_indices).flatten()
        # convert from indices to coordinates
        common_tangent_coords = [[*all_coords[_]] for _ in common_tangent_indices_flatten]
        return common_tangent_coords

    def construct(self):
        # create two axes, one for G-X and one for T-X
        g_x_axes = Axes(
            x_range=[0, 1, 0.1],
            y_range=[-1300, -300, 200],
            x_length=8,
            y_length=8,
            tips=False,
            axis_config={"include_numbers": True},
            x_axis_config={"label_direction": UP}
        ).shift(UP * 4.5)
        g_y_axis_label = g_x_axes.get_y_axis_label(
            Tex(r"Gibbs Energy (J/g)").rotate(PI / 2)
        ).next_to(g_x_axes.y_axis, LEFT)

        t_x_axes = Axes(
            x_range=[0, 1, 0.1],
            y_range=[500, 1200, 100],
            x_length=8,
            y_length=8,
            tips=False,
            axis_config={"include_numbers": True}
        ).shift(DOWN * 4.5)

        t_x_axis_label = t_x_axes.get_x_axis_label(
            Tex(r"$W_{\text{Ag}}$")
        ).next_to(t_x_axes.x_axis, DOWN)
        t_y_axis_label = t_x_axes.get_y_axis_label(
            Tex(r"Temperature ($^\circ$C)").rotate(PI / 2)
        ).next_to(t_x_axes.y_axis, LEFT)

        self.add(g_x_axes, t_x_axes, g_y_axis_label, t_x_axis_label, t_y_axis_label)
        self.wait()

        # get all the temperature in descending order
        temps = self.gibbs_df["temp"].unique()[::-1]

        current_temp = temps[0]
        # fcc_curve, liquid_curve = self.get_relevant_mobjects(current_temp, g_x_axes)
        # self.play(
        #     LaggedStart(
        #         Create(fcc_curve),
        #         Create(liquid_curve),
        #         lag_ratio=0.3
        #     )
        # )
        # self.wait()
        temp_tracker = ValueTracker(current_temp)

        def get_all_mobjects():
            fcc_x_arr, fcc_y_arr = self.interpolate_gibbs_data(temp_tracker.get_value(), "FCC")
            liquid_x_arr, liquid_y_arr = self.interpolate_gibbs_data(temp_tracker.get_value(), "LIQUID")

            fcc_spline_fit = CubicSpline(fcc_x_arr, fcc_y_arr)
            liquid_spline_fit = CubicSpline(liquid_x_arr, liquid_y_arr)
            fcc_slope_func = fcc_spline_fit.derivative()
            liquid_slope_func = liquid_spline_fit.derivative()

            fcc_phase_curve = g_x_axes.plot(
                lambda x: fcc_spline_fit(x),
                color=ORANGE, x_range=(g_x_axes.x_range[0], g_x_axes.x_range[1], 1e-3), use_smoothing=False
            )

            liquid_phase_curve = g_x_axes.plot(
                lambda x: liquid_spline_fit(x),
                color=BLUE, x_range=(g_x_axes.x_range[0], g_x_axes.x_range[1], 1e-3), use_smoothing=False
            )
            # determine the number of intersections between the fcc and liquid curves
            detect_intersect_range = np.arange(0, 1 + 1e-6, 1e-6)
            fcc_y_vals = fcc_spline_fit(detect_intersect_range)
            liquid_y_vals = liquid_spline_fit(detect_intersect_range)
            diff_y_vals = fcc_y_vals - liquid_y_vals

            # get the number of intersections
            num_intersections = (np.diff(diff_y_vals > 0) != 0).sum()

            # define a few systems of equations to solve the common tangent points
            # find the common tangent points between FCC and LIQUID
            def common_tangent_fcc_liquid(x):
                """Systems of non-linear functions to solve for the common tangent points"""
                return [fcc_slope_func(x[0]) - liquid_slope_func(x[1]),
                        (liquid_spline_fit(x[1]) - fcc_spline_fit(x[0])) / (x[1] - x[0]) - fcc_slope_func(x[0])]

            # find the common tangent points between FCC and FCC itself
            def common_tangent_fcc_only(x):
                """Systems of non-linear functions to solve for the common tangent points"""
                return [fcc_slope_func(x[0]) - fcc_slope_func(x[1]),
                        (fcc_spline_fit(x[1]) - fcc_spline_fit(x[0])) / (x[1] - x[0]) - fcc_slope_func(x[0])]

            # if the fcc curve is above the liquid curve in the entire x range,
            # there is no intersection and only liquid phase
            if sum(diff_y_vals > 0) == len(diff_y_vals):
                print(f"\nOnly liquid @{temp_tracker.get_value()}")
            elif num_intersections == 1:
                print(f"\nOne intersection @{temp_tracker.get_value()}")
                # get the x, y coordinates of the tangent points
                fcc_tangent_x, liquid_tangent_x = fsolve(common_tangent_fcc_liquid, [0.05, 0.1])
                fcc_tangent_y, liquid_tangent_y = fcc_spline_fit(fcc_tangent_x), liquid_spline_fit(liquid_tangent_x)
                print(f"\ntangent 1: {fcc_tangent_x, fcc_tangent_y}; tangent 2: {liquid_tangent_x, liquid_tangent_y}")
            elif num_intersections == 2:
                print(f"\nTwo intersections @{temp_tracker.get_value()}")
                # get the x, y coordinates of the tangent points
                fcc_tangent_x_1, liquid_tangent_x_1 = fsolve(common_tangent_fcc_liquid, [0.05, 0.1])
                fcc_tangent_y_1, liquid_tangent_y_1 = fcc_spline_fit(fcc_tangent_x_1), liquid_spline_fit(liquid_tangent_x_1)
                print(f"\ntangent 1: {fcc_tangent_x_1, fcc_tangent_y_1}; tangent 2: {liquid_tangent_x_1, liquid_tangent_y_1}")

                # get the x, y coordinates of the tangent points
                fcc_tangent_x_2, liquid_tangent_x_2 = fsolve(common_tangent_fcc_liquid, [0.9, 0.95])
                fcc_tangent_y_2, liquid_tangent_y_2 = fcc_spline_fit(fcc_tangent_x_2), liquid_spline_fit(
                    liquid_tangent_x_2)
                print(
                    f"\ntangent 3: {fcc_tangent_x_2, fcc_tangent_y_2}; tangent 4: {liquid_tangent_x_2, liquid_tangent_y_2}")
            elif sum(diff_y_vals < 0) == len(diff_y_vals):
                print(f"\nOnly fcc @{temp_tracker.get_value()}")

                # get the x, y coordinates of the tangent points
                fcc_tangent_x_1, fcc_tangent_x_2 = fsolve(common_tangent_fcc_only, [0.05, 0.95])
                fcc_tangent_y_1, fcc_tangent_y_2 = fcc_spline_fit(fcc_tangent_x_1), fcc_spline_fit(fcc_tangent_x_2)
                print(f"\ntangent 1: {fcc_tangent_x_1, fcc_tangent_y_1}; tangent 2: {fcc_tangent_x_2, fcc_tangent_y_2}")

            return VGroup(
                fcc_phase_curve, liquid_phase_curve
            )

        all_mobs = always_redraw(get_all_mobjects)

        self.play(
            LaggedStart(
                Create(all_mobs),
                lag_ratio=0.3
            )
        )
        self.wait()

        temp_hline = always_redraw(
            lambda:
            Line(
                start=t_x_axes.c2p(t_x_axes.x_range[0], temp_tracker.get_value()),
                end=t_x_axes.c2p(t_x_axes.x_range[1], temp_tracker.get_value()),
                color=WHITE, stroke_opacity=0.5
            )
        )
        self.play(
            Create(temp_hline)
        )
        self.wait()

        self.play(
            temp_tracker.animate(rate_func=linear, run_time=5).set_value(temps[-1])
        )
        self.wait()

        # for current_temp in temps[1:]:
        #     current_fcc_curve, current_liquid_curve = self.get_relevant_mobjects(current_temp, g_x_axes)
        #     self.play(
        #         Transform(fcc_curve, current_fcc_curve),
        #         Transform(liquid_curve, current_liquid_curve),
        #         temp_tracker.animate.set_value(current_temp),
        #         rate_func=linear
        #     )
