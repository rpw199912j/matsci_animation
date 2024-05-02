import math
import numpy as np
import pandas as pd

from manim import *
from itertools import combinations
from scipy.optimize import fsolve
from scipy.interpolate import CubicSpline
from scipy.spatial import ConvexHull

# config.frame_size = (1080, 1920)
config.frame_size = (2160, 3840)
config.frame_rate = 60
config.flush_cache = True
config.disable_caching = True


# cd .\thermodynamics\
# manim -p .\phase_diagram_cu_ag.py CuAgPhaseDiagram
class CuAgPhaseDiagram(Scene):
    def __init__(self):
        super().__init__()
        self.gibbs_df = pd.read_csv("../data/gibbs_energy/binary_Cu_Ag_fixed_ref.csv")
        self.g_x_axes = None
        self.paths = [
            [],  # alpha solidus
            [],  # alpha liquidus
            [],  # beta liquidus
            [],  # beta solidus
            [],  # alpha solvus
            []  # beta solvus
        ]

    @staticmethod
    def extract_gibbs_data(input_df, phase_name: str) -> tuple:
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

    @staticmethod
    def draw_phase_diagram_point(coord, ax, current_temp):
        if coord is None:
            return Dot().set_opacity(0)
        # AHA: be aware of object reference and modifying list in-place
        return Dot(
            point=ax.c2p(coord[0], current_temp)
        )

    @staticmethod
    def draw_phase_diagram_bound(path, ax):
        if not path:
            return Dot().set_opacity(0)
        # if not empty, draw the path
        path_arr = np.array(path)
        bound_path = ax.plot_line_graph(
            x_values=path_arr[:, 0],
            y_values=path_arr[:, 1],
            line_color=BLACK,
            add_vertex_dots=False,
            stroke_width=2
        )["line_graph"]

        return bound_path

    @staticmethod
    def draw_tangent_point(coord, ax):
        if coord is None:
            return Dot().set_opacity(0)

        return Dot(
            point=ax.c2p(*coord)
        )

    @staticmethod
    def get_line(point_left, point_right, ax):
        if point_left is None or point_right is None:
            return Dot().set_opacity(0)

        # unpack the coordinates
        x1, y1 = point_left
        x2, y2 = point_right

        # calculate the slope and the intercept
        m = (y2 - y1) / (x2 - x1)
        b = y1 - m * x1

        # get the common tangent line
        tangent_line = DashedVMobject(
            vmobject=ax.plot(
                lambda x: m * x + b,
                color=BLACK,
                stroke_width=2, stroke_opacity=0.5
            )
        )
        return tangent_line

    @staticmethod
    def make_phase_bound(arr_1, arr_2, extra_points=None, reverse_2nd=False):
        arr_3 = []
        if extra_points:
            arr_3 = extra_points
        if reverse_2nd:
            arr_2 = arr_2[::-1]

        phase_bound = arr_1 + arr_2 + arr_3
        return phase_bound

    @staticmethod
    def get_poly(phase_bound, ax, **kwargs):
        # convert the phase boundary array into the ax basis
        ax_coords = [
            ax.c2p(*coords) for coords in phase_bound
        ]
        return Polygon(
            *ax_coords,
            fill_opacity=0,
            **kwargs
        )

    def get_phase_regions(self, ax):
        alpha_fcc = self.make_phase_bound(
            self.paths[0], self.paths[4],
            extra_points=[[ax.x_range[0], ax.y_range[0]]]
        )
        alpha_fcc_liquid = self.make_phase_bound(
            self.paths[0], self.paths[1],
            reverse_2nd=True
        )
        liquid = self.make_phase_bound(
            self.paths[1], self.paths[2],
            extra_points=[[ax.x_range[1], ax.y_range[1]],
                          [ax.x_range[0], ax.y_range[1]]],
            reverse_2nd=True
        )
        beta_fcc_liquid = self.make_phase_bound(
            self.paths[2], self.paths[3],
            reverse_2nd=True
        )
        beta_fcc = self.make_phase_bound(
            self.paths[3], self.paths[5],
            extra_points=[[ax.x_range[1], ax.y_range[0]]]
        )
        alpha_beta_fcc = self.make_phase_bound(
            self.paths[4], self.paths[5],
            reverse_2nd=True
        )
        phase_bounds = [alpha_fcc, alpha_fcc_liquid, liquid,
                        beta_fcc_liquid, beta_fcc, alpha_beta_fcc]
        # phase_bound_colors = [ORANGE, GREEN, BLUE,
        #                       YELLOW, LIGHT_BROWN, PURPLE]
        phase_bound_colors = [BLACK] * len(phase_bounds)
        regions = [self.get_poly(phase_bound, ax, color=c)
                   for phase_bound, c in zip(phase_bounds, phase_bound_colors)]
        return regions

    @staticmethod
    def plot_phase_curve(spline_fit, ax, stroke_color):
        return ax.plot(
            lambda x: spline_fit(x),
            color=stroke_color, x_range=(ax.x_range[0], ax.x_range[1], 1e-3), use_smoothing=False
        )

    def construct(self):
        # Set the color theme to be the white background
        DecimalNumber.set_default(color=BLACK)
        Tex.set_default(color=BLACK)
        MathTex.set_default(color=BLACK)
        NumberLine.set_default(color=BLACK)
        Dot.set_default(color=BLACK)
        self.camera.background_color = WHITE

        # get all the temperature in descending order
        temps = self.gibbs_df["temp"].unique()[::-1]

        current_temp = temps[0]
        temp_tracker = ValueTracker(current_temp)

        # create two axes, one for G-X and one for T-X
        # def get_g_axes():
        #     _, fcc_y_vals = self.interpolate_gibbs_data(temp_tracker.get_value(), "FCC")
        #     _, liquid_y_vals = self.interpolate_gibbs_data(temp_tracker.get_value(), "LIQUID")
        #     all_y_vals = np.append(fcc_y_vals, liquid_y_vals)
        #     y_offset = 50
        #     y_min, y_max = np.min(all_y_vals) - y_offset, np.max(all_y_vals) + y_offset
        #
        #     # get the nearest 50 factors below the y_min
        #     y_10s_min = math.floor(y_min / 50) * 50
        #     # get the nearest 50 factors above the y_max
        #     y_10s_max = int(y_max / 50) * 50
        #
        #     # get the all the numbers that are multipliers of 50 in [y_min, y_max]
        #     y_tick_vals = np.arange(y_10s_min, y_10s_max + 50, 50, dtype=int)
        #     y_tick_vals = y_tick_vals[np.logical_and(y_tick_vals >= y_min, y_tick_vals <= y_max)]
        #
        #     ax = Axes(
        #         x_range=[0, 1, 0.1],
        #         y_range=[y_min, y_max, y_max - y_min],
        #         x_length=8,
        #         y_length=8,
        #         tips=False,
        #         y_axis_config={"include_ticks": False, "numbers_to_include": y_tick_vals,
        #                        "decimal_number_config": {"num_decimal_places": 0}},
        #         x_axis_config={"label_direction": UP, "include_numbers": True}
        #     ).shift(UP * 4.5)
        #
        #     # add custom ticks
        #     ax_y_axis = ax.y_axis
        #     ticks = VGroup()
        #     for _ in y_tick_vals:
        #         ticks.add(ax_y_axis.get_tick(_, ax_y_axis.tick_size))
        #     ax_y_axis.add(ticks)
        #     ax_y_axis.ticks = ticks
        #
        #     self.g_x_axes = ax
        #     return ax

        # use a fixed G-X axes
        y_min, y_max, y_step = -100, 50, 25
        y_tick_vals = np.arange(y_min, y_max + y_step, y_step, dtype=int)
        g_x_axes = Axes(
            x_range=[0, 1, 0.1],
            y_range=[y_min, y_max, y_max - y_min],
            x_length=8,
            y_length=8,
            tips=False,
            y_axis_config={"include_ticks": False, "numbers_to_include": y_tick_vals,
                           "decimal_number_config": {"num_decimal_places": 0},
                           "numbers_to_exclude": None},
            x_axis_config={"label_direction": UP, "include_numbers": True}
        ).shift(UP * 4.5)

        # shift the x-axis to the top of the y-axis
        x_axis = g_x_axes.get_x_axis()
        x_axis.shift(
            g_x_axes.c2p(0, y_max) - g_x_axes.c2p(0, 0)
        )

        # add custom ticks
        ax_y_axis = g_x_axes.y_axis
        ticks = VGroup()
        for _ in y_tick_vals:
            ticks.add(ax_y_axis.get_tick(_, ax_y_axis.tick_size))
        ax_y_axis.add(ticks)
        ax_y_axis.ticks = ticks

        self.g_x_axes = g_x_axes

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

        # add the thermo-calc logo
        thermocalc_logo = SVGMobject(
            r"../figure/ThermoCalc_logo.svg"
        ).set(height=config["frame_height"] * 0.10).set_color("#9b193b").shift(10 * UP + 5 * LEFT)
        self.play(
            DrawBorderThenFill(thermocalc_logo)
        )

        # add the chemistry label
        chemistry_label = Tex(
            "Cu-Ag", font_size=25
        ).next_to(thermocalc_logo, DOWN, buff=0.2)
        self.play(
            Write(chemistry_label), run_time=0.5
        )
        self.wait()

        # display the current temperature
        temp_label = Tex("T:", font_size=35).shift(10*UP + LEFT)
        temp_value = always_redraw(
            lambda:
            DecimalNumber(
                number=temp_tracker.get_value(),
                num_decimal_places=0,
                font_size=35,
                group_with_commas=False
            ).next_to(
                temp_label, RIGHT, buff=0.5*DEFAULT_MOBJECT_TO_MOBJECT_BUFFER
            ).align_to(
                temp_label, DOWN
            )
        )
        temp_unit = always_redraw(
            lambda:
            Tex(r"$^\circ$C", font_size=35).next_to(
                temp_value, RIGHT, buff=0.5 * DEFAULT_MOBJECT_TO_MOBJECT_BUFFER
            ).align_to(temp_label, DOWN)
        )
        temp_display = VGroup(temp_label, temp_value, temp_unit)
        self.play(
            Write(temp_display)
        )
        self.wait()
        
        # add the phase legend
        fcc_legend = Line(
            start=g_x_axes.c2p(0.05, 40), end=g_x_axes.c2p(0.15, 40), color=ORANGE
        )
        fcc_legend_label = Tex("FCC", font_size=35).next_to(fcc_legend, RIGHT)
        liquid_legend = Line(
            start=g_x_axes.c2p(0.05, 30), end=g_x_axes.c2p(0.15, 30), color=BLUE
        )
        liquid_legend_label = Tex("LIQUID", font_size=35).next_to(liquid_legend, RIGHT)

        self.play(
            GrowFromEdge(fcc_legend, LEFT),
            GrowFromEdge(liquid_legend, LEFT)
        )
        self.play(
            Write(fcc_legend_label),
            Write(liquid_legend_label)
        )
        self.wait()

        def get_all_mobjects():
            fcc_x_arr, fcc_y_arr = self.interpolate_gibbs_data(temp_tracker.get_value(), "FCC")
            liquid_x_arr, liquid_y_arr = self.interpolate_gibbs_data(temp_tracker.get_value(), "LIQUID")

            fcc_spline_fit = CubicSpline(fcc_x_arr, fcc_y_arr)
            liquid_spline_fit = CubicSpline(liquid_x_arr, liquid_y_arr)

            # fcc_phase_curve = g_x_axes.plot(
            #     lambda x: fcc_spline_fit(x),
            #     color=ORANGE, x_range=(g_x_axes.x_range[0], g_x_axes.x_range[1], 1e-3), use_smoothing=False
            # )
            fcc_phase_curve = self.plot_phase_curve(
                fcc_spline_fit, self.g_x_axes, ORANGE
            )

            # liquid_phase_curve = g_x_axes.plot(
            #     lambda x: liquid_spline_fit(x),
            #     color=BLUE, x_range=(g_x_axes.x_range[0], g_x_axes.x_range[1], 1e-3), use_smoothing=False
            # )
            liquid_phase_curve = self.plot_phase_curve(
                liquid_spline_fit, self.g_x_axes, BLUE
            )

            # determine the number of intersections between the fcc and liquid curves
            detect_intersect_range = np.arange(0, 1 + 1e-6, 1e-6)
            fcc_y_vals = fcc_spline_fit(detect_intersect_range)
            liquid_y_vals = liquid_spline_fit(detect_intersect_range)
            diff_y_vals = fcc_y_vals - liquid_y_vals

            # get the number of intersections
            num_intersections = (np.diff(diff_y_vals > 0) != 0).sum()

            # get the common tangent points via convex hull construction
            common_tangents = self.get_common_tangent_points(fcc_spline_fit, liquid_spline_fit)
            alpha_solidus, alpha_liquidus, beta_liquidus, beta_solidus, alpha_solvus, beta_solvus = (
                None, None, None, None, None, None
            )

            # if the fcc curve is above the liquid curve in the entire x range,
            # there is no intersection and only liquid phase
            if sum(diff_y_vals > 0) == len(diff_y_vals):
                print(f"\nOnly liquid @{temp_tracker.get_value()}")
            elif num_intersections == 1:
                print(f"\nOne intersection @{temp_tracker.get_value()}")
                alpha_solidus, alpha_liquidus = common_tangents
            elif num_intersections == 2:
                print(f"\nTwo intersections @{temp_tracker.get_value()}")
                # when there are 4 common tangent points
                if len(common_tangents) == 4:
                    alpha_solidus, alpha_liquidus, beta_liquidus, beta_solidus = common_tangents
                    # check if the alpha_liquidus and beta_liquidus is close enough
                    if np.isclose(alpha_liquidus[0], beta_liquidus[0], atol=0.001):
                        # add the eutectic line
                        left_end = [alpha_solidus[0], temp_tracker.get_value()]
                        right_end = [beta_solidus[0], temp_tracker.get_value()]
                        eutectic_line = Line(
                            start=t_x_axes.c2p(*left_end),
                            end=t_x_axes.c2p(*right_end),
                            color=BLACK,
                            stroke_width=2
                        )
                        self.add(eutectic_line)
                        self.paths[-2].append(left_end)
                        self.paths[-1].append(right_end)
                # when there are 2 common tangent points (only FCC phase)
                elif len(common_tangents) == 2:
                    alpha_solvus, beta_solvus = common_tangents

            elif sum(diff_y_vals < 0) == len(diff_y_vals):
                print(f"\nOnly fcc @{temp_tracker.get_value()}")
                alpha_solvus, beta_solvus = common_tangents

            updated_points = [alpha_solidus, alpha_liquidus, beta_liquidus, beta_solidus, alpha_solvus, beta_solvus]

            tangent_points = [
                self.draw_tangent_point(_, self.g_x_axes)
                for _ in updated_points
            ]

            phase_diagram_points = [
                self.draw_phase_diagram_point(_, t_x_axes, temp_tracker.get_value())
                for _ in updated_points
            ]

            connecting_vlines = [
                DashedLine(
                    start=g_x_point,
                    end=t_x_point,
                    stroke_width=2,
                    color=BLACK
                )
                for g_x_point, t_x_point in zip(tangent_points, phase_diagram_points)
            ]

            # construct the common tangent lines
            common_tangent_lines = [
                self.get_line(point_1, point_2, self.g_x_axes)
                for point_1, point_2 in [
                    (alpha_solidus, alpha_liquidus),
                    (beta_liquidus, beta_solidus),
                    (alpha_solvus, beta_solvus)
                ]
            ]

            for path, point in zip(self.paths, updated_points):
                # only add point coordinate when it exists (i.e., not None)
                if point is not None:
                    point[1] = temp_tracker.get_value()
                    path.append(point)

            phase_diagram_bounds = [
                self.draw_phase_diagram_bound(_, t_x_axes)
                for _ in self.paths
            ]

            return VGroup(
                fcc_phase_curve, liquid_phase_curve,
                *tangent_points, *phase_diagram_points,
                *common_tangent_lines,
                *connecting_vlines,
                *phase_diagram_bounds
            )

        all_mobs = always_redraw(get_all_mobjects)

        # TODO: investigating adding an always-updating VGroup and the creation with subgroups
        self.play(
            LaggedStart(
                Create(all_mobs[0]),
                Create(all_mobs[1]),
                lag_ratio=0.3
            )
        )
        self.add(all_mobs)
        self.wait()

        temp_hline = always_redraw(
            lambda:
            Line(
                start=t_x_axes.c2p(t_x_axes.x_range[0], temp_tracker.get_value()),
                end=t_x_axes.c2p(t_x_axes.x_range[1], temp_tracker.get_value()),
                color=BLACK, stroke_opacity=0.5
            )
        )
        self.play(
            Create(temp_hline)
        )
        self.wait()

        self.play(
            temp_tracker.animate(rate_func=linear, run_time=10).set_value(temps[-1])
        )
        self.wait()
        all_mobs.clear_updaters()
        # g_x_axes.clear_updaters()
        # remove visual clutter
        self.play(
            *[Uncreate(mob) for mob in all_mobs[2:-6]]
        )
        self.wait()

        # highlight different phase regions
        phase_regions = self.get_phase_regions(t_x_axes)
        phase_labels = [r"\alpha", r"\alpha+L", "L",
                        r"L+\beta", r"\beta", r"\alpha+\beta"]
        label_positions = [[0.04, 810], [0.23, 880], [0.5, 1030],
                           [0.86, 800], [0.96, 780], [0.5, 640]]
        font_sizes = [DEFAULT_FONT_SIZE, DEFAULT_FONT_SIZE, DEFAULT_FONT_SIZE,
                      0.7 * DEFAULT_FONT_SIZE, DEFAULT_FONT_SIZE, DEFAULT_FONT_SIZE]
        for label, label_pos, fontsize, phase_region in zip(phase_labels, label_positions, font_sizes, phase_regions):
            phase_region_label = MathTex(label, font_size=fontsize).move_to(
                t_x_axes.c2p(*label_pos)
            )
            phase_region = VGroup(phase_region, phase_region_label)
            target_phase_region = phase_region.copy().scale(1.2)
            target_phase_region[0].set_color(YELLOW).set_fill(opacity=0.6)
            final_phase_region = phase_region.copy()
            phase_region.set_opacity(0)
            self.play(
                Succession(
                    Transform(phase_region, target_phase_region, run_time=1),
                    Wait(run_time=0.2),
                    Transform(phase_region, final_phase_region, run_time=0.5)
                ),
                rate_func=smooth
            )
        self.wait()
