import numpy as np
import pandas as pd

from manim import *
from scipy.optimize import fsolve
from scipy.interpolate import CubicSpline


# noinspection DuplicatedCode
class BinaryDrivingForce(ZoomedScene):
    def __init__(self):
        ZoomedScene.__init__(
            self,
            zoomed_display_width=8,
            zoomed_display_height=4,
            zoom_factor=0.15
        )

    @staticmethod
    def extract_gibbs_data(input_df, phase_name: str) -> tuple:
        """Extract the x, y data for a gibbs curve"""
        # get the x, y values for the given phase
        data_df = input_df[input_df["phase"] == phase_name]
        # sort and drop duplicates based on the x column
        data_df = data_df.sort_values(by="x", ignore_index=True).drop_duplicates(subset=["x"], ignore_index=True)
        return data_df["x"], data_df["y"]

    def get_relevant_mobjects(self,
                              fcc_x_arr, bcc_x_arr, metastable_x_arr,
                              fcc_y_arr, bcc_y_arr, metastable_y_arr,
                              ax, ax_x_min, ax_x_max, x_mid=0.25):
        """Helper function to get all relevant mobjects at a given temperature"""
        # get the cubic spline fits
        fcc_spline_fit = CubicSpline(fcc_x_arr, fcc_y_arr)
        bcc_spline_fit = CubicSpline(bcc_x_arr, bcc_y_arr)
        metastable_spline_fit = CubicSpline(metastable_x_arr, metastable_y_arr)

        fcc_curve = ax.plot(
            lambda x: fcc_spline_fit(x),
            color=BLUE, x_range=(ax_x_min, ax_x_max, 1e-4), use_smoothing=False
        )
        bcc_curve = ax.plot(
            lambda x: bcc_spline_fit(x),
            color=GREEN, x_range=(ax_x_min, ax_x_max, 1e-4), use_smoothing=False
        )

        metastable_curve = ax.plot(
            lambda x: metastable_spline_fit(x),
            color=RED, x_range=(ax_x_min, ax_x_max, 1e-4), use_smoothing=False
        )

        # get the corresponding derivative functions
        fcc_slope_func = fcc_spline_fit.derivative()
        bcc_slope_func = bcc_spline_fit.derivative()
        metastable_slope_func = metastable_spline_fit.derivative()

        # define an array of mole fractions to get the derivatives of the fitted gibbs spline curves
        x_spline_arr = np.arange(0, 1, 1e-6)
        # define the middle point on the fcc phase to calculate driving force to the other two phases
        fcc_slope_x_mid = fcc_slope_func(x_mid)

        # find the parallel tangent at BCC
        bcc_parallel_index = np.argmin(np.abs(bcc_slope_func(x_spline_arr) - fcc_slope_x_mid))
        bcc_parallel_x = x_spline_arr[bcc_parallel_index]

        # find the parallel tangent at the metastable phase
        metastable_x_spline_arr = np.arange(0.5, 1, 1e-6)
        metastable_parallel_index = np.argmin(np.abs(metastable_slope_func(metastable_x_spline_arr) - fcc_slope_x_mid))
        metastable_parallel_x = metastable_x_spline_arr[metastable_parallel_index]

        # find the common tangent points between FCC and BCC
        def common_tangent(x):
            """Systems of non-linear functions to solve for the common tangent points"""
            return [fcc_slope_func(x[0]) - bcc_slope_func(x[1]),
                    (bcc_spline_fit(x[1]) - fcc_spline_fit(x[0])) / (x[1] - x[0]) - fcc_slope_func(x[0])]

        # get the x, y coordinates of the tangent points
        fcc_tangent_x, bcc_tangent_x = fsolve(common_tangent, [0.1, 0.3])
        fcc_tangent_y, bcc_tangent_y = fcc_spline_fit(fcc_tangent_x), bcc_spline_fit(bcc_tangent_x)

        # get the common tangent line
        common_tangent_slope, common_tangent_intercept = self.get_line_from_points(fcc_tangent_x, fcc_tangent_y,
                                                                                   bcc_tangent_x, bcc_tangent_y)
        common_tangent_line = DashedVMobject(
            ax.plot(
                lambda x: common_tangent_slope * x + common_tangent_intercept,
                x_range=(ax_x_min, ax_x_max), stroke_opacity=0.4
            ),
            num_dashes=80, stroke_width=2
        )

        dot_radius = 0.05
        # get the mid_point
        y_mid = fcc_spline_fit(x_mid)
        mid_point = Dot(
            point=ax.c2p(x_mid, y_mid), radius=dot_radius, color=BLUE
        )
        # get the tangent line at the mid_point
        mid_tangent_slope, mid_tangent_intercept = self.get_tangent_line(x_mid, y_mid, fcc_slope_x_mid)
        mid_tangent_line = DashedVMobject(
            ax.plot(
                lambda x: mid_tangent_slope * x + mid_tangent_intercept,
                x_range=(ax_x_min, ax_x_max), color=BLUE
            ),
            num_dashes=80, stroke_width=2
        )

        # get the parallel tangent at BCC
        bcc_parallel_y = bcc_spline_fit(bcc_parallel_x)
        bcc_parallel_slope, bcc_parallel_intercept = self.get_tangent_line(bcc_parallel_x,
                                                                           bcc_parallel_y,
                                                                           bcc_slope_func(bcc_parallel_x))
        bcc_parallel_line = DashedVMobject(
            ax.plot(
                lambda x: bcc_parallel_slope * x + bcc_parallel_intercept,
                x_range=(ax_x_min, ax_x_max), color=GREEN
            ),
            num_dashes=80, stroke_width=2
        )
        bcc_to_fcc_y = mid_tangent_slope * bcc_parallel_x + mid_tangent_intercept
        bcc_connect_line = Line(
            start=ax.c2p(bcc_parallel_x, bcc_to_fcc_y),
            end=ax.c2p(bcc_parallel_x, bcc_parallel_y),
            color=GREEN
        )

        # get the parallel tangent at the metastable phase
        metastable_parallel_y = metastable_spline_fit(metastable_parallel_x)
        metastable_parallel_slope, metastable_parallel_intercept = self.get_tangent_line(metastable_parallel_x,
                                                                                         metastable_parallel_y,
                                                                                         metastable_slope_func(
                                                                                             metastable_parallel_x))
        metastable_parallel_line = DashedVMobject(
            ax.plot(
                lambda x: metastable_parallel_slope * x + metastable_parallel_intercept,
                x_range=(ax_x_min, ax_x_max), color=RED
            ),
            num_dashes=80, stroke_width=2
        )
        metastable_to_fcc_y = mid_tangent_slope * metastable_parallel_x + mid_tangent_intercept
        metastable_connect_line = Line(
            start=ax.c2p(metastable_parallel_x, metastable_to_fcc_y),
            end=ax.c2p(metastable_parallel_x, metastable_parallel_y),
            color=RED
        )

        # get the driving forces
        bcc_driving_force = bcc_parallel_y - bcc_to_fcc_y
        metastable_driving_force = metastable_parallel_y - metastable_to_fcc_y

        return (bcc_curve, fcc_curve, metastable_curve, common_tangent_line,
                mid_point, mid_tangent_line,
                bcc_parallel_line, bcc_connect_line,
                metastable_parallel_line, metastable_connect_line,
                bcc_driving_force, metastable_driving_force)

    @staticmethod
    def get_line_from_points(x1, y1, x2, y2):
        slope = (y2 - y1) / (x2 - x1)
        b = y1 - slope * x1
        return slope, b

    @staticmethod
    def get_tangent_line(x, y, slope) -> tuple:
        b = y - slope * x
        return slope, b

    def construct(self):
        # define the axis attributes
        x_min, x_max, x_step = 0, 1, 0.1
        y_min, y_max, y_step = -70000, -5000, 13000
        axes = Axes(
            x_range=[x_min, x_max, x_step],
            y_range=[y_min, y_max, y_step],
            x_length=8,
            y_length=6,
            tips=False
        ).shift(LEFT)

        self.add(axes)
        self.wait()

        # get the axis
        x_axis = axes.get_x_axis()
        y_axis = axes.get_y_axis()

        # add the axes labels
        x_label = axes.get_x_axis_label(Tex(r"$X_{\text{Be}}$")).shift(0.2 * RIGHT + 0.5 * DOWN)
        y_label = axes.get_y_axis_label(
            Tex("$G$ (kJ/mol)").rotate(PI / 2), edge=LEFT, direction=LEFT, buff=1
        )

        # add custom axes tick labels
        values_x = [
            (tick_val, str(tick_val))
            for tick_val in np.arange(x_min * 10, (x_max + x_step) * 10, x_step * 10) / 10
        ]
        values_y = [
            (tick_val, str(int(tick_val / 1e3)))
            for tick_val in np.arange(y_min, y_max + y_step, y_step)
        ]

        x_axis_labels = VGroup()
        y_axis_labels = VGroup()

        for x_val, x_tex in values_x:
            tex = Tex(x_tex, font_size=35)  # Convert string to tex
            tex.next_to(x_axis.n2p(x_val), UP)  # Put tex on the position
            x_axis_labels.add(tex)

        for y_val, y_tex in values_y:
            tex = Tex(y_tex, font_size=35)  # Convert string to tex
            tex.next_to(y_axis.n2p(y_val), LEFT)  # Put tex on the position
            y_axis_labels.add(tex)

        self.play(
            Write(x_label),
            Write(y_label)
        )
        self.wait()
        self.play(
            Write(x_axis_labels),
            Write(y_axis_labels),
            run_time=0.8
        )
        self.wait()

        # get the gibbs data
        gibbs_df = pd.read_csv(r"../data/gibbs_energy/binary_Cu_Be.csv")
        # get the list of temperatures
        temp_lst = gibbs_df["temp"].unique()
        # get the data at the first temperature step
        prev_temp_data = gibbs_df[gibbs_df["temp"] == temp_lst[0]]
        # get the data corresponding to each phase
        fcc_x, fcc_y = self.extract_gibbs_data(prev_temp_data, "FCC")
        bcc_x, bcc_y = self.extract_gibbs_data(prev_temp_data, "BCC")
        metastable_x, metastable_y = self.extract_gibbs_data(prev_temp_data, "BE2CU")

        (bcc_curve, fcc_curve, metastable_curve, common_tangent_line,
         mid_point, mid_tangent_line,
         bcc_parallel_line, bcc_connect_line,
         metastable_parallel_line, metastable_connect_line,
         bcc_driving_force, metastable_driving_force) = self.get_relevant_mobjects(
            fcc_x, bcc_x, metastable_x,
            fcc_y, bcc_y, metastable_y,
            axes, x_min, x_max
        )

        # show the gibbs curves
        self.play(
            LaggedStart(
                Create(bcc_curve),
                Create(fcc_curve),
                Create(metastable_curve),
                lag_ratio=0.8
            )
        )
        self.wait()

        # display the color legend
        bcc_legend = Line(
            start=axes.c2p(0.05, -10000), end=axes.c2p(0.15, -10000), color=BLUE
        )
        bcc_legend_label = Tex("BCC", font_size=35).next_to(bcc_legend, RIGHT)
        fcc_legend = Line(
            start=axes.c2p(0.05, -15000), end=axes.c2p(0.15, -15000), color=GREEN
        )
        fcc_legend_label = Tex("FCC", font_size=35).next_to(fcc_legend, RIGHT)
        metastable_legend = Line(
            start=axes.c2p(0.05, -20000), end=axes.c2p(0.15, -20000), color=RED
        )
        metastable_legend_label = Tex("Be$_2$Cu", font_size=35).next_to(metastable_legend, RIGHT)

        self.play(
            GrowFromEdge(bcc_legend, LEFT),
            GrowFromEdge(fcc_legend, LEFT),
            GrowFromEdge(metastable_legend, LEFT)
        )
        self.play(
            Write(bcc_legend_label),
            Write(fcc_legend_label),
            Write(metastable_legend_label)
        )
        self.wait()

        # show the common tangent line
        self.play(
            FadeIn(common_tangent_line, shift=UP * 0.2)
        )
        self.wait()

        # show the mid_point and the tangent line
        self.play(
            DrawBorderThenFill(mid_point)
        )
        self.play(
            Create(mid_tangent_line)
        )
        self.wait()

        # show the parallel tangent for BCC
        self.play(
            TransformFromCopy(mid_tangent_line, bcc_parallel_line),
            GrowFromEdge(bcc_connect_line, UP)
        )
        self.wait()

        # show the parallel tangent for the metastable phase
        self.play(
            TransformFromCopy(mid_tangent_line, metastable_parallel_line),
            GrowFromEdge(metastable_connect_line, UP)
        )
        self.wait()

        # define a second axis to show the driving force
        t_step = temp_lst[1] - temp_lst[0]
        t_min = temp_lst[0] - t_step
        t_max = temp_lst[-1] + t_step

        g_min, g_max = -5000, 3000
        driving_force_axis = Axes(
            x_range=[t_min, t_max, t_max - t_min],
            y_range=[g_min, g_max, g_max - g_min],
            x_length=2.5,
            y_length=2,
            tips=False
        ).to_edge(RIGHT)

        driving_force_y_axis = driving_force_axis.get_y_axis()

        # add y_ticks inplace
        y_axis = driving_force_axis.y_axis
        ticks = VGroup()
        for _ in [g_min, 0, g_max]:
            ticks.add(y_axis.get_tick(_, y_axis.tick_size))
        y_axis.add(ticks)
        y_axis.ticks = ticks

        # get the axis labels
        driving_force_x_label = driving_force_axis.get_x_axis_label(
            MathTex(r"T (^\circ C)", font_size=25)
        ).move_to(driving_force_axis.c2p(t_max, -1000))
        driving_force_y_label = driving_force_axis.get_y_axis_label(
            Tex(r"$\Delta G$ (kJ/mol)", font_size=25).rotate(PI / 2), edge=LEFT, direction=LEFT, buff=0.4
        )

        # get custom tick labels
        values_y = [
            (tick_val, str(int(tick_val / 1e3)))
            for tick_val in [g_min, 0, g_max]
        ]

        y_axis_labels = VGroup()

        for y_val, y_tex in values_y:
            tex = Tex(y_tex, font_size=25)  # Convert string to tex
            tex.next_to(driving_force_y_axis.n2p(y_val), LEFT)  # Put tex on the position
            y_axis_labels.add(tex)

        self.play(
            FadeIn(driving_force_axis),
            Write(driving_force_x_label),
            Write(driving_force_y_label),
            Write(y_axis_labels)
        )
        self.wait()

        # set a temperature tracker
        current_temp = ValueTracker(temp_lst[0])
        bcc_driving_force_tracker = ValueTracker(bcc_driving_force)
        metastable_driving_force_tracker = ValueTracker(metastable_driving_force)

        # display the temperature
        tick_vert_offset = 200
        temp_tick = always_redraw(
            lambda:
            Line(
                start=driving_force_axis.c2p(current_temp.get_value(), tick_vert_offset),
                end=driving_force_axis.c2p(current_temp.get_value(), -tick_vert_offset)
            )
        )
        temp_tick_label = always_redraw(
            lambda:
            DecimalNumber(
                number=current_temp.get_value(),
                num_decimal_places=0,
                font_size=DEFAULT_FONT_SIZE * 0.5
            ).next_to(temp_tick, UP)
        )

        temp_tick_x_min = Line(
            start=driving_force_axis.c2p(current_temp.get_value(), tick_vert_offset),
            end=driving_force_axis.c2p(current_temp.get_value(), -tick_vert_offset)
        )
        temp_tick_label_x_min = DecimalNumber(
            number=current_temp.get_value(),
            num_decimal_places=0,
            font_size=DEFAULT_FONT_SIZE * 0.5
        ).next_to(temp_tick, UP)
        self.play(
            LaggedStart(
                Create(temp_tick),
                Write(temp_tick_label),
                lag_ratio=0.8
            )
        )
        self.add(temp_tick_x_min, temp_tick_label_x_min)

        # display the driving force
        bcc_driving_force_dot = always_redraw(
            lambda:
            Dot(
                point=driving_force_axis.c2p(current_temp.get_value(), bcc_driving_force_tracker.get_value()),
                color=GREEN, radius=0.06
            )
        )
        metastable_driving_force_dot = always_redraw(
            lambda:
            Dot(
                point=driving_force_axis.c2p(current_temp.get_value(), metastable_driving_force_tracker.get_value()),
                color=RED, radius=0.06
            )
        )
        self.play(
            DrawBorderThenFill(bcc_driving_force_dot),
            DrawBorderThenFill(metastable_driving_force_dot),
        )

        self.wait()

        # add the traced paths
        bcc_driving_force_path = TracedPath(
            bcc_driving_force_dot.get_center,
            stroke_color=GREEN, stroke_width=4
        )
        metastable_driving_force_path = TracedPath(
            metastable_driving_force_dot.get_center,
            stroke_color=RED, stroke_width=4
        )
        self.add(bcc_driving_force_path, metastable_driving_force_path)

        for temp in temp_lst[1:]:
            current_temp_data = gibbs_df[gibbs_df["temp"] == temp]
            # get the data corresponding to each phase
            fcc_x, fcc_y = self.extract_gibbs_data(current_temp_data, "FCC")
            bcc_x, bcc_y = self.extract_gibbs_data(current_temp_data, "BCC")
            metastable_x, metastable_y = self.extract_gibbs_data(current_temp_data, "BE2CU")

            (current_bcc_curve, current_fcc_curve, current_metastable_curve, current_common_tangent_line,
             current_mid_point, current_mid_tangent_line,
             current_bcc_parallel_line, current_bcc_connect_line,
             current_metastable_parallel_line, current_metastable_connect_line,
             current_bcc_driving_force, current_metastable_driving_force) = self.get_relevant_mobjects(
                fcc_x, bcc_x, metastable_x,
                fcc_y, bcc_y, metastable_y,
                axes, x_min, x_max
            )

            self.play(
                Transform(bcc_curve, current_bcc_curve),
                Transform(fcc_curve, current_fcc_curve),
                Transform(metastable_curve, current_metastable_curve),
                Transform(common_tangent_line, current_common_tangent_line),
                Transform(mid_point, current_mid_point),
                Transform(mid_tangent_line, current_mid_tangent_line),
                Transform(bcc_parallel_line, current_bcc_parallel_line),
                Transform(bcc_connect_line, current_bcc_connect_line),
                Transform(metastable_parallel_line, current_metastable_parallel_line),
                Transform(metastable_connect_line, current_metastable_connect_line),
                current_temp.animate.set_value(temp),
                bcc_driving_force_tracker.animate.set_value(current_bcc_driving_force),
                metastable_driving_force_tracker.animate.set_value(current_metastable_driving_force),
                rate_func=linear
            )

        self.wait()
