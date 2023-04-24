import numpy as np
from itertools import combinations
from manim import *
from typing import List
from scipy.spatial import ConvexHull

config.disable_caching = True
config.flush_cache = True


class CommonTangent(Scene):
    @staticmethod
    def get_axis_ticks_labels(values: List[tuple], axis: Axes, fontsize=20, label_dir=UP) -> VGroup:
        axis_labels = VGroup()
        for tick_val, tick_str in values:
            tex = Tex(tick_str, font_size=fontsize)
            tex.next_to(axis.n2p(tick_val), label_dir)
            axis_labels.add(tex)
        return axis_labels

    @staticmethod
    def add_axis_ticks(values, axis, inplace=True):
        """Add axis ticks inplace"""
        ticks = VGroup()
        for _ in values:
            ticks.add(axis.get_tick(_, axis.tick_size))
        if inplace:
            axis.add(ticks)
            axis.ticks = ticks
        return ticks

    @staticmethod
    def add_nl_ticks(values, nl, inplace=True):
        """Add NumberLine ticks inplace"""
        ticks = VGroup()
        for _ in values:
            ticks.add(nl.get_tick(_, nl.tick_size))
        if inplace:
            nl.add(ticks)
            nl.ticks = ticks
        return ticks

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

    def construct(self):
        # Set the color theme to be the white background
        DecimalNumber.set_default(color=BLACK)
        Tex.set_default(color=BLACK)
        MathTex.set_default(color=BLACK)
        NumberLine.set_default(color=BLACK)
        Line.set_default(color=BLACK)
        Dot.set_default(color=BLACK)
        self.camera.background_color = WHITE

        x_min, x_max, x_step = 0, 1, 0.1
        y_min, y_max, y_step = -3000, 2000, 1000
        axes = Axes(
            x_range=[x_min, x_max, x_step],
            y_range=[y_min, y_max, y_step],
            x_length=6,
            y_length=6.5,
            axis_config={"include_ticks": False},
            tips=False
        ).shift(UP * 0.3)

        # get the x, y axes separately
        x_axis = axes.get_x_axis().copy().shift(
            axes.c2p(axes.x_range[0], axes.y_range[0]) - axes.c2p(axes.x_range[0], 0)
        )
        y_axis = axes.get_y_axis()

        # add the ticks and tick labels
        values_x = [
            (tick_val, str(tick_val))
            for tick_val in np.arange(x_min * 10, (x_max + x_step) * 10, x_step * 10) / 10
        ]
        values_y = [
            (tick_val, str(int(tick_val / 1e3)))
            for tick_val in np.arange(y_min, y_max + y_step, y_step)
        ]

        x_tick_labels = self.get_axis_ticks_labels(values_x, x_axis, fontsize=25, label_dir=DOWN)
        y_tick_labels = self.get_axis_ticks_labels(values_y, y_axis, fontsize=25, label_dir=LEFT)
        x_ticks = self.add_axis_ticks([val[0] for val in values_x], x_axis, inplace=False)
        self.add_axis_ticks([val[0] for val in values_y], axes.y_axis)

        # get the up and right bounding line
        up_bound_line = Line(
            start=axes.c2p(x_min, y_max),
            end=axes.c2p(x_max, y_max),
            stroke_width=2
        )
        right_bound_line = Line(
            start=axes.c2p(x_max, y_min),
            end=axes.c2p(x_max, y_max),
            stroke_width=2
        )

        # add the axis titles
        x_label = axes.get_x_axis_label(
            MathTex(r"X_\text{B}", font_size=35)
        ).next_to(x_tick_labels, DOWN)
        y_label = axes.get_y_axis_label(
            Tex("$g$ (kJ/mol)", font_size=35).rotate(PI / 2)
        ).next_to(y_tick_labels, LEFT)

        self.add(
            x_axis, y_axis, up_bound_line, right_bound_line, x_ticks,
            x_label,
            y_label,
            x_tick_labels,
            y_tick_labels
        )
        self.wait()

        # add a dashed hline for the G=0
        zero_hline = DashedLine(
            start=axes.c2p(x_min, 0),
            end=axes.c2p(x_max, 0),
            stroke_width=2
        )

        # add the tick labels
        self.play(
            Create(zero_hline)
        )
        self.wait()

        solid_legend = Line(
            start=axes.c2p(1.15, -500), end=axes.c2p(1.25, -500), color=ORANGE
        )
        solid_legend_label = Tex("Solid", font_size=30).next_to(solid_legend, RIGHT)
        temp_label = Tex("At fixed T", font_size=30).next_to(
            solid_legend, UP, buff=1.2*DEFAULT_MOBJECT_TO_MOBJECT_BUFFER).align_to(solid_legend, LEFT)
        liquid_legend = Line(
            start=axes.c2p(1.15, -800), end=axes.c2p(1.25, -800), color=BLUE
        )
        liquid_legend_label = Tex("Liquid", font_size=30).next_to(liquid_legend, RIGHT)

        self.play(
            Write(temp_label)
        )
        self.wait()
        self.play(
            GrowFromEdge(solid_legend, LEFT),
            GrowFromEdge(liquid_legend, LEFT)
        )
        self.play(
            Write(solid_legend_label),
            Write(liquid_legend_label)
        )
        self.wait()

        # define relevant variables
        liquid_a_ref_tracker = ValueTracker(0)
        solid_b_ref_tracker = ValueTracker(0)
        delta_mu_a = 1000
        delta_mu_b = 1500
        temp_tracker = ValueTracker(436)
        omega_solid_tracker = ValueTracker(4400)
        omega_liquid_tracker = ValueTracker(2370)

        # define a series of gibbs functions
        def get_ref_gibbs(x, a_ref, b_ref):
            gibbs = a_ref * (1 - x) + b_ref * x
            return gibbs

        def get_ideal_gibbs(x):
            if x in {0, 1}:
                return 0
            gibbs = 8.31 * temp_tracker.get_value() * (x * np.log(x) + (1 - x) * np.log(1 - x))
            return gibbs

        def get_excess_gibbs(x, interaction):
            gibbs = (1 - x) * x * interaction
            return gibbs

        @np.vectorize
        def get_solid_gibbs(x):
            ref_gibbs = get_ref_gibbs(x,
                                      liquid_a_ref_tracker.get_value() + delta_mu_a,
                                      solid_b_ref_tracker.get_value())
            ideal_gibbs = get_ideal_gibbs(x)
            excess_gibbs = get_excess_gibbs(x, omega_solid_tracker.get_value())
            return ref_gibbs + ideal_gibbs + excess_gibbs

        @np.vectorize
        def get_liquid_gibbs(x):
            ref_gibbs = get_ref_gibbs(x,
                                      liquid_a_ref_tracker.get_value(),
                                      solid_b_ref_tracker.get_value() + delta_mu_b)
            ideal_gibbs = get_ideal_gibbs(x)
            excess_gibbs = get_excess_gibbs(x, omega_liquid_tracker.get_value())
            return ref_gibbs + ideal_gibbs + excess_gibbs

        solid_gibbs = always_redraw(
            lambda:
            axes.plot(
                lambda x: get_solid_gibbs(x),
                color=ORANGE
            )
        )

        liquid_gibbs = always_redraw(
            lambda:
            axes.plot(
                lambda x: get_liquid_gibbs(x),
                color=BLUE
            )
        )

        def get_mu_label(label, x, y, direct):
            return MathTex(label, font_size=25).next_to(
                axes.c2p(x, y), direction=direct, buff=0.8*DEFAULT_MOBJECT_TO_MOBJECT_BUFFER
            )

        mu_a_solid = always_redraw(
            lambda:
            get_mu_label(
                label=r"\mu_{\text{A}}^{o\text{S}}",
                x=0,
                y=get_solid_gibbs(0),
                direct=RIGHT
            )
        )

        mu_a_liquid = always_redraw(
            lambda:
            get_mu_label(
                label=r"\mu_{\text{A}}^{o\text{L}}",
                x=0,
                y=get_liquid_gibbs(0),
                direct=RIGHT
            )
        )

        mu_b_solid = always_redraw(
            lambda:
            get_mu_label(
                label=r"\mu_{\text{B}}^{o\text{S}}",
                x=1,
                y=get_solid_gibbs(1),
                direct=LEFT
            )
        )

        mu_b_liquid = always_redraw(
            lambda:
            get_mu_label(
                label=r"\mu_{\text{B}}^{o\text{L}}",
                x=1,
                y=get_liquid_gibbs(1),
                direct=LEFT
            )
        )

        # Create the common tangent points and connecting line
        def get_tangent():
            common_tangents = self.get_common_tangent_points(get_solid_gibbs, get_liquid_gibbs)

            tangent_points = [
                self.draw_tangent_point(_, axes)
                for _ in common_tangents
            ]
            tangents_on_x = [[common_tangents[0][0], y_min], [common_tangents[1][0], y_min]]
            x_points = [
                self.draw_tangent_point(_, axes)
                for _ in tangents_on_x
            ]

            connecting_vlines = [
                DashedLine(
                    start=g_x_point,
                    end=x_point,
                    stroke_width=2,
                    color=BLACK
                )
                for g_x_point, x_point in zip(tangent_points, x_points)
            ]
            common_tangent_line = self.get_line(common_tangents[0], common_tangents[1], axes)
            return VGroup(
                common_tangent_line,
                *tangent_points,
                *connecting_vlines,
                *x_points
            )

        tangent = always_redraw(get_tangent)

        self.play(
            Create(solid_gibbs)
        )
        self.play(
            Create(liquid_gibbs)
        )
        self.wait()

        self.play(
            LaggedStart(
                Write(mu_a_solid),
                Write(mu_a_liquid),
                Write(mu_b_solid),
                Write(mu_b_liquid),
                lag_ratio=1
            )
        )
        self.wait()

        self.play(
            Create(tangent[0])
        )
        self.play(
            FadeIn(tangent[1:3])
        )
        self.play(
            Create(tangent[3]),
            Create(tangent[4])
        )
        self.play(
            FadeIn(tangent[-2:])
        )
        self.add(tangent)
        self.wait()

        self.play(
            liquid_a_ref_tracker.animate.set_value(-delta_mu_a)
        )
        self.wait()

        self.play(
            liquid_a_ref_tracker.animate.set_value(0)
        )
        self.wait()

        self.play(
            solid_b_ref_tracker.animate.set_value(-delta_mu_b)
        )
        self.wait()

        self.play(
            liquid_a_ref_tracker.animate.set_value(-1200),
            solid_b_ref_tracker.animate.set_value(300),
            run_time=2, rate_func=linear
        )
        self.wait()

        self.play(
            liquid_a_ref_tracker.animate.set_value(500),
            solid_b_ref_tracker.animate.set_value(-1200),
            run_time=2, rate_func=linear
        )
        self.wait()

