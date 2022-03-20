import numpy as np
import pandas as pd
from itertools import product
from manim import *
from impingement import get_circ_intersec


class SpecialAxes(Axes):
    def __init__(self, origin, **kwargs):
        Axes.__init__(self, **kwargs)
        self.origin = origin

    def get_origin(self) -> np.ndarray:
        origin_x, origin_y = self.origin
        return self.coords_to_point(origin_x, origin_y)


class TestNewAxes(Scene):
    # origin shifting doesn't work for now, probably better to just shift the data
    def construct(self):
        axes = SpecialAxes(
            origin=(-1, -1),
            x_range=[-3, -1, 1],
            y_range=[-3, -1, 1]
        ).add_coordinates()

        self.add(axes)


# import the plotting data
saved_2d_spacetime = pd.read_csv("./data/2d_spacetime_results.csv")
saved_2d_contours = pd.read_csv("./data/2d_spacetime_contours.csv")

# get the contour levels
contour_levels = saved_2d_contours["level"].unique()

# shift the x, y values so that the smallest (x, y) pair is at (0, 0)
shift_by = 2
contour_lines = [
    (saved_2d_contours[saved_2d_contours["level"] == level_val][["x", "y"]] + shift_by).to_numpy()
    for level_val in contour_levels
]

# define the nucleation and growth rate ratios
nuc_rate_ratios = np.logspace(start=-2, stop=2, num=100)
growth_rate_ratios = np.logspace(start=-2, stop=2, num=100)

# get the probabilities when alpha wins and those when beta wins
alpha_wins_probs: np.ndarray = saved_2d_spacetime["alpha_wins_probs"].to_numpy()
beta_wins_probs: np.ndarray = saved_2d_spacetime["beta_wins_probs"].to_numpy()

growth_rate_ratios_3d = [
    growth_ratio
    for _, growth_ratio in product(nuc_rate_ratios, growth_rate_ratios)
]
nuc_rate_ratios_3d = [
    nuc_rate_ratio
    for nuc_rate_ratio, _ in product(nuc_rate_ratios, growth_rate_ratios)
]

prob_ratios = alpha_wins_probs / beta_wins_probs

prob_ratios_reshaped = prob_ratios.reshape((len(nuc_rate_ratios), len(growth_rate_ratios)))

nuc_rate_ratios_log10 = np.log10(nuc_rate_ratios)
growth_rate_ratios_log10 = np.log10(growth_rate_ratios)
prob_ratios_log10 = np.log10(prob_ratios_reshaped).T

prob_ratios_log10 += 6


class TestParametricSurface(ThreeDScene):

    @staticmethod
    def get_z_value(x, y):
        # return x + y
        return prob_ratios_log10[x, y]

    def construct(self):
        # define the axes
        x_max, z_max = 4, 13
        y_max = x_max

        axes_3d = ThreeDAxes(
            x_range=[0, x_max, 1],
            y_range=[0, y_max, 1],
            z_range=[0, z_max, 1],
            x_length=5,
            y_length=5,
            z_length=5,
            tips=False
        )
        x_axis = axes_3d.get_x_axis()
        y_axis = axes_3d.get_y_axis()
        z_axis = axes_3d.get_z_axis()

        # get the unit length for the axes
        unit_y_length = axes_3d.c2p(0, 1, 0)[1] - axes_3d.c2p(0, 0, 0)[1]

        x_label = axes_3d.get_x_axis_label(Tex("$\\log_{10}\\left(\\frac{\\dot{R_\\alpha}}{\\dot{R_\\beta}}\\right)$"))
        y_label = axes_3d.get_y_axis_label(
            Tex("$\\log_{10}\\left(\\frac{J_\\alpha}{J_\\beta}\\right)$"),
            edge=LEFT, direction=LEFT, buff=0.4
        )
        z_label = axes_3d.get_z_axis_label(Tex(
            "$\\log_{10}\\left(\\frac{P_\\alpha}{P_\\beta}\\right)$",
            font_size=30
        ))

        # add the quasi-2d axes
        self.add(x_axis, y_axis)
        self.play(
            Write(y_label)
        )
        self.wait()
        self.play(
            Write(x_label),
        )
        self.wait()

        # create the dots
        sampling_number = 100
        sampling_unit_increment = x_max / sampling_number
        dots_on_2d_plane = [Dot(point=axes_3d.c2p((row_index + 0.5) * sampling_unit_increment,
                                                  (col_index + 0.5) * sampling_unit_increment, 0),
                                radius=(axes_3d.c2p(0.2 * sampling_unit_increment, 0, 0) - axes_3d.c2p(0, 0, 0))[0])
                            for row_index in range(sampling_number) for col_index in range(sampling_number)
                            ]
        dots_on_2d_plane = VGroup(*dots_on_2d_plane)
        self.play(
            FadeIn(dots_on_2d_plane)
        )

        # re-orient the camera to move from quasi-2d to 3d
        self.move_camera(
            phi=80 * DEGREES, theta=-90 * DEGREES,
            frame_center=axes_3d.get_center()
        )

        # add the third axis and its label
        self.play(
            LaggedStart(
                *[
                    Create(z_axis),
                    Write(z_label)
                ],
                lag_ratio=1.2
            )
        )

        # shift the position of each dot
        shift_vectors = [
            axes_3d.c2p(
                (row_index + 0.5) * sampling_unit_increment,
                (col_index + 0.5) * sampling_unit_increment,
                self.get_z_value(row_index, col_index)
            ) - axes_3d.c2p(
                (row_index + 0.5) * sampling_unit_increment,
                (col_index + 0.5) * sampling_unit_increment,
                0
            )
            for row_index in range(sampling_number) for col_index in range(sampling_number)
        ]
        self.play(
            LaggedStart(
                *[dot.animate.shift(shift_vector)
                  for dot, shift_vector in zip(dots_on_2d_plane, shift_vectors)],
                lag_ratio=0.3
            ),
            run_time=5
        )
        self.wait()

        # add the surface
        surface = Surface(
            lambda u, v: axes_3d.c2p(
                u, v, self.get_z_value(
                    int(np.round(u / (4 / 99))),
                    int(np.round(v / (4 / 99)))
                    # u, v
                )
            ),
            resolution=(sampling_number - 1, sampling_number - 1),
            u_range=[0, x_max],
            v_range=[0, x_max],
            stroke_opacity=0.2,
            fill_opacity=0.5
        )
        surface.set_fill_by_value(
            axes=axes_3d, colors=[PURE_GREEN, "#f06b0c"], axis=2
        )

        self.play(
            Create(surface),
            FadeOut(dots_on_2d_plane)
        )
        self.wait()

        # add the contour lines
        contours = [
            axes_3d.plot_line_graph(
                x_values=contour_coords[:, 0],
                y_values=contour_coords[:, 1],
                add_vertex_dots=False
            ).set_color(
                interpolate_color(PURE_GREEN, "#f06b0c", _ / len(contour_lines))
            )
            for _, contour_coords in enumerate(contour_lines)
        ]

        # recreate the contours lines to be shifted upward
        contours_shifted = [
            axes_3d.plot_line_graph(
                x_values=contour_coords[:, 0],
                y_values=contour_coords[:, 1],
                add_vertex_dots=False
            ).set_color(
                interpolate_color(PURE_GREEN, "#f06b0c", _ / len(contour_lines))
            )
            for _, contour_coords in enumerate(contour_lines)
        ]
        contours_shifted = [
            contour.shift(
                axes_3d.c2p(0, 0, _) - axes_3d.c2p(0, 0, 0)
            )
            for _, contour in enumerate(contours_shifted)
        ]

        # add the slicing planes
        slicing_planes = [
            Square(side_length=y_max * unit_y_length).move_to(
                axes_3d.c2p(x_max / 2, y_max / 2, _)
            ).set_opacity(0.2).set_stroke(opacity=0.2).set_shade_in_3d()
            for _ in range(len(contours_shifted))
        ]

        prev_slice = contours[0]
        slice_shifted = contours_shifted[0]
        slicing_plane = slicing_planes[0]

        self.play(
            DrawBorderThenFill(slicing_plane)
        )
        self.play(
            Create(slice_shifted)
        )
        self.wait()
        self.play(
            TransformFromCopy(
                slice_shifted,
                prev_slice
            )
        )
        for contour, contour_shifted, slice_plane in zip(contours[1:], contours_shifted[1:],
                                                         slicing_planes[1:]):
            self.play(
                TransformFromCopy(
                    prev_slice,
                    contour
                ),
                Transform(
                    slice_shifted,
                    contour_shifted
                ),
                Transform(
                    slicing_plane,
                    slice_plane
                )
            )
            prev_slice = contour
        self.wait()

        # hide the surface and reorient the camera
        self.play(
            FadeOut(surface, slice_shifted, slicing_plane)
        )
        self.wait()

        # move back to starting quasi 2D view
        self.move_camera(
            phi=0 * DEGREES, theta=-90 * DEGREES,
            added_anims=[
                FadeOut(z_axis, z_label)
            ]
        )
        self.wait()

        # rotate the camera around z-axis
        # self.begin_ambient_camera_rotation(
        #     rate=0.2
        # )
        # self.wait(5)
        # self.stop_ambient_camera_rotation()
        # self.wait()


class FocusOn3D(Transform):
    """Shrink a spotlight to a position.

    Parameters
    ----------
    focus_point
        The point at which to shrink the spotlight. If it is a :class:`.~Mobject` its center will be used.
    opacity
        The opacity of the spotlight.
    color
        The color of the spotlight.
    run_time
        The duration of the animation.
    kwargs : Any
        Additional arguments to be passed to the :class:`~.Succession` constructor

    Examples
    --------
    .. manim:: UsingFocusOn

        class UsingFocusOn(Scene):
            def construct(self):
                dot = Dot(color=YELLOW).shift(DOWN)
                self.add(Tex("Focusing on the dot below:"), dot)
                self.play(FocusOn(dot))
                self.wait()
    """

    def __init__(
            self,
            focus_point,
            opacity: float = 0.2,
            color: str = GREY,
            run_time: float = 2,
            **kwargs
    ) -> None:
        self.focus_point = focus_point
        self.color = color
        self.opacity = opacity
        remover = True
        starting_dot = Dot(
            radius=config["frame_x_radius"] + config["frame_y_radius"],
            stroke_width=0,
            fill_color=self.color,
            fill_opacity=0,
        ).rotate(PI / 2, axis=X_AXIS)
        super().__init__(starting_dot, run_time=run_time, remover=remover, **kwargs)

    def create_target(self) -> Dot:
        little_dot = Dot(radius=0).rotate(PI / 2, axis=X_AXIS)
        little_dot.set_fill(self.color, opacity=self.opacity)
        little_dot.add_updater(lambda d: d.move_to(self.focus_point))
        return little_dot


# slide on re-framing the time zone calculation as a competition to the time axis
class TestCompetitivePhaseTransformation(ThreeDScene):
    def construct(self):
        x_max = 3
        y_max = 3
        z_max = 4.5
        axes = ThreeDAxes(
            x_range=[-x_max, x_max, 0.5],
            y_range=[-y_max, y_max, 0.5],
            z_range=[0, z_max, 0.5],
            x_length=10,
            y_length=10,
            z_length=6.75
        ).shift(-3 * Z_AXIS)

        # get the individual axis
        x_axis = axes.get_x_axis()
        y_axis = axes.get_y_axis()
        z_axis = axes.get_z_axis()

        axes.x_axis = x_axis.rotate(PI / 2, axis=X_AXIS)
        unit_y_length = axes.c2p(0, 1, 0)[1] - axes.c2p(0, 0, 0)[1]
        unit_z_length = axes.c2p(0, 0, 1)[2] - axes.c2p(0, 0, 0)[2]
        unit_y_vector = axes.c2p(0, 1, 0) - axes.c2p(0, 0, 0)
        unit_z_vector = axes.c2p(0, 0, 1) - axes.c2p(0, 0, 0)

        # add the axes labels
        x_label = axes.get_x_axis_label(Tex("$x$")).rotate(PI / 2, axis=X_AXIS)
        y_label = axes.get_y_axis_label(Tex("$y$")).rotate(PI / 2, axis=X_AXIS)
        z_label = axes.get_z_axis_label(Tex("$t$"))

        # create the quasi-2d view (1D SpaceTime)
        self.set_camera_orientation(
            100 * DEGREES, -90 * DEGREES,
        )
        self.play(
            LaggedStart(
                *[
                    Create(z_axis),
                    Write(z_label),
                    Create(x_axis),
                    Write(x_label)
                ],
                lag_ratio=1.1
            ),
            run_time=3
        )
        self.wait()

        # add the point of interest
        t_n = 4
        point_of_interest = Dot3D(
            point=axes.c2p(
                0, 0, t_n
            )
        )
        point_of_interest_label = Tex("$t_N$").move_to(axes.c2p(0.25, 0, t_n)).rotate(PI / 2, axis=X_AXIS)
        self.play(
            FadeIn(point_of_interest),
            FocusOn3D(point_of_interest, run_time=1.6)
        )
        self.wait()
        self.play(
            Write(point_of_interest_label)
        )
        self.wait()
        self.play(
            Unwrite(point_of_interest_label)
        )
        self.wait()

        # define the alpha phase growth rate and the beta phase growth rate
        alpha_rate = 0.5
        beta_rate = 0.25

        # define the inverted alpha and beta regions
        alpha_tot_region_left = DashedLine(
            start=axes.c2p(0, 0, t_n),
            end=axes.c2p(-alpha_rate * t_n, 0, 0),
            color=ORANGE
        )
        alpha_tot_region_right = DashedLine(
            start=axes.c2p(0, 0, t_n),
            end=axes.c2p(alpha_rate * t_n, 0, 0),
            color=ORANGE
        )
        self.play(
            Create(alpha_tot_region_left),
            Create(alpha_tot_region_right)
        )
        self.wait()

        beta_tot_region_left = DashedLine(
            start=axes.c2p(0, 0, t_n),
            end=axes.c2p(-beta_rate * t_n, 0, 0),
            color=GREEN
        )
        beta_tot_region_right = DashedLine(
            start=axes.c2p(0, 0, t_n),
            end=axes.c2p(beta_rate * t_n, 0, 0),
            color=GREEN
        )
        self.play(
            Create(beta_tot_region_left),
            Create(beta_tot_region_right)
        )
        self.wait()

        # highlight Region 1 and 2
        region_1_left = Polygon(
            axes.c2p(0, 0, t_n),
            axes.c2p(-alpha_rate * t_n, 0, 0),
            axes.c2p(-beta_rate * t_n, 0, 0),
            color=ORANGE,
            fill_opacity=0.7,
            stroke_width=0.3
        )
        region_1_right = Polygon(
            axes.c2p(0, 0, t_n),
            axes.c2p(alpha_rate * t_n, 0, 0),
            axes.c2p(beta_rate * t_n, 0, 0),
            color=ORANGE,
            fill_opacity=0.7,
            stroke_width=0.3
        )

        region_1_left_label = Tex("$1$").move_to(axes.c2p(-0.95, 0, 1.5)).rotate(PI / 2, axis=X_AXIS)
        region_1_right_label = Tex("$1$").move_to(axes.c2p(0.95, 0, 1.5)).rotate(PI / 2, axis=X_AXIS)

        self.play(
            DrawBorderThenFill(region_1_left),
            DrawBorderThenFill(region_1_right)
        )
        self.play(
            Write(region_1_left_label),
            Write(region_1_right_label)
        )
        self.wait()
        self.play(
            FadeOut(region_1_left, region_1_right, region_1_left_label, region_1_right_label)
        )
        self.wait()

        region_2 = Polygon(
            axes.c2p(0, 0, t_n),
            axes.c2p(-beta_rate * t_n, 0, 0),
            axes.c2p(beta_rate * t_n, 0, 0),
            color=GREEN,
            fill_opacity=0.7,
            stroke_width=0.3
        )
        region_2_label = Tex("$2$").move_to(axes.c2p(0, 0, 1.5)).rotate(PI / 2, axis=X_AXIS)
        self.play(
            DrawBorderThenFill(region_2)
        )
        self.play(
            Write(region_2_label)
        )
        self.wait()
        self.play(
            FadeOut(region_2, region_2_label)
        )
        self.wait()

        # darken the Region 1 and 2 boundaries
        self.play(
            alpha_tot_region_left.animate.set_stroke(opacity=0.5),
            alpha_tot_region_right.animate.set_stroke(opacity=0.5),
            beta_tot_region_left.animate.set_stroke(opacity=0.5),
            beta_tot_region_right.animate.set_stroke(opacity=0.5)
        )
        self.wait()

        # define the alpha and beta phase growth cone
        alpha_x_tracker = ValueTracker(1)
        alpha_y_tracker = ValueTracker(0)
        alpha_t_tracker = ValueTracker(1.5)

        cone_t_limit = 8 * unit_z_length

        alpha_cone_1d_left = always_redraw(
            lambda:
            Line(
                start=axes.c2p(alpha_x_tracker.get_value(), alpha_y_tracker.get_value(), alpha_t_tracker.get_value()),
                end=axes.c2p(alpha_x_tracker.get_value() - cone_t_limit * alpha_rate,
                             alpha_y_tracker.get_value(), cone_t_limit + alpha_t_tracker.get_value()),
                color=ORANGE
            )
        )
        alpha_cone_1d_right = always_redraw(
            lambda:
            Line(
                start=axes.c2p(alpha_x_tracker.get_value(), alpha_y_tracker.get_value(), alpha_t_tracker.get_value()),
                end=axes.c2p(alpha_x_tracker.get_value() + cone_t_limit * alpha_rate,
                             alpha_y_tracker.get_value(), cone_t_limit + alpha_t_tracker.get_value()),
                color=ORANGE
            )
        )

        beta_x_tracker = ValueTracker(-0.3)
        beta_y_tracker = ValueTracker(0)
        beta_t_tracker = ValueTracker(1)

        cone_t_limit = 8 * unit_z_length

        beta_cone_1d_left = always_redraw(
            lambda:
            Line(
                start=axes.c2p(beta_x_tracker.get_value(), beta_y_tracker.get_value(), beta_t_tracker.get_value()),
                end=axes.c2p(beta_x_tracker.get_value() - cone_t_limit * beta_rate,
                             beta_y_tracker.get_value(), cone_t_limit + beta_t_tracker.get_value()),
                color=GREEN
            )
        )
        beta_cone_1d_right = always_redraw(
            lambda:
            Line(
                start=axes.c2p(beta_x_tracker.get_value(), beta_y_tracker.get_value(), beta_t_tracker.get_value()),
                end=axes.c2p(beta_x_tracker.get_value() + cone_t_limit * beta_rate,
                             beta_y_tracker.get_value(), cone_t_limit + beta_t_tracker.get_value()),
                color=GREEN
            )
        )
        self.play(
            Create(beta_cone_1d_left),
            Create(beta_cone_1d_right),
            run_time=1.5
        )
        self.play(
            Create(alpha_cone_1d_left),
            Create(alpha_cone_1d_right),
            run_time=1.5
        )
        self.wait()

        # get the impingement line
        impingement_line = always_redraw(
            lambda: DashedLine(
                start=axes.c2p(((
                                        alpha_rate * beta_rate
                                ) / (
                                        alpha_rate + beta_rate
                                ) * (
                                        (
                                                alpha_x_tracker.get_value() / alpha_rate) + (
                                                beta_x_tracker.get_value() / beta_rate) +
                                        alpha_t_tracker.get_value() - beta_t_tracker.get_value()
                                )),
                               0,
                               ((
                                        beta_rate * beta_t_tracker.get_value() +
                                        alpha_rate * alpha_t_tracker.get_value() -
                                        beta_x_tracker.get_value() + alpha_x_tracker.get_value()
                                ) / (
                                        alpha_rate + beta_rate)),
                               ),
                end=axes.c2p(((
                                      alpha_rate * beta_rate
                              ) / (
                                      alpha_rate + beta_rate
                              ) * (
                                      (
                                              alpha_x_tracker.get_value() / alpha_rate) + (
                                              beta_x_tracker.get_value() / beta_rate) +
                                      alpha_t_tracker.get_value() - beta_t_tracker.get_value()
                              )),
                             0, cone_t_limit)
            ).set_stroke(color=(ORANGE if ((
                                                   alpha_rate * beta_rate
                                           ) / (
                                                   alpha_rate + beta_rate
                                           ) * (
                                                   (
                                                           alpha_x_tracker.get_value() / alpha_rate) + (
                                                           beta_x_tracker.get_value() / beta_rate) +
                                                   alpha_t_tracker.get_value() - beta_t_tracker.get_value()
                                           )) <= 0 else GREEN))
        )
        self.play(
            Create(impingement_line)
        )

        # move the alpha cone
        self.play(
            alpha_x_tracker.animate.set_value(0.4),
            alpha_t_tracker.animate.set_value(0.5),
            run_time=2
        )
        self.wait()

        # reframing the problem as the race to the time axis
        # highlight the time axis
        self.play(
            z_axis.animate.set_color(color=YELLOW),
            z_label.animate.set_color(color=YELLOW),
            run_time=0.5
        )
        self.wait()

        # add the alpha and beta dot on the time axis
        alpha_time_dot = always_redraw(
            lambda:
            Dot3D(
                point=axes.c2p(
                    0, 0, alpha_t_tracker.get_value() + np.abs(alpha_x_tracker.get_value() / alpha_rate)
                ),
                color=ORANGE
            ).set_shade_in_3d(False)
        )

        beta_time_dot = always_redraw(
            lambda:
            Dot3D(
                point=axes.c2p(
                    0, 0, beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)
                ),
                color=GREEN
            ).set_shade_in_3d(False)
        )

        self.play(
            FadeIn(alpha_time_dot)
        )
        self.play(
            FadeIn(beta_time_dot)
        )
        self.wait()

        # go back to the case when beta wins
        self.play(
            alpha_x_tracker.animate.set_value(1),
            alpha_t_tracker.animate.set_value(1.5),
            run_time=2
        )
        self.wait()

        # move the alpha cone again and fade out the impingement line
        self.play(
            alpha_x_tracker.animate.set_value(0.4),
            alpha_t_tracker.animate.set_value(0.5),
            run_time=2
        )
        self.wait()
        self.play(
            FadeOut(impingement_line)
        )
        self.wait()

        # add the time indicator for alpha_time_dot and beta_time_dot
        alpha_time_height_hline = always_redraw(
            lambda:
            Line(
                start=axes.c2p(-0.1, 0, alpha_t_tracker.get_value() + np.abs(alpha_x_tracker.get_value() / alpha_rate)),
                end=axes.c2p(-0.15, 0, alpha_t_tracker.get_value() + np.abs(alpha_x_tracker.get_value() / alpha_rate)),
                color=ORANGE
            )
        )
        alpha_time_height_vline = always_redraw(
            lambda:
            DashedLine(
                start=axes.c2p(-0.125, 0,
                               alpha_t_tracker.get_value() + np.abs(alpha_x_tracker.get_value() / alpha_rate)),
                end=axes.c2p(-0.125, 0, 0),
                color=ORANGE
            )
        )
        alpha_time_height_label = always_redraw(
            lambda:
            Tex("$t_a$", color=ORANGE).rotate(
                PI / 2, axis=X_AXIS
            ).move_to(
                axes.c2p(-0.3, 0,
                         (alpha_t_tracker.get_value() + np.abs(alpha_x_tracker.get_value() / alpha_rate)) / 2),
                aligned_edge=Z_AXIS
            )
        )

        self.play(
            Create(alpha_time_height_hline), run_time=0.5
        )
        self.play(
            Create(alpha_time_height_vline), run_time=0.7
        )
        self.play(
            Write(alpha_time_height_label)
        )
        self.wait()

        beta_time_height_hline = always_redraw(
            lambda:
            Line(
                start=axes.c2p(0.1, 0, beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)),
                end=axes.c2p(0.15, 0, beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)),
                color=GREEN
            )
        )
        beta_time_height_vline = always_redraw(
            lambda:
            DashedLine(
                start=axes.c2p(0.125, 0,
                               beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)),
                end=axes.c2p(0.125, 0, 0),
                color=GREEN
            )
        )
        beta_time_height_label = always_redraw(
            lambda:
            Tex("$t_b$", color=GREEN).rotate(
                PI / 2, axis=X_AXIS
            ).move_to(
                axes.c2p(0.3, 0,
                         (beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)) / 2),
                aligned_edge=Z_AXIS
            )
        )

        self.play(
            Create(beta_time_height_hline), run_time=0.5
        )
        self.play(
            Create(beta_time_height_vline), run_time=0.7
        )
        self.play(
            Write(beta_time_height_label)
        )
        self.wait()

        # let t_a = t_b
        self.play(
            alpha_x_tracker.animate.set_value(0.85),
            run_time=2
        )
        self.wait()

        # stop updating time_height objects
        alpha_time_height_label.clear_updaters()
        alpha_time_height_hline.clear_updaters()
        alpha_time_height_vline.clear_updaters()
        beta_time_height_label.clear_updaters()
        beta_time_height_hline.clear_updaters()
        beta_time_height_vline.clear_updaters()

        # write out the formula for finding the alpha wins boundary
        self.play(
            beta_time_height_label.animate.move_to(
                axes.c2p(-4, 0, 3)
            ),
            alpha_time_height_label.animate.set_opacity(0),
            alpha_time_height_hline.animate.set_opacity(0),
            alpha_time_height_vline.animate.set_opacity(0)
        )
        self.wait()

        t_b_formula = Tex(
            "$t_b$", " $=t_\\beta+\\frac{|x_\\beta|}{\\dot{R_\\beta}}$"
        ).move_to(
            beta_time_height_label, aligned_edge=-X_AXIS
        ).rotate(PI / 2, axis=X_AXIS)

        t_b_formula.submobjects[0].set_opacity(0)

        self.play(
            Write(t_b_formula)
        )
        self.wait()

        # show the t_beta segments
        t_beta_vline_1 = DashedLine(
            start=axes.c2p(beta_x_tracker.get_value(), 0, beta_t_tracker.get_value()),
            end=axes.c2p(beta_x_tracker.get_value(), 0, 0)
        )
        t_beta_vline_2 = DashedLine(
            start=axes.c2p(beta_x_tracker.get_value(), 0, beta_t_tracker.get_value()),
            end=axes.c2p(beta_x_tracker.get_value(), 0,
                         beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate))
        )
        t_beta_hline = DashedLine(
            start=axes.c2p(beta_x_tracker.get_value(), 0, beta_t_tracker.get_value()),
            end=axes.c2p(0, 0, beta_t_tracker.get_value())
        )
        self.play(
            Create(t_beta_vline_1)
        )
        self.wait()
        self.play(
            t_beta_vline_1.animate.set_opacity(0.5),
            Create(t_beta_vline_2)
        )
        self.wait()
        self.play(
            Create(t_beta_hline)
        )
        self.wait()
        self.play(
            FadeOut(t_beta_vline_1, t_beta_vline_2, t_beta_hline,
                    beta_time_height_hline, beta_time_height_vline)
        )
        self.wait()

        # show the t_alpha segments
        self.play(
            alpha_time_height_label.animate.set_opacity(1),
            alpha_time_height_hline.animate.set_opacity(1),
            alpha_time_height_vline.animate.set_opacity(1)
        )
        self.play(
            alpha_time_height_label.animate.move_to(
                axes.c2p(-3.95, 0, 2.2)
            )
        )
        self.wait()

        t_a_formula_1 = Tex(
            "$t_a$", " $=t_\\alpha+\\frac{|x_\\alpha|}{\\dot{R_\\alpha}}$"
        ).move_to(
            alpha_time_height_label, aligned_edge=-X_AXIS
        ).rotate(PI / 2, axis=X_AXIS)

        t_a_formula_1.submobjects[0].set_opacity(0)

        self.play(
            Write(t_a_formula_1)
        )
        self.wait()

        t_alpha_vline = DashedLine(
            start=axes.c2p(alpha_x_tracker.get_value(), 0, alpha_t_tracker.get_value()),
            end=axes.c2p(alpha_x_tracker.get_value(), 0, 0)
        )
        t_alpha_hline = DashedLine(
            start=axes.c2p(alpha_x_tracker.get_value(), 0, alpha_t_tracker.get_value()),
            end=axes.c2p(0, 0, alpha_t_tracker.get_value())
        )
        self.play(
            Create(t_alpha_vline)
        )
        self.wait()
        self.play(
            t_alpha_vline.animate.set_opacity(0.5),
            Create(t_alpha_hline)
        )
        self.wait()

        # drop the alpha subscript
        t_a_formula_2 = Tex(
            "$t_a$", " $=t+\\frac{|x|}{\\dot{R_\\alpha}}$"
        ).move_to(
            alpha_time_height_label, aligned_edge=-X_AXIS
        ).rotate(PI / 2, axis=X_AXIS)

        t_a_formula_2.submobjects[0].set_opacity(0)

        self.play(
            TransformMatchingShapes(t_a_formula_1, t_a_formula_2)
        )
        self.wait()

        # get the boundary formula
        t_bound_1d_formula = Tex(
            "$t=-\\frac{1}{\\dot{R_\\alpha}}|x|+\\left(t_\\beta+\\frac{|x_\\beta|}{\\dot{R_\\beta}}\\right)$",
            color=BLUE
        ).move_to(
            axes.c2p(-3.9, 0, 1.2), aligned_edge=-X_AXIS
        ).rotate(PI / 2, axis=X_AXIS)

        self.play(
            Write(t_bound_1d_formula)
        )
        self.wait()
        self.play(
            FadeOut(
                t_alpha_vline, t_alpha_hline, alpha_cone_1d_left, alpha_cone_1d_right,
                alpha_time_height_hline, alpha_time_height_vline, alpha_time_dot
            )
        )
        self.wait()
        alpha_cone_1d_left.suspend_updating()
        alpha_cone_1d_right.suspend_updating()

        # move the boundary formula
        self.play(
            FadeOut(alpha_time_height_label, beta_time_height_label, t_a_formula_2, t_b_formula),
            t_bound_1d_formula.animate.shift(
                axes.c2p(0.2, 0, 2.2) - axes.c2p(0, 0, 1.2)
            )
        )
        self.wait()

        # draw the alpha wins region in 1d spacetime
        alpha_wins_1d_no_intersect = Line(
            start=axes.c2p(
                0, 0, beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)
            ),
            end=axes.c2p(
                alpha_rate * (beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)),
                0, 0
            ),
            color=BLUE
        )
        alpha_wins_1d_with_intersect_1 = Line(
            start=axes.c2p(
                0, 0, beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)
            ),
            end=axes.c2p(
                - alpha_rate * (beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)),
                0, 0
            ),
            color=BLUE
        )

        self.play(
            Create(alpha_wins_1d_no_intersect)
        )
        self.wait()
        self.play(
            Create(alpha_wins_1d_with_intersect_1)
        )
        self.wait()

        alpha_wins_1d_with_intersect_2 = Line(
            start=axes.c2p(
                beta_x_tracker.get_value(), 0, beta_t_tracker.get_value()
            ),
            end=axes.c2p(
                beta_x_tracker.get_value() - alpha_rate * beta_t_tracker.get_value(),
                0, 0
            ),
            color=BLUE
        )
        self.play(
            ReplacementTransform(alpha_wins_1d_with_intersect_1, alpha_wins_1d_with_intersect_2),
            run_time=1.5
        )
        self.wait()

        # define the alpha wins region in 1d spacetime
        alpha_wins_1d_region_2 = always_redraw(
            lambda:
            Polygon(
                axes.c2p(beta_x_tracker.get_value(), 0, beta_t_tracker.get_value()),
                axes.c2p(0, 0, beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)),
                axes.c2p(
                    alpha_rate * (beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)),
                    0, 0),
                axes.c2p(beta_x_tracker.get_value() - alpha_rate * beta_t_tracker.get_value(), 0, 0),
                color=BLUE,
                fill_opacity=0.7
            )
        )
        self.play(
            DrawBorderThenFill(alpha_wins_1d_region_2),
            FadeOut(beta_time_dot)
        )
        self.remove(alpha_wins_1d_with_intersect_2, alpha_wins_1d_no_intersect)
        self.wait()

        # write out the ratio
        alpha_numerator = Tex(
            "$||\\Omega_\\alpha||$"
        ).move_to(
            axes.c2p(1.75, 0, 3)
        ).rotate(PI / 2, axis=X_AXIS)

        self.play(
            TransformFromCopy(alpha_wins_1d_region_2, alpha_numerator)
        )
        self.wait()

        # highlight the total time cone volume
        tot_time_cone_1d = Polygon(
            axes.c2p(0, 0, t_n),
            axes.c2p(alpha_rate * t_n, 0, 0),
            axes.c2p(-alpha_rate * t_n, 0, 0),
            color=ORANGE,
            fill_opacity=0.7
        )
        self.play(
            DrawBorderThenFill(tot_time_cone_1d)
        )
        self.wait()

        alpha_divider = Line(
            start=axes.c2p(1.3, 0, 2.8),
            end=axes.c2p(2.2, 0, 2.8)
        )
        self.play(
            Create(alpha_divider)
        )

        alpha_denominator = Tex(
            "$||\\Omega_{\\text{tot}}||$"
        ).move_to(
            axes.c2p(1.75, 0, 2.6), aligned_edge=Z_AXIS
        ).rotate(PI / 2, axis=X_AXIS)

        self.play(
            ReplacementTransform(tot_time_cone_1d, alpha_denominator)
        )
        self.wait()

        # add the double integral sign
        double_int = Tex(
            "$\\iint$"
        ).scale(1.5).move_to(
            axes.c2p(1.25, 0, 2.8), aligned_edge=X_AXIS
        ).rotate(PI / 2, axis=X_AXIS)
        int_variables = Tex(
            "$dx\\,dt$"
        ).move_to(
            axes.c2p(2.3, 0, 2.8), aligned_edge=-X_AXIS
        ).rotate(PI / 2, axis=X_AXIS)

        self.play(
            Write(double_int),
            Write(int_variables),
        )
        self.wait()

        # move the alpha wins in Region 2
        self.play(
            beta_x_tracker.animate.set_value(-0.2),
            beta_t_tracker.animate.set_value(2.2)
        )
        self.wait()
        self.play(
            beta_x_tracker.animate.set_value(0),
        )
        self.wait()
        self.play(
            beta_x_tracker.animate.set_value(-0.8),
            beta_t_tracker.animate.set_value(0)
        )
        self.wait()
        self.play(
            beta_x_tracker.animate.set_value(-0.3),
            beta_t_tracker.animate.set_value(1)
        )
        self.wait()

        self.play(
            FadeOut(double_int, int_variables, alpha_divider, alpha_numerator, alpha_denominator)
        )
        self.wait()

        # transition into real 3d space for 2D spacetime
        self.move_camera(
            90 * DEGREES, -100 * DEGREES,
        )
        self.play(
            LaggedStart(
                *[
                    Create(y_axis),
                    Write(y_label),
                ],
                lag_ratio=1.1
            ),
            run_time=2
        )
        self.wait()

        # add the alpha wins boundary
        # write out the distance
        generalize_to_2d = Tex(
            "$|x|$", " $=$", " $\\sqrt{x^2+y^2}$"
        ).move_to(
            axes.c2p(-3.9, 0, 1.5), aligned_edge=-X_AXIS
        ).rotate(PI / 2, axis=X_AXIS)

        generalize_to_3d = Tex(
            "$|x|$", " $=$", " $\\sqrt{x^2+y^2+z^2}$"
        ).move_to(
            axes.c2p(-3.9, 0, 1.5), aligned_edge=-X_AXIS
        ).rotate(PI / 2, axis=X_AXIS)

        self.play(
            Write(generalize_to_2d)
        )
        self.wait()
        self.play(
            TransformMatchingTex(generalize_to_2d, generalize_to_3d)
        )
        self.wait()

        # remove all the equations displayed
        self.play(
            FadeOut(generalize_to_2d, generalize_to_3d, t_bound_1d_formula)
        )

        # draw the extra allow region allowed due to faster growth speed
        alpha_wins_2d_with_intersect_1 = Line(
            start=axes.c2p(
                0, 0, beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)
            ),
            end=axes.c2p(
                - alpha_rate * (beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)),
                0, 0
            ),
            color=BLUE
        )
        self.play(
            Create(alpha_wins_2d_with_intersect_1)
        )
        self.wait()

        extra_allowed_region = Polygon(
            axes.c2p(beta_x_tracker.get_value(), 0, beta_t_tracker.get_value()),
            axes.c2p(beta_x_tracker.get_value() - beta_rate * (
                    (
                            (alpha_rate - beta_rate) * np.abs(beta_x_tracker.get_value())
                    ) / (beta_rate * (alpha_rate + beta_rate))),
                     0,
                     (alpha_rate - beta_rate) * np.abs(beta_x_tracker.get_value()) / (
                             beta_rate * (alpha_rate + beta_rate)) + beta_t_tracker.get_value()
                     ),
            axes.c2p(
                - alpha_rate * (beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)),
                0, 0
            ),
            axes.c2p(beta_x_tracker.get_value() - alpha_rate * beta_t_tracker.get_value(), 0, 0),
            color=YELLOW_D,
            fill_opacity=0.7
        )
        self.play(
            DrawBorderThenFill(extra_allowed_region)
        )
        self.wait()

        # show two growing circles
        t_tracker = ValueTracker(0)
        # move the alpha growth cone to the boundary
        alpha_x_tracker.set_value(
            - alpha_rate * (beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate))
        )
        alpha_t_tracker.set_value(0)
        # AHA: always resume updating after the new values have been set
        alpha_cone_1d_left.resume_updating()
        alpha_cone_1d_right.resume_updating()
        self.play(
            Create(alpha_cone_1d_left),
            Create(alpha_cone_1d_right)
        )
        self.wait()
        # alpha_cone_1d_left.suspend_updating()
        # alpha_cone_1d_right.suspend_updating()

        # add the plane
        plane = always_redraw(
            lambda:
            Square(side_length=5 * unit_y_length).move_to(axes.c2p(0, 0, 0)).shift(
                t_tracker.get_value() * unit_z_vector
            ).set_opacity(0.1).set_stroke(opacity=0.5)
        )
        # add the growing alpha circle
        circ_alpha = always_redraw(
            lambda:
            Circle(radius=axes.c2p(t_tracker.get_value() * alpha_rate, 0, 0)[0],
                   color=ORANGE, stroke_width=4).move_to(
                axes.c2p(
                    - alpha_rate * (beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)),
                    0, 0
                )).shift(
                t_tracker.get_value() * unit_z_vector
            )
        )
        # add the growing beta circle
        beta_cone_time = beta_t_tracker.get_value()
        circ_beta = always_redraw(
            lambda:
            Circle(radius=axes.c2p((t_tracker.get_value() - beta_cone_time) * beta_rate, 0, 0)[0],
                   color=GREEN, stroke_width=4).move_to(
                axes.c2p(
                    beta_x_tracker.get_value(),
                    0, 0
                )).shift(
                t_tracker.get_value() * unit_z_vector
            ) if t_tracker.get_value() > beta_cone_time else Circle(radius=0).move_to(
                axes.c2p(
                    beta_x_tracker.get_value(),
                    0, 0
                )).shift(
                t_tracker.get_value() * unit_z_vector
            )
        )

        self.add(circ_alpha, circ_beta)
        self.play(
            DrawBorderThenFill(plane)
        )
        self.wait()

        # move the camera to a higher altitute
        self.move_camera(
            70 * DEGREES, -100 * DEGREES,
            frame_center=axes.c2p(0, 0, 2),
            added_anims=[
                x_axis.animate.rotate(- PI / 2, axis=X_AXIS),
                x_label.animate.rotate(- PI / 2, axis=X_AXIS),
                y_label.animate.rotate(- PI / 2, axis=X_AXIS)
            ]
        )
        self.wait()

        self.play(
            t_tracker.animate.set_value(
                (alpha_rate - beta_rate) * np.abs(beta_x_tracker.get_value()) / (
                        beta_rate * (alpha_rate + beta_rate)) + beta_t_tracker.get_value()
            ),
            run_time=3,
            rate_func=linear
        )
        self.wait()

        # add the impingement points
        intersec_points = always_redraw(
            lambda: get_circ_intersec(
                - alpha_rate * (beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)), 0,
                beta_x_tracker.get_value(), 0,
                t_tracker.get_value() * alpha_rate, (t_tracker.get_value() - beta_cone_time) * beta_rate,
                axes, t_tracker.get_value(),
                radius=0.04, color=YELLOW
            )
        )
        intersec_point_1_path = TracedPath(intersec_points.submobjects[0].get_center, stroke_color=PURPLE,
                                           stroke_width=4)
        intersec_point_2_path = TracedPath(intersec_points.submobjects[1].get_center, stroke_color=PURPLE,
                                           stroke_width=4)
        self.play(
            FadeIn(intersec_points, intersec_point_1_path, intersec_point_2_path)
        )

        self.play(
            t_tracker.animate.set_value(
                beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)
            ),
            run_time=2,
            rate_func=linear
        )
        self.wait()

        # reset the time to zero
        # self.play(
        #     t_tracker.animate.set_value(
        #         0
        #     ),
        #     run_time=1
        # )
        # self.wait()

        # # fade out the alpha wins 1d region 2
        # self.play(
        #     FadeOut(alpha_wins_1d_region_2)
        # )
        # self.wait()
        #
        # top down view for circle-circle impingement
        # self.move_camera(
        #     0 * DEGREES, -90 * DEGREES,
        #     frame_center=axes.c2p(0, 0, 2),
        #     added_anims=[
        #         FadeOut(z_axis, z_label)
        #     ]
        # )
        # self.wait()
        #
        # self.play(
        #     t_tracker.animate.set_value(
        #         beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)
        #     ),
        #     run_time=3,
        #     rate_func=linear
        # )
        # self.wait()

        # # return to quasi-2D view
        # self.move_camera(
        #     100 * DEGREES, -90 * DEGREES,
        #     frame_center=np.array([0, 0, 0]),
        #     added_anims=[
        #         FadeIn(z_axis, z_label),
        #         FadeOut(y_axis, y_label),
        #         x_axis.animate.rotate(PI / 2, axis=X_AXIS),
        #         x_label.animate.rotate(PI / 2, axis=X_AXIS),
        #     ]
        # )
        # self.wait()
