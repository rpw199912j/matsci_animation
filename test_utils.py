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
saved_1d_spacetime = pd.read_csv("./data/1d_spacetime_results.csv")
saved_1d_contours = pd.read_csv("./data/1d_spacetime_contours.csv")
saved_2d_spacetime = pd.read_csv("./data/2d_spacetime_results.csv")
saved_2d_contours = pd.read_csv("./data/2d_spacetime_contours.csv")
saved_3d_spacetime = pd.read_csv("./data/3d_spacetime_results.csv")
saved_3d_contours = pd.read_csv("./data/3d_spacetime_contours.csv")

# get the contour levels
contour_levels_1d = saved_1d_contours["level"].unique()
contour_levels_2d = saved_2d_contours["level"].unique()
contour_levels_3d = saved_3d_contours["level"].unique()

# shift the x, y values so that the smallest (x, y) pair is at (0, 0)
shift_by = 2
contour_lines_1d = [
    (saved_1d_contours[saved_1d_contours["level"] == level_val][["x", "y"]] + shift_by).to_numpy()
    for level_val in contour_levels_1d
]
contour_lines_2d = [
    (saved_2d_contours[saved_2d_contours["level"] == level_val][["x", "y"]] + shift_by).to_numpy()
    for level_val in contour_levels_2d
]
contour_lines_3d = [
    (saved_3d_contours[saved_3d_contours["level"] == level_val][["x", "y"]] + shift_by).to_numpy()
    for level_val in contour_levels_3d
]

# define the nucleation and growth rate ratios
nuc_rate_ratios = np.logspace(start=-2, stop=2, num=100)
growth_rate_ratios = np.logspace(start=-2, stop=2, num=100)

# get the probabilities when alpha wins and those when beta wins
alpha_wins_probs_1d: np.ndarray = saved_1d_spacetime["alpha_wins_probs"].to_numpy()
beta_wins_probs_1d: np.ndarray = saved_1d_spacetime["beta_wins_probs"].to_numpy()
alpha_wins_probs_2d: np.ndarray = saved_2d_spacetime["alpha_wins_probs"].to_numpy()
beta_wins_probs_2d: np.ndarray = saved_2d_spacetime["beta_wins_probs"].to_numpy()
alpha_wins_probs_3d: np.ndarray = saved_3d_spacetime["alpha_wins_probs"].to_numpy()
beta_wins_probs_3d: np.ndarray = saved_3d_spacetime["beta_wins_probs"].to_numpy()

prob_ratios_1d = alpha_wins_probs_1d / beta_wins_probs_1d
prob_ratios_2d = alpha_wins_probs_2d / beta_wins_probs_2d
prob_ratios_3d = alpha_wins_probs_3d / beta_wins_probs_3d

prob_ratios_reshaped_1d = prob_ratios_1d.reshape((len(nuc_rate_ratios), len(growth_rate_ratios)))
prob_ratios_reshaped_2d = prob_ratios_2d.reshape((len(nuc_rate_ratios), len(growth_rate_ratios)))
prob_ratios_reshaped_3d = prob_ratios_3d.reshape((len(nuc_rate_ratios), len(growth_rate_ratios)))

nuc_rate_ratios_log10 = np.log10(nuc_rate_ratios)
growth_rate_ratios_log10 = np.log10(growth_rate_ratios)

prob_ratios_log10_1d = np.log10(prob_ratios_reshaped_1d).T
prob_ratios_log10_2d = np.log10(prob_ratios_reshaped_2d).T
prob_ratios_log10_3d = np.log10(prob_ratios_reshaped_3d).T

prob_ratios_log10_1d += 8
prob_ratios_log10_2d += 8
prob_ratios_log10_3d += 8


class TestParametricSurface(ThreeDScene):

    @staticmethod
    def get_z_value_1d(x, y):
        return prob_ratios_log10_1d[x, y]

    @staticmethod
    def get_z_value_2d(x, y):
        return prob_ratios_log10_2d[x, y]

    @staticmethod
    def get_z_value_3d(x, y):
        return prob_ratios_log10_3d[x, y]

    def construct(self):
        # define the axes
        x_max, z_max = 4, 17
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

        # add the custom labels
        values_x = [
            (_, tick_label) for _, tick_label in enumerate(["-2", "-1", "0", "1", "2"])
        ]
        values_z = [
            (_, str(tick_label))
            for _, tick_label in enumerate(np.arange(-8, 8 + 1))
        ]

        x_axis_labels = VGroup()
        y_axis_labels = VGroup()
        z_axis_labels = VGroup()

        for x_val, x_tex in values_x:
            tex = Tex(x_tex)  # Convert string to tex
            tex.next_to(x_axis.n2p(x_val), DOWN)  # Put tex on the position
            x_axis_labels.add(tex)

        for y_val, y_tex in values_x:
            tex = Tex(y_tex)  # Convert string to tex
            tex.next_to(y_axis.n2p(y_val), LEFT)  # Put tex on the position
            y_axis_labels.add(tex)

        for z_val, z_tex in values_z:
            tex = Tex(z_tex, font_size=20)  # Convert string to tex
            tex.next_to(z_axis.n2p(z_val), LEFT).rotate(PI / 2, axis=X_AXIS)  # Put tex on the position
            z_axis_labels.add(tex)

        # get the unit length for the axes
        unit_y_length = axes_3d.c2p(0, 1, 0)[1] - axes_3d.c2p(0, 0, 0)[1]

        x_label = axes_3d.get_x_axis_label(Tex("$\\log_{10}\\left(\\frac{\\dot{R_\\alpha}}{\\dot{R_\\beta}}\\right)$"))
        y_label = axes_3d.get_y_axis_label(
            Tex("$\\log_{10}\\left(\\frac{J_\\alpha}{J_\\beta}\\right)$"),
            edge=LEFT, direction=LEFT, buff=0.6
        )
        z_label = axes_3d.get_z_axis_label(Tex(
            "$\\log_{10}\\left(\\frac{P(\\alpha)}{P(\\beta)}\\right)$",
            font_size=30
        ))

        # add the quasi-2d axes
        self.add(x_axis, y_axis, x_axis_labels, y_axis_labels)
        self.play(
            Write(y_label)
        )
        self.wait()
        self.play(
            Write(x_label),
        )
        self.wait()

        # add the example point
        example_point = Dot(
            point=axes_3d.c2p(2, 2)
        )
        example_point_vline = DashedLine(
            start=axes_3d.c2p(2, 2),
            end=axes_3d.c2p(2, 0)
        )
        example_point_hline = DashedLine(
            start=axes_3d.c2p(2, 2),
            end=axes_3d.c2p(0, 2)
        )
        example_point_label = Tex(
            r"$\frac{P(\alpha)}{P(\beta)}$"
        ).move_to(axes_3d.c2p(2.5, 2.5, 0))

        self.play(
            FadeIn(example_point)
        )
        self.play(
            Create(example_point_vline),
            Create(example_point_hline),
            run_time=0.5
        )
        self.play(
            Write(example_point_label)
        )
        self.wait()

        self.play(
            FadeOut(example_point, example_point_vline, example_point_hline, example_point_label)
        )
        self.wait()

        # create the dots
        sampling_number = 100
        sampling_unit_increment = x_max / sampling_number
        # 1D Space-time
        dots_on_plane_1d = [Dot(point=axes_3d.c2p((row_index + 0.5) * sampling_unit_increment,
                                                  (col_index + 0.5) * sampling_unit_increment, 0),
                                radius=(axes_3d.c2p(0.5 * sampling_unit_increment, 0, 0) - axes_3d.c2p(0, 0, 0))[0]
                                ).set_color(interpolate_color(PURE_GREEN, ORANGE,
                                                              (self.get_z_value_1d(row_index, col_index) + np.abs(
                                                                  np.min(prob_ratios_log10_1d)
                                                              )) / (np.max(prob_ratios_log10_1d) + np.abs(
                                                                  np.min(prob_ratios_log10_1d)
                                                              ))))
                            for row_index in range(sampling_number) for col_index in range(sampling_number)
                            ]
        dots_on_plane_1d = VGroup(*dots_on_plane_1d)

        self.play(
            FadeIn(dots_on_plane_1d)
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
                    Write(z_label),
                    Write(z_axis_labels)
                ],
                lag_ratio=1.2
            )
        )

        # shift the position of each dot
        shift_vectors_1d = [
            axes_3d.c2p(
                (row_index + 0.5) * sampling_unit_increment,
                (col_index + 0.5) * sampling_unit_increment,
                self.get_z_value_1d(row_index, col_index)
            ) - axes_3d.c2p(
                (row_index + 0.5) * sampling_unit_increment,
                (col_index + 0.5) * sampling_unit_increment,
                0
            )
            for row_index in range(sampling_number) for col_index in range(sampling_number)
        ]
        shift_vectors_1d_to_2d = [
            axes_3d.c2p(
                (row_index + 0.5) * sampling_unit_increment,
                (col_index + 0.5) * sampling_unit_increment,
                self.get_z_value_2d(row_index, col_index)
            ) - axes_3d.c2p(
                (row_index + 0.5) * sampling_unit_increment,
                (col_index + 0.5) * sampling_unit_increment,
                self.get_z_value_1d(row_index, col_index)
            )
            for row_index in range(sampling_number) for col_index in range(sampling_number)
        ]
        shift_vectors_2d_to_3d = [
            axes_3d.c2p(
                (row_index + 0.5) * sampling_unit_increment,
                (col_index + 0.5) * sampling_unit_increment,
                self.get_z_value_3d(row_index, col_index)
            ) - axes_3d.c2p(
                (row_index + 0.5) * sampling_unit_increment,
                (col_index + 0.5) * sampling_unit_increment,
                self.get_z_value_2d(row_index, col_index)
            )
            for row_index in range(sampling_number) for col_index in range(sampling_number)
        ]
        shift_vectors_3d_to_1d = [
            axes_3d.c2p(
                (row_index + 0.5) * sampling_unit_increment,
                (col_index + 0.5) * sampling_unit_increment,
                self.get_z_value_1d(row_index, col_index)
            ) - axes_3d.c2p(
                (row_index + 0.5) * sampling_unit_increment,
                (col_index + 0.5) * sampling_unit_increment,
                self.get_z_value_3d(row_index, col_index)
            )
            for row_index in range(sampling_number) for col_index in range(sampling_number)
        ]
        self.play(
            LaggedStart(
                *[dot.animate.shift(shift_vector)
                  for dot, shift_vector in zip(dots_on_plane_1d, shift_vectors_1d)],
                lag_ratio=0.3
            ),
            run_time=5
        )
        self.wait()

        # move the camera
        self.move_camera(
            phi=80 * DEGREES, theta=-120 * DEGREES,
            frame_center=axes_3d.get_center()
        )
        # SECTION HERE

        # transition into 2d
        self.play(
            *[dot.animate.shift(shift_vector) for dot, shift_vector in zip(dots_on_plane_1d,
                                                                           shift_vectors_1d_to_2d)],
            run_time=1,
        )
        self.wait()
        # SECTION HERE
        # transition into 3d
        self.play(
            *[dot.animate.shift(shift_vector) for dot, shift_vector in zip(dots_on_plane_1d,
                                                                           shift_vectors_2d_to_3d)],
            run_time=1,
        )
        self.wait()
        # SECTION HERE

        # # SECTION HERE
        # rotate the camera !!!
        self.begin_ambient_camera_rotation(
            rate=-33 * DEGREES,
        )
        self.wait(10)
        self.stop_ambient_camera_rotation()
        self.set_camera_orientation(
            phi=80 * DEGREES, theta=-90 * DEGREES,
            frame_center=axes_3d.get_center()
        )
        # SECTION HERE

        # transition back to 1d
        self.play(
            *[dot.animate.shift(shift_vector) for dot, shift_vector in zip(dots_on_plane_1d,
                                                                           shift_vectors_3d_to_1d)],
            run_time=1,
        )
        self.wait()

        # add the surface
        surface = Surface(
            lambda u, v: axes_3d.c2p(
                u, v, self.get_z_value_1d(int(np.round(u / (4 / 99))), int(np.round(v / (4 / 99))))
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
            FadeOut(dots_on_plane_1d)
        )
        self.wait()

        # add the contour lines
        contours_1d = [
            axes_3d.plot_line_graph(
                x_values=contour_coords[:, 0],
                y_values=contour_coords[:, 1],
                add_vertex_dots=False
            ).set_color(
                interpolate_color(PURE_GREEN, "#f06b0c", _ / len(contour_lines_1d))
            )
            for _, contour_coords in enumerate(contour_lines_1d)
        ]
        # add the contour labels
        contours_1d_labels = [
            Tex(
                f"{_ - 4}", font_size=25
            ).add_background_rectangle().move_to(
                axes_3d.c2p(
                    contour_coords[contour_coords.shape[0] // 2, 0],
                    contour_coords[contour_coords.shape[0] // 2, 1],
                    0
                )
            )
            for _, contour_coords in enumerate(contour_lines_1d)
        ]
        contours_2d_labels = [
            Tex(
                f"{_ - 6}", font_size=25
            ).add_background_rectangle().move_to(
                axes_3d.c2p(
                    contour_coords[contour_coords.shape[0] // 2, 0],
                    contour_coords[contour_coords.shape[0] // 2, 1],
                    0
                )
            )
            for _, contour_coords in enumerate(contour_lines_2d)
        ]
        contours_3d_labels = [
            Tex(
                f"{_ - 8}", font_size=25
            ).add_background_rectangle().move_to(
                axes_3d.c2p(
                    contour_coords[contour_coords.shape[0] // 2, 0],
                    contour_coords[contour_coords.shape[0] // 2, 1],
                    0
                )
            )
            for _, contour_coords in enumerate(contour_lines_3d)
        ]

        # recreate the contours_1d lines to be shifted upward
        contours_shifted_1d = [
            axes_3d.plot_line_graph(
                x_values=contour_coords[:, 0],
                y_values=contour_coords[:, 1],
                add_vertex_dots=False
            ).set_color(
                interpolate_color(PURE_GREEN, "#f06b0c", _ / len(contour_lines_1d))
            )
            for _, contour_coords in enumerate(contour_lines_1d)
        ]
        contours_shifted_1d = [
            contour.shift(
                axes_3d.c2p(0, 0, _ + 4) - axes_3d.c2p(0, 0, 0)
            )
            for _, contour in enumerate(contours_shifted_1d)
        ]

        # add the slicing planes
        slicing_planes = [
            Square(side_length=y_max * unit_y_length).move_to(
                axes_3d.c2p(x_max / 2, y_max / 2, _ + 4)
            ).set_opacity(0.2).set_stroke(opacity=0.2).set_shade_in_3d()
            for _ in range(len(contours_shifted_1d))
        ]

        prev_slice = contours_1d[0]
        slice_shifted = contours_shifted_1d[0]
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
        for contour, contour_shifted, slice_plane in zip(contours_1d[1:], contours_shifted_1d[1:],
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
                FadeOut(z_axis, z_label, z_axis_labels)
            ]
        )
        self.wait()

        self.play(
            FadeIn(*contours_1d_labels)
        )
        self.wait()

        # add two horizontal lines
        hline_1 = DashedLine(
            start=axes_3d.c2p(0, 2.5, 0),
            end=axes_3d.c2p(4, 2.5, 0),
        )
        hline_2 = DashedLine(
            start=axes_3d.c2p(0, 1.5, 0),
            end=axes_3d.c2p(4, 1.5, 0),
        )
        v_line_1d = DashedLine(
            start=axes_3d.c2p(1.55, 2.5, 0),
            end=axes_3d.c2p(1.55, 1.5, 0),
        )


        self.play(
            Create(hline_1),
            Create(hline_2)
        )
        # SECTION HERE
        self.play(
            Create(v_line_1d)
        )
        # SECTION HERE
        self.play(
            FadeOut(v_line_1d)
        )
        # SECTION HERE

        self.play(
            FadeOut(*contours_1d_labels)
        )
        self.wait()

        # SECTION HERE

        # transform into the countour lines in 2d and 3d space-time
        contours_2d = [
            axes_3d.plot_line_graph(
                x_values=contour_coords[:, 0],
                y_values=contour_coords[:, 1],
                add_vertex_dots=False
            ).set_color(
                interpolate_color(PURE_GREEN, "#f06b0c", _ / len(contour_lines_2d))
            )
            for _, contour_coords in enumerate(contour_lines_2d)
        ]

        contours_3d = [
            axes_3d.plot_line_graph(
                x_values=contour_coords[:, 0],
                y_values=contour_coords[:, 1],
                add_vertex_dots=False
            ).set_color(
                interpolate_color(PURE_GREEN, "#f06b0c", _ / len(contour_lines_3d))
            )
            for _, contour_coords in enumerate(contour_lines_3d)
        ]

        self.play(
            *[
                ReplacementTransform(contour_1d, countour_2d)
                for contour_1d, countour_2d in zip(contours_1d, contours_2d[:len(contours_1d)])
            ],
            *[
                FadeIn(countour_2d)
                for countour_2d in contours_2d[len(contours_1d):]
            ]
        )
        self.wait()

        self.play(
            FadeIn(*contours_2d_labels)
        )
        self.wait()

        v_line_2d = DashedLine(
            start=axes_3d.c2p(1.85, 2.5, 0),
            end=axes_3d.c2p(1.85, 1.5, 0),
        )
        # SECTION HERE
        self.play(
            Create(v_line_2d)
        )
        # SECTION HERE
        self.play(
            FadeOut(v_line_2d)
        )

        self.play(
            FadeOut(*contours_2d_labels)
        )
        self.wait()
        # SECTION HERE

        self.play(
            *[
                ReplacementTransform(contour_2d, countour_3d)
                for contour_2d, countour_3d in zip(contours_2d, contours_3d[:len(contours_2d)])
            ],
            *[
                FadeIn(countour_3d)
                for countour_3d in contours_3d[len(contours_2d):]
            ]
        )
        self.wait()
        self.play(
            FadeIn(*contours_3d_labels)
        )
        self.wait()
        v_line_3d = DashedLine(
            start=axes_3d.c2p(1.9, 2.5, 0),
            end=axes_3d.c2p(1.9, 1.5, 0),
        )
        # SECTION HERE
        self.play(
            Create(v_line_3d)
        )
        # SECTION HERE
        self.wait()

        # add the observation statement
        observation = Tex(
            r"$\frac{P(\alpha)}{P(\beta)}$", r"$\propto$",
            "$\\begin{cases} \\frac{J_\\alpha\\cdot \\dot{R_\\alpha}}{J_\\beta\\cdot \\dot{R_\\beta}} & \\text{1D}\\\\ \\frac{J_\\alpha\\cdot {\\dot{R_\\alpha}}^2}{J_\\beta\\cdot {\\dot{R_\\beta}}^2} & \\text{2D}\\\\ \\frac{J_\\alpha\\cdot {\\dot{R_\\alpha}}^3}{J_\\beta\\cdot {\\dot{R_\\beta}}^3} & \\text{3D}\\end{cases}$",
        ).scale(0.85).to_corner(UR)

        self.play(
            Write(observation)
        )
        self.wait()


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
        alpha_cone_1d_left.suspend_updating()
        alpha_cone_1d_right.suspend_updating()

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
                radius=0.06, color=YELLOW
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

        # reset the time to the time of impingement
        self.play(
            t_tracker.animate.set_value(
                (alpha_rate - beta_rate) * np.abs(beta_x_tracker.get_value()) / (
                        beta_rate * (alpha_rate + beta_rate)) + beta_t_tracker.get_value()
            ),
            FadeOut(intersec_point_1_path, intersec_point_2_path),
            run_time=1
        )
        self.wait()

        # redefine the traced paths
        intersec_point_1_path = TracedPath(intersec_points.submobjects[0].get_center, stroke_color=PURPLE,
                                           stroke_width=4)
        intersec_point_2_path = TracedPath(intersec_points.submobjects[1].get_center, stroke_color=PURPLE,
                                           stroke_width=4)
        self.play(
            FadeIn(intersec_point_1_path, intersec_point_2_path),
            run_time=0.5
        )

        # # fade out the alpha wins 1d region 2
        # self.play(
        #     FadeOut(alpha_wins_1d_region_2)
        # )
        # self.wait()
        #
        # top down view for circle-circle impingement
        self.move_camera(
            0 * DEGREES, -90 * DEGREES,
            frame_center=axes.c2p(0, 0, 2),
            added_anims=[
                FadeOut(z_axis, z_label)
            ]
        )
        self.wait()

        self.play(
            t_tracker.animate.set_value(
                beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)
            ),
            run_time=3,
            rate_func=linear
        )
        self.wait()

        # fadeout all the unused mobjects
        self.play(
            FadeOut(circ_alpha),
            FadeOut(
                intersec_points, intersec_point_1_path, intersec_point_2_path
            )
        )
        self.wait()
        intersec_points.suspend_updating()
        intersec_point_1_path.suspend_updating()
        intersec_point_2_path.suspend_updating()
        circ_alpha.suspend_updating()

        # move the camera back to the 3D view with a higher altitute
        self.move_camera(
            70 * DEGREES, -100 * DEGREES,
            added_anims=[
                FadeIn(z_axis, z_label)
            ]
        )
        self.wait()

        # VISUALIZE THE CONES AND INTERSECTION OF CIRCLES
        # set the time tracker back to 0
        self.play(
            t_tracker.animate.set_value(0),
            run_time=3,
            rate_func=linear
        )
        self.wait()

        alpha_cone_radius = (axes.c2p(
            alpha_rate * (beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)), 0, 0
        ) - axes.c2p(0, 0, 0))[0]
        alpha_cone_height = (axes.c2p(
            0, 0, beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)
        ) - axes.c2p(0, 0, 0))[2]
        alpha_win_cone = Cone(
            base_radius=alpha_cone_radius,
            height=alpha_cone_height,
            stroke_color=BLUE,
            fill_color=BLUE,
            fill_opacity=0.3,
            stroke_opacity=0,
            resolution=50
        ).move_to(
            axes.c2p(0, 0, 0), aligned_edge=-Z_AXIS
        )

        beta_cone_radius = (axes.c2p(
            np.abs(beta_x_tracker.get_value()), 0, 0
        ) - axes.c2p(
            0, 0, 0))[0]
        beta_cone_height = (axes.c2p(
            0, 0, beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)
        ) - axes.c2p(
            0, 0, beta_t_tracker.get_value()))[2]
        beta_cone = Cone(
            base_radius=beta_cone_radius,
            height=beta_cone_height,
            stroke_color=GREEN_E,
            fill_color=GREEN_E,
            fill_opacity=0.2,
            stroke_opacity=0,
            show_base=True,
            direction=-Z_AXIS,
            resolution=50
        ).move_to(
            axes.c2p(
                beta_x_tracker.get_value(), 0, beta_t_tracker.get_value()
            ), aligned_edge=-Z_AXIS
        )

        self.play(
            Create(alpha_win_cone),
            Create(beta_cone)
        )
        self.play(
            FadeOut(alpha_wins_1d_region_2, extra_allowed_region)
        )
        self.wait()

        alpha_wins_circle = always_redraw(
            lambda:
            Circle(radius=axes.c2p(
                (
                        alpha_rate * (beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate))
                ) - t_tracker.get_value() * alpha_rate,
                0, 0)[0],
                   color=BLUE, stroke_width=4).move_to(
                axes.c2p(0, 0, 0)).shift(
                t_tracker.get_value() * unit_z_vector
            )
        )

        circ_overlap = always_redraw(
            lambda:
            Intersection(alpha_wins_circle, circ_beta,
                         color=PURPLE, fill_color=YELLOW, fill_opacity=1, stroke_width=2).shift(
                axes.c2p(0, 0, t_tracker.get_value()) - np.array([0, 0, 0])
            )
        )

        self.play(
            FadeIn(alpha_wins_circle, circ_overlap)
        )

        # move the time to t_beta
        self.play(
            t_tracker.animate.set_value(beta_t_tracker.get_value()),
            run_time=2,
            rate_func=linear
        )
        self.wait()
        # move the time to impingement time
        self.play(
            t_tracker.animate.set_value(
                (alpha_rate - beta_rate) * np.abs(beta_x_tracker.get_value()) / (
                        beta_rate * (alpha_rate + beta_rate)) + beta_t_tracker.get_value()
            ),
            run_time=2,
            rate_func=linear
        )
        self.wait()
        # move the time to impingement time
        self.play(
            t_tracker.animate.set_value(
                beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)
            ),
            run_time=2,
            rate_func=linear
        )
        self.wait()

        # fade out all the unused mobjects
        self.play(
            FadeOut(circ_overlap)
        )
        circ_overlap.clear_updaters()

        # SHOW THE 3D SPACETIME VERSION
        # move the time back to zero
        self.play(
            t_tracker.animate.set_value(0),
            run_time=2,
            rate_func=linear
        )
        self.wait()

        # define the spheres
        alpha_wins_sphere = always_redraw(
            lambda:
            Sphere(
                center=axes.c2p(0, 0, t_tracker.get_value()),
                radius=axes.c2p(
                    (
                            alpha_rate * (beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate))
                    ) - t_tracker.get_value() * alpha_rate,
                    0, 0)[0],
                stroke_width=1,
            ).set_color(BLUE).set_opacity(0.5)
        )

        sphere_beta = always_redraw(
            lambda:
            Sphere(
                center=axes.c2p(beta_x_tracker.get_value(), 0, t_tracker.get_value()),
                radius=axes.c2p((t_tracker.get_value() - beta_cone_time) * beta_rate, 0, 0)[0],
            ).set_color(GREEN_E).set_opacity(0.5)
            if t_tracker.get_value() > beta_cone_time else Sphere(
                center=axes.c2p(beta_x_tracker.get_value(), 0, t_tracker.get_value()),
                radius=0,
                stroke_opacity=0, fill_opacity=0
            )
        )

        self.play(
            Create(alpha_wins_sphere),
            Create(sphere_beta)
        )
        self.play(
            t_tracker.animate.set_value(
                beta_t_tracker.get_value() + np.abs(beta_x_tracker.get_value() / beta_rate)
            ),
            run_time=2,
            rate_func=linear
        )
        self.wait()

        # fade out unused mobjects
        self.play(
            FadeOut(
                plane, sphere_beta, alpha_wins_sphere,
                alpha_wins_circle, circ_beta, alpha_win_cone, beta_cone,
                alpha_wins_2d_with_intersect_1
            )
        )
        self.wait()

        # return to quasi-2D view
        self.move_camera(
            100 * DEGREES, -90 * DEGREES,
            frame_center=np.array([0, 0, 0]),
            added_anims=[
                FadeOut(y_axis, y_label, alpha_cone_1d_left, alpha_cone_1d_right),
                x_axis.animate.rotate(PI / 2, axis=X_AXIS),
                x_label.animate.rotate(PI / 2, axis=X_AXIS)
            ]
        )
        self.wait()

        # ALPHA WINS WHEN BETA IN REGION 1
        # move the beta cone position and time
        self.play(
            beta_x_tracker.animate.set_value(-1.1),
            beta_t_tracker.animate.set_value(0.75)
        )
        self.wait()

        # define the alpha win boundary
        alpha_wins_1d_region_1 = Polygon(
            axes.c2p(beta_x_tracker.get_value(), 0, beta_t_tracker.get_value()),
            axes.c2p(beta_x_tracker.get_value() - alpha_rate * beta_t_tracker.get_value(), 0, 0),
            axes.c2p(t_n * alpha_rate, 0, 0),
            axes.c2p(0, 0, t_n),
            axes.c2p(
                beta_x_tracker.get_value() + beta_rate * (
                        (alpha_rate * t_n - beta_rate * beta_t_tracker.get_value() - np.abs(
                            beta_x_tracker.get_value())) / (
                                alpha_rate - beta_rate
                        ) - beta_t_tracker.get_value()
                ),
                0,
                (alpha_rate * t_n - beta_rate * beta_t_tracker.get_value() - np.abs(beta_x_tracker.get_value())) / (
                        alpha_rate - beta_rate
                )
            ),
            color=BLUE, fill_opacity=0.7
        )
        self.play(
            DrawBorderThenFill(alpha_wins_1d_region_1)
        )
        self.wait()

        self.play(
            FadeOut(alpha_wins_1d_region_1),
            FadeOut(beta_cone_1d_left, beta_cone_1d_right)
        )
        self.wait()
        beta_cone_1d_left.suspend_updating()
        beta_cone_1d_right.suspend_updating()

        # BETA WINS WHEN ALPHA IN REGION 2
        # move the alpha cone into Region 2
        alpha_x_tracker.set_value(-0.3)
        alpha_t_tracker.set_value(0.8)
        alpha_cone_1d_left.resume_updating()
        alpha_cone_1d_right.resume_updating()
        self.play(
            Create(alpha_cone_1d_left),
            Create(alpha_cone_1d_right)
        )
        self.wait()

        # create the beta phase win polygon
        beta_wins_1d_region_2 = always_redraw(
            lambda:
            Polygon(
                axes.c2p(0, 0, alpha_t_tracker.get_value() + np.abs(alpha_x_tracker.get_value() / alpha_rate)),
                axes.c2p(
                    beta_rate * (alpha_t_tracker.get_value() + np.abs(alpha_x_tracker.get_value() / alpha_rate)),
                    0, 0),
                axes.c2p(
                    alpha_x_tracker.get_value() - beta_rate * alpha_t_tracker.get_value(),
                    0, 0),
                axes.c2p(
                    alpha_x_tracker.get_value(),
                    0, alpha_t_tracker.get_value()),
                color=BLUE, fill_opacity=0.7
            )
        )
        self.play(
            DrawBorderThenFill(beta_wins_1d_region_2)
        )
        self.wait()

        # BETA WINS WHEN ALPHA IN REGION 1
        # move the alpha cone into Region 1
        self.play(
            alpha_x_tracker.animate.set_value(-1.1),
            run_time=2
        )
        self.wait()

        # highlight the forbidden region
        nothing_wins_region = Polygon(
            axes.c2p(
                alpha_x_tracker.get_value(),
                0, alpha_t_tracker.get_value()),
            axes.c2p(
                alpha_x_tracker.get_value() - beta_rate * alpha_t_tracker.get_value(),
                0, 0),
            axes.c2p(-beta_rate * t_n, 0, 0),
            axes.c2p(
                alpha_x_tracker.get_value() + alpha_rate * (
                        (alpha_rate * alpha_t_tracker.get_value() - beta_rate * t_n + np.abs(
                            alpha_x_tracker.get_value())) / (
                                alpha_rate - beta_rate
                        ) - alpha_t_tracker.get_value()
                ),
                0,
                (alpha_rate * alpha_t_tracker.get_value() - beta_rate * t_n + np.abs(alpha_x_tracker.get_value())) / (
                        alpha_rate - beta_rate
                )
            ),
            color=YELLOW_D, fill_opacity=0.4
        )
        self.play(
            DrawBorderThenFill(nothing_wins_region)
        )
        self.wait()

        # show the beta growth cone
        beta_x_tracker.set_value(-0.95)
        beta_t_tracker.set_value(0.7)

        beta_cone_1d_left.resume_updating()
        beta_cone_1d_right.resume_updating()
        self.play(
            Create(beta_cone_1d_left),
            Create(beta_cone_1d_right)
        )
        self.wait()


class TestTimeConeVolumeSummary(Scene):
    def construct(self):
        summary_table_1d_str = r"""
        \begin{table}[h]
    \centering
    \begin{tabular}{|c|c|}
        \hline
        $||\Omega||$& 1D Space-time\\
        \hline
        $\alpha$ wins & $\frac{1}{6}{t_N}^2\left(5\dot{R_\alpha}-2\dot{R_\beta}\right)$\\
        \hline
        $\beta$ wins & $\frac{1}{6}{t_N}^2\frac{4\dot{R_\alpha}\dot{R_\beta}-{\dot{R_\beta}}^2}{\dot{R_\alpha}}$\\
        \hline
        Nothing wins&$\frac{1}{6}{t_N}^2\frac{\left(\dot{R_\alpha}-\dot{R_\beta}\right)^2}{\dot{R_\alpha}}$\\
        \hline
        Sum &$\dot{R_\alpha}{t_N}^2$\\
        \hline
    \end{tabular}
\end{table}
        """
        summary_table_1d = Tex(summary_table_1d_str)
        self.play(
            Write(summary_table_1d)
        )
        self.wait()
        self.play(
            FadeOut(summary_table_1d)
        )
        self.wait()

        summary_table_3d_str = r"""
        \begin{table}[h]
    \centering
    \begin{tabular}{|c|c|}
        \hline
        $||\Omega||$& 3D Space-time\\
        \hline
        $\alpha$ wins & $\frac{\pi  {t_N}^4 \left(35 {\dot{R_\alpha}}^7+105 {\dot{R_\alpha}}^6 {\dot{R_\beta}}+105 {\dot{R_\alpha}}^5 {\dot{R_\beta}}^2+17 {\dot{R_\alpha}}^4
   {\dot{R_\beta}}^3-54 {\dot{R_\alpha}}^3 {\dot{R_\beta}}^4-54 {\dot{R_\alpha}}^2 {\dot{R_\beta}}^5-8 {\dot{R_\alpha}} {\dot{R_\beta}}^6-6 {\dot{R_\beta}}^7\right)}{105 {\dot{R_\alpha}}
   ({\dot{R_\alpha}}+{\dot{R_\beta}})^3}$ \\
        \hline
        $\beta$ wins & $\frac{\pi  {\dot{R_\beta}}^3 {t_N}^4 \left(2835 {\dot{R_\alpha}}^7-3780 {\dot{R_\alpha}}^6 {\dot{R_\beta}}+1890 {\dot{R_\alpha}}^5 {\dot{R_\beta}}^2-375
   {\dot{R_\alpha}}^4 {\dot{R_\beta}}^3-65 {\dot{R_\alpha}}^3 {\dot{R_\beta}}^4+66 {\dot{R_\alpha}}^2 {\dot{R_\beta}}^5-12 {\dot{R_\alpha}} {\dot{R_\beta}}^6+{\dot{R_\beta}}^7\right)}{210 {\dot{R_\alpha}}^3
   ({\dot{R_\beta}}-3 {\dot{R_\alpha}})^4}$\\
        \hline
        Nothing wins& $\frac{\pi  {\dot{R_\beta}}^3 {t_N}^4 ({\dot{R_\alpha}}-{\dot{R_\beta}})^3 \left(81 {\dot{R_\alpha}}^7+378 {\dot{R_\alpha}}^6 {\dot{R_\beta}}+864 {\dot{R_\alpha}}^5
   {\dot{R_\beta}}^2-219 {\dot{R_\alpha}}^4 {\dot{R_\beta}}^3+245 {\dot{R_\alpha}}^3 {\dot{R_\beta}}^4-6 {\dot{R_\alpha}} {\dot{R_\beta}}^6+{\dot{R_\beta}}^7\right)}{210 {\dot{R_\alpha}}^3
   ({\dot{R_\beta}}-3 {\dot{R_\alpha}})^4 ({\dot{R_\alpha}}+{\dot{R_\beta}})^3}$\\
        \hline
        Sum & $\frac{\pi}{3}{\dot{R_\alpha}}^3{t_N}^4$\\
        \hline
    \end{tabular}
\end{table}
        """
        summary_table_3d = Tex(summary_table_3d_str, font_size=30)
        self.play(
            Write(summary_table_3d)
        )
        self.wait()

        self.play(
            FadeOut(summary_table_3d)
        )
        self.wait()

        summary_table_2d_str = r"""
                \begin{table}[h]
            \centering
            \begin{tabular}{|c|c|}
                \hline
                $||\Omega||$& 2D Space-time\\
                \hline
                $\alpha$ wins & $f(\dot{R_\alpha}, \dot{R_\beta}, t_N)$ \\
                \hline
                $\beta$ wins & $f(\dot{R_\alpha}, \dot{R_\beta}, t_N)$\\
                \hline
                Nothing wins& $f(\dot{R_\alpha}, \dot{R_\beta}, t_N)$\\
                \hline
                Sum & $\frac{\pi}{3}{\dot{R_\alpha}}^2{t_N}^3$\\
                \hline
            \end{tabular}
        \end{table}
                """
        summary_table_2d = Tex(summary_table_2d_str, font_size=40)
        self.play(
            Write(summary_table_2d)
        )
        self.wait()

        self.play(
            summary_table_2d.animate.shift(LEFT * 3)
        )
        self.wait()

        # ADD THE CIRCLES
        circ_left = Circle(radius=1.5, color=ORANGE).shift(RIGHT * 2.5)
        circ_right = Circle(radius=1, color=GREEN).shift(RIGHT * 4.5)
        self.play(
            GrowFromCenter(circ_left),
            GrowFromCenter(circ_right)
        )
        self.wait()

        # add the lines
        circ_left_radius = DashedLine(
            start=circ_left.get_center(),
            end=circ_left.get_top()
        )
        circ_left_radius_label = BraceLabel(
            circ_left_radius, "R",
            brace_direction=LEFT
        )
        circ_right_radius = DashedLine(
            start=circ_right.get_center(),
            end=circ_right.get_top()
        )
        circ_right_radius_label = BraceLabel(
            circ_right_radius, "r",
            brace_direction=RIGHT
        )
        self.play(
            Create(circ_left_radius),
            Create(circ_right_radius),
            run_time=0.5,
            rate_func=smooth
        )
        self.play(
            Write(circ_left_radius_label),
            Write(circ_right_radius_label)
        )
        self.wait()

        # add the distance in between
        circ_distance = DashedLine(
            start=circ_left.get_center(),
            end=circ_right.get_center()
        )
        circ_distance_label = BraceLabel(
            circ_distance, "d"
        )
        self.play(
            Create(circ_distance),
            run_time=0.5,
            rate_func=smooth
        )
        self.play(
            Write(circ_distance_label)
        )
        self.wait()

        circ_overlap_formula = Tex(
            "$r^2 \\cos ^{-1}\\left(\\frac{d^2+r^2-R^2}{2 d r}\\right)+$",
            "$R^2 \\cos ^{-1}\\left(\\frac{d^2-r^2+R^2}{2 d R}\\right)-$",
            "$\\frac{1}{2} \\sqrt{(d+r-R) (d-r+R) (-d+r+R) (d+r+R)}$",
            font_size=30
        ).shift(DOWN * 2.5)
        self.play(
            Write(circ_overlap_formula)
        )
        self.wait()


class TestMakeSenseOfTheseEquations(Scene):
    def construct(self):
        question = Tex(
            "How can we make sense of those formulas?", font_size=50
        )

        self.play(Write(question), run_time=0.5)
        self.wait()

        self.play(
            question.animate.scale(0.6).to_corner(UL)
        )
        self.wait()

        num_line = NumberLine(x_range=[0, 1], length=10)
        self.play(
            Create(num_line)
        )
        self.wait()

        sep_tick_1_tracker = ValueTracker(0.2)
        sep_tick_2_tracker = ValueTracker(0.6)
        p_nothing = always_redraw(
            lambda:
            VGroup(
                Line(
                    start=np.array([0, 0.1, 0]),
                    end=np.array([0, -0.1, 0]),
                    stroke_width=2
                ).move_to(num_line.n2p(sep_tick_1_tracker.get_value())),
                Line(
                    start=np.array([0, 0.1, 0]),
                    end=np.array([0, -0.1, 0]),
                    stroke_width=2
                ).move_to(num_line.n2p(0)).set_opacity(0)
            )
        )

        p_alpha = always_redraw(
            lambda:
            VGroup(
                Line(
                    start=np.array([0, 0.1, 0]),
                    end=np.array([0, -0.1, 0]),
                    stroke_width=2
                ).move_to(num_line.n2p(sep_tick_1_tracker.get_value())),
                Line(
                    start=np.array([0, 0.1, 0]),
                    end=np.array([0, -0.1, 0]),
                    stroke_width=2
                ).move_to(num_line.n2p(sep_tick_2_tracker.get_value()))
            )
        )
        p_beta = always_redraw(
            lambda:
            VGroup(
                Line(
                    start=np.array([0, 0.1, 0]),
                    end=np.array([0, -0.1, 0]),
                    stroke_width=2
                ).move_to(num_line.n2p(sep_tick_2_tracker.get_value())),
                Line(
                    start=np.array([0, 0.1, 0]),
                    end=np.array([0, -0.1, 0]),
                    stroke_width=2
                ).move_to(num_line.n2p(1))
            )
        )

        p_nothing_label = always_redraw(
            lambda:
            BraceLabel(p_nothing, "P(\\text{nothing})", brace_direction=UP)
        )
        p_alpha_label = always_redraw(
            lambda:
            BraceLabel(p_alpha, "P(\\alpha)", brace_direction=UP)
        )
        p_beta_label = always_redraw(
            lambda:
            BraceLabel(p_beta, "P(\\beta)", brace_direction=UP)
        )
        p_tot_label = BraceLabel(num_line, "1")

        self.play(
            Create(p_nothing), run_time=0.5
        )
        self.play(
            Write(p_nothing_label)
        )
        self.play(
            Create(p_alpha), run_time=0.5
        )
        self.play(
            Write(p_alpha_label)
        )
        self.play(
            Create(p_beta), run_time=0.5
        )
        self.play(
            Write(p_beta_label)
        )
        self.wait()
        self.play(
            Write(p_tot_label), run_time=0.5
        )
        self.wait()
        self.play(
            Unwrite(p_tot_label), run_time=0.5
        )
        self.wait()

        self.play(
            sep_tick_2_tracker.animate.set_value(0.8)
        )
        self.wait()
        self.play(
            sep_tick_2_tracker.animate.set_value(0.4)
        )
        self.wait()

        self.play(
            num_line.animate.shift(UP)
        )
        self.wait()

        # write out the proposition
        proposition_3 = Tex(
            r"$\frac{P(\alpha)}{P(\beta)}=\frac{\text{\# $\alpha$ nuc. events in $||\Omega_\alpha||$}}{\text{\# $\beta$ nuc. events in $||\Omega_\beta||$}}$",
            r" $=\frac{J_\alpha\cdot ||\Omega_\alpha||}{J_\beta\cdot ||\Omega_\beta||}$",
            r" $\propto\frac{J_\alpha\cdot{\dot{R_\alpha}}^n}{J_\beta\cdot{\dot{R_\beta}}^n}$"
        )
        proposition_2 = Tex(
            r"$\frac{P(\alpha)}{P(\beta)}=\frac{\text{\# $\alpha$ nuc. events in $||\Omega_\alpha||$}}{\text{\# $\beta$ nuc. events in $||\Omega_\beta||$}}$",
            r" $=\frac{J_\alpha\cdot ||\Omega_\alpha||}{J_\beta\cdot ||\Omega_\beta||}$"
        ).align_to(proposition_3, direction=LEFT)
        proposition_1 = Tex(
            r"$\frac{P(\alpha)}{P(\beta)}=\frac{\text{\# $\alpha$ nuc. events in $||\Omega_\alpha||$}}{\text{\# $\beta$ nuc. events in $||\Omega_\beta||$}}$"
        ).align_to(proposition_3, direction=LEFT)
        self.play(
            Write(proposition_1)
        )
        self.wait()
        self.play(
            TransformMatchingTex(proposition_1, proposition_2)
        )
        self.wait()
        self.play(
            TransformMatchingTex(proposition_2, proposition_3)
        )
        self.wait()

        # add the time cone visualization

        axes = Axes(
            x_range=[-3, 3],
            y_range=[0, 4],
            x_length=6,
            y_length=4
        ).shift(DOWN * 2).scale(0.7)

        self.play(
            Create(axes)
        )
        self.wait()

        dot = Dot(
            point=axes.c2p(0, 3)
        )
        time_cone_left = DashedLine(
            start=axes.c2p(0, 3),
            end=axes.c2p(-2, 0)
        )
        time_cone_right = DashedLine(
            start=axes.c2p(0, 3),
            end=axes.c2p(2, 0)
        )

        division_point_tracker = ValueTracker(0)

        self.play(
            Create(dot)
        )
        self.play(
            Create(time_cone_left),
            Create(time_cone_right)
        )
        self.wait()

        # get the division of the time cone
        alpha_wins_portion = always_redraw(
            lambda:
            Polygon(
                axes.c2p(0, 3),
                axes.c2p(2, 0),
                axes.c2p(division_point_tracker.get_value(), 0),
                color=ORANGE, fill_opacity=0.7
            )
        )
        beta_wins_portion = always_redraw(
            lambda:
            Polygon(
                axes.c2p(0, 3),
                axes.c2p(-2, 0),
                axes.c2p(division_point_tracker.get_value(), 0),
                color=GREEN, fill_opacity=0.7
            )
        )
        self.play(
            DrawBorderThenFill(alpha_wins_portion)
        )
        self.wait()
        self.play(
            DrawBorderThenFill(beta_wins_portion)
        )
        self.wait()

        alpha_points = [
            Dot(point=axes.c2p(1, 1), color=ORANGE),
            Dot(point=axes.c2p(1.7, 0.2), color=ORANGE),
            Dot(point=axes.c2p(0.5, 0.6), color=ORANGE),
            Dot(point=axes.c2p(0.2, 2), color=ORANGE)
        ]
        beta_points = [
            Dot(point=axes.c2p(-1, 1), color=GREEN),
            Dot(point=axes.c2p(-0.6, 0.5), color=GREEN),
        ]

        self.play(
            LaggedStart(
                *[
                    FadeIn(point) for point in alpha_points[:2]
                ],
                lag_ratio=0.4,
                run_time=0.5
            ),
            LaggedStart(
                *[
                    FadeIn(point) for point in beta_points[:2]
                ],
                lag_ratio=0.4,
                run_time=0.5
            )
        )
        self.wait()

        self.play(
            LaggedStart(
                *[
                    FadeIn(point) for point in alpha_points[2:]
                ],
                lag_ratio=0.4,
                run_time=0.25
            )
        )
        self.wait()

        self.play(
            division_point_tracker.animate.set_value(0.5)
        )
        self.play(
            division_point_tracker.animate.set_value(-0.7)
        )
        self.wait()
