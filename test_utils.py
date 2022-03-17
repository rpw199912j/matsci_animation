import pandas as pd
from itertools import product
from manim import *


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
