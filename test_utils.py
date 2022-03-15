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


class TestParametricSurface(ThreeDScene):
    def construct(self):
        # define the axes
        axes_3d = ThreeDAxes(
            x_range=[0, 3, 1],
            y_range=[0, 3, 1],
            z_range=[0, 6, 1],
            x_length=4,
            y_length=4,
            z_length=4,
            tips=False
        )
        x_axis = axes_3d.get_x_axis()
        y_axis = axes_3d.get_y_axis()
        z_axis = axes_3d.get_z_axis()

        x_label = axes_3d.get_x_axis_label(Tex("$\\log_{10}\\left(\\frac{\\dot{R_\\alpha}}{\\dot{R_\\beta}}\\right)$"))
        y_label = axes_3d.get_y_axis_label(Tex("$\\log_{10}\\left(\\frac{J_\\alpha}{J_\\beta}\\right)$"))
        z_label = axes_3d.get_z_axis_label(Tex(
            "$\\log_{10}\\left(\\frac{P_\\alpha}{P_\\beta}\\right)$",
            font_size=30
        ))

        # add the quasi-2d axes
        self.add(x_axis, y_axis, x_label, y_label)
        self.wait()

        # create the dots
        sampling_unit_increment = 3 / 10
        dots_on_2d_plane = [Dot(point=axes_3d.c2p((row_index + 0.5) * sampling_unit_increment,
                                                  (col_index + 0.5) * sampling_unit_increment, 0),
                                radius=(axes_3d.c2p(0.1 * sampling_unit_increment, 0, 0) - axes_3d.c2p(0, 0, 0))[0])
                            for row_index in range(10) for col_index in range(10)
                            ]
        dots_on_2d_plane = VGroup(*dots_on_2d_plane)
        self.play(
            FadeIn(dots_on_2d_plane)
        )

        # re-orient the camera to move from quasi-2d to 3d
        axes_3d_shift_vector = -5 * Y_AXIS - Z_AXIS
        self.move_camera(
            phi=80 * DEGREES, theta=-90 * DEGREES,
            added_anims=[
                VGroup(
                    x_axis, y_axis, x_label, y_label, dots_on_2d_plane
                ).animate.shift(axes_3d_shift_vector)
            ]
        )

        # add the third axis and its label
        self.play(
            LaggedStart(
                *[
                    Create(z_axis.shift(axes_3d_shift_vector)),
                    Write(z_label.shift(axes_3d_shift_vector))
                ],
                lag_ratio=1.2
            )
        )

        # shift the position of each dot
        shift_vectors = [
            axes_3d.c2p(
                (row_index + 0.5) * sampling_unit_increment,
                (col_index + 0.5) * sampling_unit_increment,
                (row_index + 0.5) * sampling_unit_increment + (col_index + 0.5) * sampling_unit_increment
            ) - axes_3d.c2p(
                (row_index + 0.5) * sampling_unit_increment,
                (col_index + 0.5) * sampling_unit_increment,
                0
            )
            for row_index in range(10) for col_index in range(10)
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
                u, v, u + v
            ),
            resolution=(10, 10),
            u_range=[0, 3],
            v_range=[0, 3]
        )
        surface.set_style(fill_opacity=1, stroke_opacity=0)
        surface.set_fill_by_value(
            axes=axes_3d, colors=[(GREEN_E, 2), (RED, 5)], axis=2
        )
        self.play(
            Create(surface)
        )
        self.wait()

        self.begin_ambient_camera_rotation(

        )
        self.wait(5)
        self.stop_ambient_camera_rotation()
        self.wait()

