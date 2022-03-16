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

    @staticmethod
    def get_z_value(x, y):
        return x + y

    def construct(self):
        # define the axes
        x_max, z_max = 4, 8
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
        sampling_number = 10
        sampling_unit_increment = x_max / sampling_number
        dots_on_2d_plane = [Dot(point=axes_3d.c2p((row_index + 0.5) * sampling_unit_increment,
                                                  (col_index + 0.5) * sampling_unit_increment, 0),
                                radius=(axes_3d.c2p(0.1 * sampling_unit_increment, 0, 0) - axes_3d.c2p(0, 0, 0))[0])
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
                self.get_z_value((row_index + 0.5) * sampling_unit_increment,
                                 (col_index + 0.5) * sampling_unit_increment)
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
                u, v, self.get_z_value(u, v)
            ),
            resolution=(sampling_number, sampling_number),
            u_range=[0, x_max],
            v_range=[0, x_max],
            stroke_opacity=0
        )
        surface.set_fill_by_value(
            axes=axes_3d, colors=[PURE_GREEN, "#f06b0c"], axis=2
        )
        surface.set_style(fill_opacity=1)

        self.play(
            Create(surface),
            FadeOut(dots_on_2d_plane)
        )
        self.wait()

        # delete
        test_plot = axes_3d.plot(
            lambda u: 1,
            x_range=[0, x_max]
        )
        self.play(
            Create(test_plot)
        )
        self.wait()

        # rotate the camera
        # self.renderer.camera._frame_center.move_to(z_axis)
        self.begin_ambient_camera_rotation(
            rate=0.2
        )
        self.wait(5)
        self.stop_ambient_camera_rotation()
        self.wait()
