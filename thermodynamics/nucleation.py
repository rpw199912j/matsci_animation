import numpy as np
from manim import *


class CriticalRadius(Scene):
    def construct(self):
        # define the plotting axes
        x_max = 1.2
        y_max = 1
        axes = Axes(
            x_range=[0, x_max, 0.2],
            y_range=[-y_max, y_max, 0.5],
            x_length=8,
            y_length=6,
            tips=False
        )

        # get the unit y length
        unit_y_length = (axes.c2p(0, 1) - axes.c2p(0, 0))[1]

        x_label = axes.get_x_axis_label(Tex("$r$"))
        y_label = axes.get_y_axis_label(Tex("$\Delta G$"))

        # define the variables
        num_of_lines = 100
        delta_G = -0.18
        interfacial_energy = 0.054
        r_crit = 2 * interfacial_energy / np.abs(delta_G)
        dot_radius = 0.3
        y_shift = dot_radius / (unit_y_length * 2.86)
        print(unit_y_length)

        # define the volume Gibbs free energy term
        volume_G = axes.plot(
            lambda r:
            4 / 3 * np.pi * r ** 3 * delta_G,
            x_range=(0, x_max), color=RED
        )

        # define the surface Gibbs free energy term
        surface_G = axes.plot(
            lambda r:
            4 * np.pi * r ** 2 * interfacial_energy,
            x_range=(0, x_max), color=GREEN
        )

        # define the total energy
        total_G = axes.plot(
            lambda r:
            4 / 3 * np.pi * r ** 3 * delta_G + 4 * np.pi * r ** 2 * interfacial_energy,
            x_range=(0, x_max), color=YELLOW
        )

        self.play(
            Create(axes)
        )
        self.play(
            Write(x_label),
            Write(y_label)
        )
        self.wait()

        self.play(
            Create(volume_G)
        )
        self.wait()

        self.play(
            Create(surface_G)
        )
        self.wait()

        # create the vertical line
        volume_G_lines = axes.get_vertical_lines_to_graph(
            volume_G, x_range=[0, x_max], color=BLUE,
            num_lines=num_of_lines, line_config={"dashed_ratio": 0.9}
        )
        self.play(
            Create(volume_G_lines), run_time=2
        )
        self.wait()

        shift_vectors = [
            axes.c2p(
                x_coord, 4 * np.pi * x_coord ** 2 * interfacial_energy
            ) - axes.c2p(
                x_coord, 0
            )
            for x_coord in np.linspace(0, x_max, num_of_lines)
        ]

        # shift all the lines in sin_curve_lines
        self.play(
            LaggedStart(
                *[vline.animate.shift(shift_vector)
                  for vline, shift_vector in zip(volume_G_lines, shift_vectors)],
                lag_ratio=0.3
            ),
            run_time=3
        )
        self.wait()

        # morph the surface_G line to the total_G line
        self.play(
            TransformFromCopy(surface_G, total_G)
        )
        self.wait()

        # remove the surface_G, volume_G and vertical lines
        self.play(
            FadeOut(surface_G, volume_G, *volume_G_lines)
        )
        self.wait()

        # add a circle
        dot_x_tracker = ValueTracker(r_crit)
        dot = always_redraw(
            lambda:
            Dot(
                point=axes.c2p(
                    dot_x_tracker.get_value(),
                    4 / 3 * np.pi * dot_x_tracker.get_value() ** 3 * delta_G + (
                            4 * np.pi * dot_x_tracker.get_value() ** 2 * interfacial_energy +
                            dot_x_tracker.get_value() / (unit_y_length * 2.84)
                    ),
                ),
                radius=dot_x_tracker.get_value() / unit_y_length
            )
        )
        self.play(
            DrawBorderThenFill(dot)
        )
        self.wait()

        self.play(
            dot_x_tracker.animate.set_value(0),
            run_time=2
        )
        self.wait()

        self.play(
            dot_x_tracker.animate.set_value(r_crit),
            run_time=2
        )
        self.wait()

        self.play(
            dot_x_tracker.animate.set_value(x_max),
            run_time=2
        )
        self.wait()
