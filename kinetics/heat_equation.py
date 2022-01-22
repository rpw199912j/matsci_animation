import numpy as np
from manim import *


class FourierSeries(Scene):
    def construct(self):
        # define the plotting axes
        axes = Axes(
            x_range=[0, 1, 0.5],
            y_range=[-1.5, 1.5, 0.5],
            x_length=8,
            y_length=6,
            tips=False
        )

        self.add(axes)

        # define the variable
        num_sum_terms = 50
        length = 0.5

        # define the Fourier series solution
        heat_eq_sol = axes.plot(
            lambda x, j=num_sum_terms, l=length:
            4 / PI * np.sum([1 / (2 * _ + 1) * np.sin((2 * _ + 1) * PI * x / l) for _ in range(j + 1)]),
            x_range=(0, length * 2, 1E-3)
        )

        lines = axes.get_vertical_lines_to_graph(
            heat_eq_sol, x_range=[0, length * 2], num_lines=100, color=BLUE
        )

        self.play(
            Create(heat_eq_sol)
        )
        self.wait()

        self.play(
            Create(lines)
        )
        self.wait()
