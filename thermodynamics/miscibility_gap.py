import numpy as np
from manim import *


class GibbsFreeEnergy(Scene):
    def construct(self):
        # define the plotting axis
        axes = Axes(
            x_range=[0, 1, 0.1],
            y_range=[-1000, 500, 100],
            x_length=6,
            y_length=8,
            tips=False
        )

        self.add(axes)

        # define the Gibbs Free energy as a function of temperature
        temperature = ValueTracker(360)
        omega = ValueTracker(5200)

        gibbs_curve = always_redraw(
            lambda:
            axes.plot(
                lambda x: omega.get_value() * x * (1 - x) + 8.31 * temperature.get_value() * (
                            x * np.log(x) + (1 - x) * np.log(1 - x)),
                x_range=[1E-6, 1-1E-6]
            )
        )

        self.play(
            Create(gibbs_curve)
        )
        self.wait()

        self.play(
            temperature.animate.set_value(0), run_time=4, rate_func=linear
        )
        self.wait()
