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

        # define the variables
        num_of_lines = 200
        num_sum_terms = 0
        length = 0.5

        # define the Fourier series solution
        heat_eq_sol = axes.plot(
            lambda x, j=num_sum_terms, l=length:
            4 / PI * np.sum([1 / (2 * _ + 1) * np.sin((2 * _ + 1) * PI * x / l) for _ in range(j + 1)]),
            x_range=(0, length * 2, 1E-3)
        )

        self.play(
            Create(heat_eq_sol)
        )
        self.wait()

        for _ in range(4):
            # increment the number of terms to sum by 1
            num_sum_terms += 1

            # define the sine curve to be added
            sin_curve_to_add = axes.plot(
                lambda x, j=num_sum_terms, l=length:
                4 / PI / (2 * j + 1) * np.sin((2 * j + 1) * PI * x / l),
                x_range=(0, length * 2, 1E-3), color=ORANGE
            )

            sin_curve_lines = axes.get_vertical_lines_to_graph(
                sin_curve_to_add, x_range=[0, length * 2], num_lines=num_of_lines, color=BLUE,
                line_config={"dashed_ratio": 0.9}
            )

            self.play(
                Create(sin_curve_to_add)
            )
            self.wait()
            self.play(
                Create(sin_curve_lines)
            )
            self.wait()

            # determine for each sin_curve_line, the vector it should be shifted by
            shift_vectors = [
                axes.c2p(
                    x_coord, 4 / PI * np.sum([1 / (2 * _ + 1) * np.sin((2 * _ + 1) * PI * x_coord / length)
                                              for _ in range(num_sum_terms)])
                ) - axes.c2p(
                    x_coord, 0
                ) for x_coord in np.linspace(0, length * 2, num_of_lines)
            ]
            # shift all the lines in sin_curve_lines
            self.play(
                LaggedStart(
                    *[vline.animate.shift(shift_vector)
                      for vline, shift_vector in zip(sin_curve_lines, shift_vectors)],
                    lag_ratio=0.3
                ),
                run_time=2
            )
            self.wait()

            # define the new Fourier series solution with one more term added
            heat_eq_sol_new = axes.plot(
                lambda x, j=num_sum_terms, l=length:
                4 / PI * np.sum([1 / (2 * _ + 1) * np.sin((2 * _ + 1) * PI * x / l) for _ in range(j + 1)]),
                x_range=(0, length * 2, 1E-3)
            )

            # add a smooth transition from the previous solution to the current one
            self.play(
                TransformFromCopy(heat_eq_sol, heat_eq_sol_new)
            )
            self.wait()

            # remove the previous solution, and the new sine curve and the corresponding vertical lines
            self.play(
                FadeOut(heat_eq_sol, sin_curve_to_add, sin_curve_lines)
            )
            self.wait()

            # assign the new solution as the previous solution for the next iteration
            heat_eq_sol = heat_eq_sol_new

        num_sum_terms_tracker = ValueTracker(num_sum_terms)

        heat_eq_sol_updated = always_redraw(
            lambda:
            axes.plot(
                lambda x:
                4 / PI * np.sum([1 / (2 * _ + 1) * np.sin((2 * _ + 1) * PI * x / length)
                                 for _ in range(int(np.floor(num_sum_terms_tracker.get_value())) + 1)]),
                x_range=(0, length * 2, 1E-3)
            )
        )

        # AHA: need to be careful about using Transform on updater object
        self.add(heat_eq_sol_updated)
        self.remove(heat_eq_sol)

        self.play(
            num_sum_terms_tracker.animate.set_value(200), run_time=5, rate_func=linear
        )
        self.wait()
