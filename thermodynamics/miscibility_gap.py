import numpy as np
from manim import *
from scipy.optimize import fsolve


def common_tangent(x: list, omega, temperature, gas_constant=8.31):
    return [((
                     omega * x[1] * (1 - x[1]) + gas_constant * temperature * (x[1] * np.log(x[1]) +
                                                                               (1 - x[1]) * np.log(1 - x[1]))
             ) - (
                     omega * x[0] * (1 - x[0]) + gas_constant * temperature * (x[0] * np.log(x[0]) +
                                                                               (1 - x[0]) * np.log(1 - x[0]))
             )) / (x[1] - x[0]) - (
                    omega * (1 - 2 * x[1]) + gas_constant * temperature * (np.log(x[1]) - np.log(1 - x[1]))
            ),
            (
                    omega * (1 - 2 * x[1]) + gas_constant * temperature * (np.log(x[1]) - np.log(1 - x[1]))
            ) - (
                    omega * (1 - 2 * x[0]) + gas_constant * temperature * (np.log(x[0]) - np.log(1 - x[0]))
            )
            ]


def get_y_coord_on_gibbs(x: float, omega: float, temperature: float, gas_constant: float = 8.31) -> float:
    return omega * x * (1 - x) + gas_constant * temperature * (x * np.log(x) +
                                                               (1 - x) * np.log(1 - x))


class GibbsFreeEnergy(Scene):
    def construct(self):
        # define the plotting axes
        axes = Axes(
            x_range=[0, 1, 0.1],
            y_range=[-1200, 1400, 200],
            x_length=6,
            y_length=8,
            tips=False
        )

        temperature_y_label = axes.get_y_axis_label(
            Tex("$G$"), edge=LEFT, direction=LEFT, buff=0.4
        )
        composition_label = axes.get_x_axis_label(Tex("$X$"))

        # define some variables
        gas_constant = 8.31
        temperature = 190
        omega = 5200

        # define enthalpy of mixing
        h_mix = axes.plot(
            lambda x: omega * x * (1 - x),
            x_range=[1E-6, 1 - 1E-6]
        )

        # define -T * entropy of mixing
        combined_s_mixed = axes.plot(
            lambda x: gas_constant * temperature * (x * np.log(x) + (1 - x) * np.log(1 - x)),
            x_range=[1E-6, 1 - 1E-6]
        )

        num_of_lines = 40
        hlines_combined_s_mixed = axes.get_vertical_lines_to_graph(
            combined_s_mixed, x_range=[1E-6, 1 - 1E-6], num_lines=num_of_lines, color=BLUE
        )

        # define gibbs free energy
        gibbs_free_energy = axes.plot(
            lambda x: omega * x * (1 - x) + gas_constant * temperature * (x * np.log(x) + (1 - x) * np.log(1 - x)),
            x_range=[1E-6, 1 - 1E-6]
        )

        self.add(axes)
        self.play(
            LaggedStart(
                *[Write(label) for label in [temperature_y_label, composition_label]],
                lag_ratio=1
            )
        )
        self.wait()

        self.play(
            Create(h_mix)
        )
        self.play(
            Create(combined_s_mixed)
        )
        self.wait()

        self.play(
            Create(hlines_combined_s_mixed, run_time=2)
        )
        self.wait()



        shift_up_vectors = [
            axes.c2p(
                x_coord, omega * x_coord * (1 - x_coord)
            ) - axes.c2p(
                x_coord, 0
            ) for x_coord in np.linspace(1E-6, 1 - 1E-6, num_of_lines)
        ]

        self.play(
            LaggedStart(
                *[hline.animate.shift(shift_up_vector)
                  for hline, shift_up_vector in zip(hlines_combined_s_mixed, shift_up_vectors)],
                lag_ratio=0.3
            ),
            run_time=3
        )
        self.wait()

        self.play(
            TransformFromCopy(
                h_mix, gibbs_free_energy
            )
        )
        self.wait()


class MiscibilityGap(Scene):
    def construct(self):
        # define the plotting axis for Gibbs free energy vs. composition
        axes_gibbs = Axes(
            x_range=[0, 1, 0.1],
            y_range=[-600, 300, 100],
            x_length=6,
            y_length=3,
            tips=False
        ).shift(UP * 2)

        gibbs_y_label = axes_gibbs.get_y_axis_label(
            Tex("$G$"), edge=LEFT, direction=LEFT, buff=0.4
        )

        # define the plotting axis for temperature vs. composition
        axes_temperature = Axes(
            x_range=[0, 1, 0.1],
            y_range=[180, 320, 20],
            x_length=6,
            y_length=3.5,
            tips=False
        ).shift(DOWN * 2)

        temperature_y_label = axes_temperature.get_y_axis_label(
            Tex("$T$"), edge=LEFT, direction=LEFT, buff=0.4
        )
        composition_label = axes_temperature.get_x_axis_label(Tex("$X$"))

        self.add(axes_gibbs, axes_temperature)
        self.play(
            LaggedStart(
                *[Write(label) for label in [temperature_y_label, composition_label, gibbs_y_label]],
                lag_ratio=1
            )
        )
        self.wait()

        # define the Gibbs free energy as a function of composition at a given temperature and omega
        gas_constant = 8.31
        temperature_tracker = ValueTracker(190)
        omega_tracker = ValueTracker(5200)

        gibbs_curve = always_redraw(
            lambda:
            axes_gibbs.plot(
                lambda x: omega_tracker.get_value() * x * (1 - x) + gas_constant * temperature_tracker.get_value() * (
                        x * np.log(x) + (1 - x) * np.log(1 - x)),
                x_range=[1E-6, 1 - 1E-6]
            )
        )

        common_tangent_1 = always_redraw(
            lambda:
            Dot(point=axes_gibbs.c2p(
                fsolve(common_tangent,
                       np.array([0.1, 0.9]),
                       args=(omega_tracker.get_value(), temperature_tracker.get_value())
                       )[0],
                get_y_coord_on_gibbs(
                    fsolve(common_tangent,
                           np.array([0.1, 0.9]),
                           args=(omega_tracker.get_value(), temperature_tracker.get_value())
                           )[0],
                    omega_tracker.get_value(),
                    temperature_tracker.get_value()
                )
            ), color=RED)
        )

        common_tangent_2 = always_redraw(
            lambda:
            Dot(point=axes_gibbs.c2p(
                fsolve(common_tangent,
                       np.array([0.1, 0.9]),
                       args=(omega_tracker.get_value(), temperature_tracker.get_value())
                       )[1],
                get_y_coord_on_gibbs(
                    fsolve(common_tangent,
                           np.array([0.1, 0.9]),
                           args=(omega_tracker.get_value(), temperature_tracker.get_value())
                           )[1],
                    omega_tracker.get_value(),
                    temperature_tracker.get_value()
                )
            ), color=RED)
        )

        common_tangent_trace_1 = always_redraw(
            lambda:
            Dot(
                point=axes_temperature.c2p(
                    fsolve(common_tangent,
                           np.array([0.1, 0.9]),
                           args=(omega_tracker.get_value(), temperature_tracker.get_value())
                           )[0],
                    temperature_tracker.get_value()
                ),
                color=BLUE
            )
        )

        common_tangent_trace_2 = always_redraw(
            lambda:
            Dot(
                point=axes_temperature.c2p(
                    fsolve(common_tangent,
                           np.array([0.1, 0.9]),
                           args=(omega_tracker.get_value(), temperature_tracker.get_value())
                           )[1],
                    temperature_tracker.get_value()
                ),
                color=BLUE
            )
        )

        trace_1_path = TracedPath(common_tangent_trace_1.get_center, stroke_color=BLUE, stroke_width=4)
        trace_2_path = TracedPath(common_tangent_trace_2.get_center, stroke_color=BLUE, stroke_width=4)

        hline_trace_1_to_2 = always_redraw(
            lambda:
            DashedLine(
                start=common_tangent_1.get_center(),
                end=common_tangent_2.get_center(),
                color=RED
            )
        )

        vline_trace_1 = always_redraw(
            lambda:
            DashedLine(
                start=common_tangent_1.get_center(),
                end=common_tangent_trace_1.get_center()
            ).set_color_by_gradient(RED, BLUE)
        )

        vline_trace_2 = always_redraw(
            lambda:
            DashedLine(
                start=common_tangent_2.get_center(),
                end=common_tangent_trace_2.get_center()
            ).set_color_by_gradient(RED, BLUE)
        )

        self.play(
            Create(gibbs_curve)
        )
        self.wait()

        self.play(
            Create(hline_trace_1_to_2)
        )

        self.play(
            FadeIn(common_tangent_1, common_tangent_2)
        )
        self.wait()

        self.play(
            Create(vline_trace_1), Create(vline_trace_2)
        )
        self.play(
            FadeIn(common_tangent_trace_1, common_tangent_trace_2)
        )
        self.add(trace_1_path, trace_2_path)
        self.wait()

        self.play(
            temperature_tracker.animate.set_value(312.9), run_time=3, rate_func=linear
        )
        self.wait()
