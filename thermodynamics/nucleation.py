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
        delta_G = -0.44
        interfacial_energy = 0.16
        r_crit = 2 * interfacial_energy / np.abs(delta_G)
        delta_G_crit = 4 / 3 * np.pi * r_crit ** 3 * delta_G + 4 * np.pi * r_crit ** 2 * interfacial_energy

        # define the volume Gibbs free energy term
        volume_G = axes.plot(
            lambda r:
            4 / 3 * np.pi * r ** 3 * delta_G,
            x_range=(0, x_max), color=RED
        )
        volume_G_label = MathTex(
            r"\Delta G_{\text{vol}}",
            color=RED
        ).move_to(
            axes.c2p(0.4, -0.5)
        )

        # define the surface Gibbs free energy term
        surface_G = axes.plot(
            lambda r:
            4 * np.pi * r ** 2 * interfacial_energy,
            x_range=(0, x_max), color=GREEN
        )
        surface_G_label = MathTex(
            r"\Delta G_{\text{sur}}",
            color=GREEN
        ).move_to(
            axes.c2p(0.35, 0.5)
        )

        # define the total energy
        total_G = axes.plot(
            lambda r:
            4 / 3 * np.pi * r ** 3 * delta_G + 4 * np.pi * r ** 2 * interfacial_energy,
            x_range=(0, x_max), color=YELLOW
        )
        total_G_label = MathTex(
            r"\Delta G_{\text{tot}}",
            color=YELLOW
        ).move_to(
            axes.c2p(1, -0.25)
        )

        total_G_crit = axes.plot(
            lambda r:
            4 / 3 * np.pi * r ** 3 * delta_G + 4 * np.pi * r ** 2 * interfacial_energy,
            x_range=(0, r_crit), color=YELLOW
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
        self.play(
            Write(volume_G_label)
        )
        self.wait()

        self.play(
            Create(surface_G)
        )
        self.play(
            Write(surface_G_label)
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
        self.play(
            Write(total_G_label)
        )
        self.wait()

        # remove the surface_G, volume_G and vertical lines
        self.play(
            FadeOut(surface_G, volume_G, *volume_G_lines, surface_G_label, volume_G_label)
        )
        self.wait()

        # add a tangent line
        tangent_frac_tracker = ValueTracker(0)
        tangent = always_redraw(
            lambda:
            TangentLine(
                total_G, alpha=tangent_frac_tracker.get_value(),
                color=PURPLE_C, length=2, stroke_width=7
            )
        )

        self.play(
            Create(tangent)
        )
        self.wait()

        self.play(
            tangent_frac_tracker.animate.set_value(1),
            run_time=2
        )
        self.wait()

        self.play(
            tangent_frac_tracker.animate.set_value(
                total_G_crit.get_arc_length() / total_G.get_arc_length()
            )
        )
        self.wait()

        self.play(
            FadeOut(tangent), run_time=0.5
        )
        self.wait()

        # highlight the r_crit and G_crit
        lines_crit_1 = DashedLine(
            start=axes.c2p(r_crit, delta_G_crit),
            end=axes.c2p(0, delta_G_crit)
        )

        lines_crit_2 = DashedLine(
            start=axes.c2p(r_crit, delta_G_crit),
            end=axes.c2p(r_crit, 0)
        )

        self.play(
            Create(lines_crit_1),
            Create(lines_crit_2)
        )
        self.wait()

        # add the critical radius and Delta G label
        delta_G_crit_label = Tex(
            r"$\Delta G^*$"
        ).next_to(lines_crit_1.get_left(), direction=LEFT, buff=0.1)

        r_crit_label = Tex(
            r"$r^*$"
        ).next_to(lines_crit_2.get_bottom(), direction=DOWN, buff=0.1)

        self.play(
            Write(r_crit_label),
            Write(delta_G_crit_label)
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


# noinspection DuplicatedCode
class ChangeVariables(Scene):
    @staticmethod
    def get_r_crit(interfacial_energy, delta_G):
        r_crit = 2 * interfacial_energy / np.abs(delta_G)
        return r_crit

    @staticmethod
    def get_delta_G_crit(interfacial_energy, delta_G):
        r_crit = 2 * interfacial_energy / np.abs(delta_G)
        delta_G_crit = 4 / 3 * np.pi * r_crit ** 3 * delta_G + (
                4 * np.pi * r_crit ** 2 * interfacial_energy
        )
        return delta_G_crit

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
        ).shift(LEFT * 2)

        # get the axes labels
        x_label = axes.get_x_axis_label(Tex("$r$"))
        y_label = axes.get_y_axis_label(Tex("$\Delta G$"))

        # define the variables
        delta_G_tracker = ValueTracker(-0.44)
        interfacial_energy_tracker = ValueTracker(0.16)

        # define the volume Gibbs free energy term
        volume_G = always_redraw(
            lambda:
            axes.plot(
                lambda r:
                4 / 3 * np.pi * r ** 3 * delta_G_tracker.get_value(),
                x_range=(0, x_max), color=RED
            )
        )

        # define the surface Gibbs free energy term
        surface_G = always_redraw(
            lambda:
            axes.plot(
                lambda r:
                4 * np.pi * r ** 2 * interfacial_energy_tracker.get_value(),
                x_range=(0, x_max), color=GREEN
            )
        )

        # define the total energy
        total_G = always_redraw(
            lambda:
            axes.plot(
                lambda r:
                4 / 3 * np.pi * r ** 3 * delta_G_tracker.get_value() + (
                        4 * np.pi * r ** 2 * interfacial_energy_tracker.get_value()
                ),
                x_range=(0, x_max), color=YELLOW
            )
        )

        # highlight the r_crit and G_crit
        lines_crit_1 = always_redraw(
            lambda:
            DashedLine(
                start=axes.c2p(self.get_r_crit(interfacial_energy_tracker.get_value(),
                                               delta_G_tracker.get_value()),
                               self.get_delta_G_crit(
                                   interfacial_energy_tracker.get_value(),
                                   delta_G_tracker.get_value()
                               )),
                end=axes.c2p(0, self.get_delta_G_crit(
                    interfacial_energy_tracker.get_value(),
                    delta_G_tracker.get_value()
                ))
            )
        )

        lines_crit_2 = always_redraw(
            lambda:
            DashedLine(
                start=axes.c2p(self.get_r_crit(interfacial_energy_tracker.get_value(),
                                               delta_G_tracker.get_value()),
                               self.get_delta_G_crit(
                                   interfacial_energy_tracker.get_value(),
                                   delta_G_tracker.get_value()
                               )),
                end=axes.c2p(self.get_r_crit(interfacial_energy_tracker.get_value(),
                                             delta_G_tracker.get_value()),
                             0),
            )
        )

        # add the critical radius and Delta G label
        delta_G_crit_label = always_redraw(
            lambda:
            Tex(
                r"$\Delta G^*$"
            ).next_to(lines_crit_1.get_left(), direction=LEFT, buff=0.1)
        )

        r_crit_label = always_redraw(
            lambda:
            Tex(
                r"$r^*$"
            ).next_to(lines_crit_2.get_bottom(), direction=DOWN, buff=0.1)
        )

        self.add(
            axes, x_label, y_label, volume_G, surface_G, total_G,
            lines_crit_1, lines_crit_2, r_crit_label, delta_G_crit_label
        )
        self.wait()

        # write out the total delta G expression
        total_G_eq = MathTex(
            r"\Delta G", "&=",
            r"\frac{4}{3}\pi r^3\left(G_\beta-G_\alpha\right)\\",
            "&+",
            r"4\pi r^2\gamma"
        ).to_corner(UR)

        total_G_eq[0].set_color(YELLOW)
        total_G_eq[2].set_color(RED)
        total_G_eq[4].set_color(GREEN)

        bounding_box_1 = SurroundingRectangle(
            total_G_eq[2][6:], color=BLUE
        )
        arrow_1 = Arrow(
            start=bounding_box_1.get_corner(UR) + 0.2 * UP,
            end=bounding_box_1.get_corner(DR) + 0.2 * DOWN
        ).next_to(bounding_box_1, direction=RIGHT, buff=0.5 * DEFAULT_MOBJECT_TO_MOBJECT_BUFFER)

        bounding_box_2 = SurroundingRectangle(
            total_G_eq[4][4:], color=BLUE
        )
        arrow_2 = Arrow(
            start=bounding_box_2.get_corner(DR) + 0.2 * DOWN,
            end=bounding_box_2.get_corner(UR) + 0.2 * UP
        ).next_to(bounding_box_2, direction=RIGHT, buff=0.5 * DEFAULT_MOBJECT_TO_MOBJECT_BUFFER)

        self.play(Write(total_G_eq))
        self.wait()

        self.play(Create(bounding_box_1))
        self.wait()

        # increase the magnitude of delta_G
        self.play(
            FadeIn(arrow_1, run_time=0.5),
            delta_G_tracker.animate(run_time=2).set_value(-0.7),
        )
        self.play(
            FadeOut(arrow_1, run_time=0.5)
        )
        self.wait()

        # increase the magnitude of interfacial_energy
        self.play(
            Transform(bounding_box_1, bounding_box_2)
        )
        self.wait()

        self.play(
            FadeIn(arrow_2, run_time=0.5),
            interfacial_energy_tracker.animate(run_time=2).set_value(0.25)
        )
        self.play(
            FadeOut(arrow_2, run_time=0.5)
        )
        self.wait()
