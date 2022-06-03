import numpy as np
import pandas as pd
from manim import *
from scipy import interpolate


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


class MathDerivation(Scene):
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

        # define the derivative of total energy wrt r
        total_G_deriv = axes.plot_derivative_graph(
            total_G, color=BLUE
        )

        total_G_deriv_label = MathTex(
            r"\frac{\partial}{\partial r}\left(\Delta G_{\text{tot}}\right)",
            color=BLUE
        ).move_to(
            axes.c2p(0.6, -0.5)
        )

        self.add(
            axes, x_label, y_label, total_G, total_G_label
        )
        self.wait()

        self.play(
            Create(total_G_deriv)
        )
        self.play(
            Write(total_G_deriv_label)
        )
        self.wait()

        # highlight the point where derivative is equal to 0
        deriv_0 = Dot(
            point=axes.c2p(r_crit, 0), color=PURPLE
        )
        self.play(
            LaggedStart(
                DrawBorderThenFill(deriv_0),
                Flash(deriv_0),
                lag_ratio=0.5
            )

        )
        self.wait()

        # shift everything to the left
        self.play(
            *[mob.animate.shift(LEFT * 2.5) for mob in self.mobjects]
        )
        self.wait()

        # write out the equations
        total_G_eq = MathTex(
            r"\Delta G", "&=",
            r"\frac{4}{3}", r"\pi", "r^3", r"\left(G_\beta-G_\alpha\right)\\",
            "&+",
            "4", r"\pi", "r^2", r"\gamma"
        ).to_corner(UR)

        total_G_deriv_eq_1 = MathTex(
            r"\frac{\partial}{\partial r}", r"\left(", r"\Delta G", r"\right)", "&=",
            "4", r"\pi", "r^2", r"\left(G_\beta-G_\alpha\right)\\",
            "&+",
            "8", r"\pi", "r", r"\gamma"
        ).to_corner(UR)

        total_G_deriv_eq_2 = MathTex(
            "0", "&=",
            "4", r"\pi", "r^2", r"\left(G_\beta-G_\alpha\right)\\",
            "&+",
            "8", r"\pi", "r", r"\gamma"
        ).to_corner(UR).align_to(
            total_G_deriv_eq_1, UP
        )

        total_G_deriv_eq_3 = MathTex(
            "-", "8", r"\pi", "r", r"\gamma", "&=",
            "4", r"\pi", "r^2", r"\left(G_\beta-G_\alpha\right)\\",
        ).to_corner(UR).align_to(
            total_G_deriv_eq_2, UP
        )

        total_G_deriv_eq_4 = MathTex(
            "-", "8", r"\pi", "r", r"\gamma", "&=",
            "4", r"\pi", "r^2", r"\left(G_\beta-G_\alpha\right)\\",
            r"\frac{-8\pi r \gamma}{4\pi r^2}", "&=", r"\left(G_\beta-G_\alpha\right)\\"
        ).to_corner(UR).align_to(
            total_G_deriv_eq_3, UP
        )

        total_G_deriv_eq_5 = MathTex(
            "-", "8", r"\pi", "r", r"\gamma", "&=",
            "4", r"\pi", "r^2", r"\left(G_\beta-G_\alpha\right)\\",
            r"\frac{-2 \gamma}{r}", "&=", r"\left(G_\beta-G_\alpha\right)\\"
        ).to_corner(UR).align_to(
            total_G_deriv_eq_3, UP
        )

        total_G_deriv_eq_6 = MathTex(
            "-", "8", r"\pi", "r", r"\gamma", "&=",
            "4", r"\pi", "r^2", r"\left(G_\beta-G_\alpha\right)\\",
            r"\frac{-2 \gamma}{r}", "&=", r"\left(G_\beta-G_\alpha\right)\\",
            "r^*", "&=", r"\frac{-2 \gamma}{G_\beta-G_\alpha}"
        ).to_corner(UR).align_to(
            total_G_deriv_eq_3, UP
        )

        self.play(
            Write(total_G_eq)
        )
        self.wait()

        self.play(
            TransformMatchingTex(total_G_eq, total_G_deriv_eq_1)
        )
        self.wait()

        self.play(
            TransformMatchingTex(total_G_deriv_eq_1, total_G_deriv_eq_2)
        )
        self.wait()

        self.play(
            TransformMatchingTex(total_G_deriv_eq_2, total_G_deriv_eq_3)
        )
        self.wait()

        self.play(
            Write(total_G_deriv_eq_4[-3:])
        )
        self.wait()

        self.play(
            TransformMatchingShapes(total_G_deriv_eq_4[-3:], total_G_deriv_eq_5[-3:])
        )
        self.wait()

        self.play(
            Write(total_G_deriv_eq_6[-3:])
        )
        self.wait()


class CrossOver(ZoomedScene):
    def __init__(self):
        ZoomedScene.__init__(
            self,
            zoomed_display_width=5,
            zoomed_display_height=3
        )

    def construct(self):
        # define the axes
        x_min, x_max, x_step = 990, 1030, 5
        y_min, y_max, y_step = -46500, -43000, 500
        axes = Axes(
            x_range=[x_min, x_max, x_step],
            y_range=[y_min, y_max, y_step],
            x_length=8,
            y_length=6,
            tips=False
        )

        # add the axes labels
        x_label = axes.get_x_axis_label(Tex("$T$"))
        y_label = axes.get_y_axis_label(Tex("$G$"), edge=LEFT, direction=LEFT)

        self.add(axes)
        self.wait()

        self.play(
            Write(x_label),
            Write(y_label)
        )
        self.wait()

        # add the thermo-calc logo
        thermocalc_logo = SVGMobject(
            r"C:\Users\rpw19\PycharmProjects\matsci_animation\figure\ThermoCalc_logo.svg"
        ).set(height=config["frame_height"] * 0.15).set_color("#9b193b").to_corner(UL)
        self.play(
            DrawBorderThenFill(thermocalc_logo)
        )
        self.wait()

        # read in the gibbs energy data
        gibbs_df = pd.read_csv(
            r"C:\Users\rpw19\PycharmProjects\matsci_animation\data\gibbs_energy\gibbs_energy_per_mole.csv",
            names=["temp", "bcc", "fcc"]
        ).drop(index=0)
        gibbs_df = gibbs_df.astype(float)
        gibbs_df = gibbs_df[gibbs_df["temp"].between(x_min, x_max)]
        bcc_gibbs = gibbs_df[["temp", "bcc"]].to_numpy()
        fcc_gibbs = gibbs_df[["temp", "fcc"]].to_numpy()

        # fit a spline curve for bcc
        bcc_fit = interpolate.interp1d(
            bcc_gibbs[:, 0],
            bcc_gibbs[:, 1]
        )

        bcc_curve = axes.plot(
            lambda x: bcc_fit(x),
            x_range=[x_min, x_max],
            color=GREEN, stroke_width=2
        )

        # fit a spline curve for fcc
        fcc_fit = interpolate.interp1d(
            fcc_gibbs[:, 0],
            fcc_gibbs[:, 1]
        )

        fcc_curve = axes.plot(
            lambda x: fcc_fit(x),
            x_range=[x_min, x_max],
            color=RED, stroke_width=2
        )

        # find the intersection points
        intersect_x_lst = np.linspace(x_min, x_max, num=10000)
        bcc_y_lst = bcc_fit(intersect_x_lst)
        fcc_y_lst = fcc_fit(intersect_x_lst)
        intersect_x_idx = np.argmin(np.abs(bcc_y_lst - fcc_y_lst))
        intersect_x = intersect_x_lst[intersect_x_idx]

        # bcc_curve = axes.plot_line_graph(
        #     x_values=bcc_gibbs[:, 0],
        #     y_values=bcc_gibbs[:, 1],
        #     add_vertex_dots=False,
        #     line_color=GREEN,
        #     stroke_width=2
        # )
        # fcc_curve = axes.plot_line_graph(
        #     x_values=fcc_gibbs[:, 0],
        #     y_values=fcc_gibbs[:, 1],
        #     add_vertex_dots=False,
        #     line_color=RED,
        #     stroke_width=2
        # )

        self.play(
            Create(bcc_curve)
        )
        self.play(
            Create(fcc_curve)
        )

        # add the zoomed frame
        zoomed_camera = self.zoomed_camera
        zoomed_display = self.zoomed_display
        zoomed_camera_frame = zoomed_camera.frame
        zoomed_display_frame = zoomed_display.display_frame

        zoomed_camera_frame.set_color(YELLOW)
        zoomed_display_frame.set_color(YELLOW)

        # add the moving dots
        x_tracker = ValueTracker(x_min)
        bcc_dot = always_redraw(
            lambda:
            Dot(
                point=axes.i2gp(x_tracker.get_value(), bcc_curve),
                radius=0.03, color=GREEN
            )
        )
        fcc_dot = always_redraw(
            lambda:
            Dot(
                point=axes.i2gp(x_tracker.get_value(), fcc_curve),
                radius=0.03, color=RED
            )
        )

        zoomed_camera_frame.add_updater(
            lambda x: x.move_to(
                (bcc_dot.get_center() + fcc_dot.get_center()) / 2
            )
        )

        # add the G_beta - G_alpha labels
        delta_gibbs_label = MathTex(
            r"G_\beta", "-", r"G_\alpha", "="
        ).move_to(
            axes.c2p(1025, -45000)
        )

        delta_gibbs_label[0].set_color(GREEN)
        delta_gibbs_label[2].set_color(RED)

        delta_gibbs = DecimalNumber(
            bcc_fit(x_tracker.get_value()) - fcc_fit(x_tracker.get_value()),
            num_decimal_places=2,
            include_sign=True
        ).next_to(
            delta_gibbs_label, RIGHT
        )

        delta_gibbs_label_unit = Tex("J").next_to(delta_gibbs, RIGHT)

        delta_gibbs.add_updater(
            lambda d: d.set_value(
                bcc_fit(x_tracker.get_value()) - fcc_fit(x_tracker.get_value())
            )
        )

        self.play(
            Create(zoomed_camera_frame)
        )
        self.activate_zooming()
        self.play(
            self.get_zoomed_display_pop_out_animation()
        )
        self.wait()

        self.play(
            FadeIn(bcc_dot),
            FadeIn(fcc_dot)
        )
        self.wait()

        self.play(
            LaggedStart(
                Write(delta_gibbs_label),
                Write(delta_gibbs),
                Write(delta_gibbs_label_unit),
                lag_ratio=0.2
            )
        )
        self.wait()

        self.play(
            x_tracker.animate.set_value(intersect_x)
        )
        self.wait()

        self.play(
            x_tracker.animate.set_value(x_max)
        )
        self.wait()
