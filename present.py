import numpy as np
from manim import *
from rigid_mechanics_manim_physics import *

DOT_TO_LINE = 200


class IntroSlide(Scene):
    def construct(self):
        # title
        title = Text("Using the Time Cone Method to\nModel Competitive Phase Transformation", font_size=50)
        title.shift(UP * 1.2)

        # umich_mse_logo
        logo = ImageMobject("./figure/umich_mse_logo.png")
        logo.scale(1.2)
        logo.align_to(title, LEFT).shift(DOWN * 2)

        # name and affiliation
        name = Text("Peiwen Ren", font_size=40, slant=ITALIC)
        affiliation = Text("Sun Research Group\nUniversity of Michigan\nMaterials Science and Engineering\n09/23/2021",
                           font_size=30)
        name.next_to(logo, RIGHT).align_to(logo, UP)
        affiliation.next_to(name, DOWN).align_to(name, LEFT)

        # page_number
        page_1 = Tex("1", font_size=40)
        page_1.to_edge(DR)

        self.add(title, logo, name, affiliation, page_1)
        self.wait(1)
        self.play(FadeOut(logo), FadeOut(name), FadeOut(affiliation), title.animate.shift(DOWN * 1.2), run_time=2)
        self.wait(1)
        self.play(title.animate.set_color_by_t2c(t2c={"[8:23]": BLUE, "[28:-1]": RED}), run_time=1)
        self.wait()


class PhysicalPic(Scene):
    def construct(self):
        rectangle = ScreenRectangle()
        rectangle_mask = ScreenRectangle(stroke_width=200, stroke_color=BLACK, height=6, aspect_ratio=4 * 20.5 / 9 / 6)

        circles = [Circle(color=YELLOW, fill_color=YELLOW, fill_opacity=1, radius=1E-6, arc_center=shift_vec)
                   for shift_vec in [UR+0.5*RIGHT, 2*UL, -0.5*LEFT, DR+RIGHT, DL+LEFT]]

        phase_fraction_axes = Axes(
            x_range=[0, 1000, 100],
            y_range=[0, 1.2, 0.2],
            x_length=7,
            y_length=3,
            y_axis_config={
                "numbers_to_include": np.arange(0, 1 + 0.2, 0.2)
            }
        )

        phase_fraction_axes.to_edge(DOWN, buff=0.01)

        phase_fraction_x_label = phase_fraction_axes.get_x_axis_label(Tex("$t$"))
        phase_fraction_y_label = phase_fraction_axes.get_y_axis_label(
            Tex("$f_{\\text{transformed}}$").scale(0.7).rotate(PI / 2),
            edge=LEFT, direction=LEFT, buff=0.4
        )

        fraction_curve = phase_fraction_axes.plot(
            lambda t: 1 - np.exp(
                -np.pi / 3 * 0.012 * (0.001 ** 3) * (t ** 4)),
            color=YELLOW
        )

        jmak_formula_1 = Tex(
            "$f_{\\text{transformed}}$", color=YELLOW
        )
        jmak_formula_2 = Tex(
            "$=1-e^{-\\frac{\\pi}{3}J\\dot{R}^3t^4}$", color=YELLOW
        )

        jmak_formula_1.shift(4.5 * RIGHT + UP)
        jmak_formula_2.align_to(jmak_formula_1, LEFT).next_to(jmak_formula_1, DOWN)

        jmak_formula = VGroup(jmak_formula_1, jmak_formula_2)
        jmak_label = Tex(
            "JMAK Equation"
        ).next_to(jmak_formula, UP)

        nucleation_rate = Tex(
            "$J:\\,$", "nucleation rate", font_size=40
        ).next_to(jmak_formula, DOWN, buff=LARGE_BUFF)
        growth_rate = Tex(
            "$\\dot{R}:\\,$", "growth rate", font_size=40
        ).next_to(nucleation_rate, DOWN).align_to(nucleation_rate, LEFT)

        self.add(*circles)
        self.add_foreground_mobjects(rectangle, rectangle_mask)
        self.play(
            LaggedStart(
                *[circle.animate(rate_func=linear, run_time=5).scale(3E6) for circle in circles],
                lag_ratio=0.1
            )
        )
        self.wait()

        self.play(
            LaggedStart(
                *[circle.animate(rate_func=linear, run_time=1).scale(1/3E6) for circle in circles[::-1]],
                lag_ratio=0.1
            )
        )
        self.wait()

        self.play(
            *[obj.animate.shift(2*UP) for obj in [*circles, rectangle, rectangle_mask]]
        )
        self.wait()

        self.add_foreground_mobjects(phase_fraction_axes, phase_fraction_x_label, phase_fraction_y_label)
        self.play(Create(phase_fraction_axes), Write(phase_fraction_x_label), Write(phase_fraction_y_label),
                  run_time=2)
        self.wait()

        self.add_foreground_mobjects(fraction_curve)
        self.play(
            LaggedStart(
                *[circle.animate(rate_func=linear, run_time=5).scale(2.6E6) for circle in circles],
                lag_ratio=0.1
            ),
            Create(fraction_curve, rate_func=linear, run_time=6)
        )
        self.wait()

        self.play(
            *[obj.animate.shift(2*LEFT) for obj in self.mobjects]
        )
        self.wait()

        self.add_foreground_mobjects(jmak_label)
        self.play(Write(jmak_label))
        self.wait()

        self.add_foreground_mobjects(jmak_formula_1, jmak_formula_2)
        self.play(
            LaggedStart(
                Write(jmak_formula_1), Write(jmak_formula_2),
                lag_ratio=1
            )
        )
        self.wait()

        self.add_foreground_mobjects(nucleation_rate, growth_rate)
        self.play(
            LaggedStart(
                Write(nucleation_rate),
                Write(growth_rate),
                lag_ratio=1.1
            )
        )
        self.wait()


class Motivation(Scene):
    def construct(self):
        slide_title = Tex(
            "Motivation", font_size=40
        )
        slide_title.to_corner(UL)

        page_2 = Tex("$2$", font_size=40)
        page_2.to_edge(DR)

        self.add(slide_title, page_2)
        self.wait()

        motivation_statement = Tex(
            "In predictive synthesis, we want to know what is\\\\the optimal processing pathway\\\\to transition from a precursor phase to a desired phase?",
            font_size=45
        )
        motivation_statement.shift(UP * 2)

        self.play(Write(motivation_statement), run_time=3)
        self.wait()

        ttt_diagram = ImageMobject("./figure/ttt_diagram.png")
        ttt_diagram.scale(1.2)
        ttt_diagram.to_edge(RIGHT, buff=0.6).shift(DOWN * 1.7)
        ttt_diagram_label = Tex(
            "Time-Temperature-Transition Diagram", font_size=30
        )
        ttt_diagram_label.next_to(ttt_diagram, UP)
        ttt_diagram_line = Line(
            LEFT * 1.5, RIGHT * 1.5
        ).set_color(RED).move_to(ttt_diagram.get_center()).shift(UP * 0.09)

        phase_fraction = ImageMobject("./figure/phase_fraction.png")
        phase_fraction.scale(1.2)
        phase_fraction.shift(DOWN * 1.7).next_to(ttt_diagram, LEFT, buff=0.8)
        phase_fraction_label = Tex(
            "(Volume) Phase Fraction Transformed", font_size=30
        )
        phase_fraction_label.next_to(phase_fraction, UP).shift(LEFT * 0.7)

        question_mark = Tex(
            "$?$", font_size=50
        )
        question_mark.shift(DOWN * 1.7).next_to(phase_fraction, LEFT, buff=2)

        arrow_1 = Arrow(
            start=ttt_diagram.get_left(),
            end=phase_fraction.get_right()
        )
        arrow_1.move_to(
            (ttt_diagram.get_left() + phase_fraction.get_right()) / 2
        )

        arrow_2 = arrow_1.copy()
        arrow_2.move_to(
            (phase_fraction.get_left() + question_mark.get_right()) / 2 + (
                    phase_fraction.get_left() - question_mark.get_right()) / 4
        )

        time_cone_method_label = Tex(
            "Time Cone Method", font_size=30, color=BLUE
        )
        time_cone_method_label.move_to(question_mark.get_center())

        self.play(FadeIn(ttt_diagram), Write(ttt_diagram_label))
        self.wait()
        self.play(Create(ttt_diagram_line))
        self.wait()

        self.play(Create(arrow_1))
        self.play(FadeIn(phase_fraction), Write(phase_fraction_label))
        self.wait()

        self.play(Create(arrow_2))
        self.play(Write(question_mark))
        self.wait()
        self.play(ReplacementTransform(question_mark, time_cone_method_label))
        self.wait()

        self.play(arrow_2.animate.rotate(PI))
        self.play(arrow_1.animate.rotate(PI))
        self.wait()


class TimeCone1D2D(ThreeDScene):
    def construct(self):
        slide_title = Tex(
            "What is a time cone?", font_size=40
        )
        slide_title.to_corner(UL)

        page_3 = Tex("$3$", font_size=40)
        page_3.to_edge(DR)

        self.add_fixed_in_frame_mobjects(slide_title, page_3)
        self.wait()

        axes = ThreeDAxes(
            x_range=[-7, 7],
            y_range=[-7, 7],
            z_range=[-7, 7],
            x_length=7,
            y_length=7,
            z_length=7
        )

        axes_scale_factor = 1.7
        axes.scale(axes_scale_factor)
        axes.shift(2 * DOWN + 2 * LEFT)

        x_label = axes.get_x_axis_label(Tex("$x$"))
        x_label.add_updater(
            lambda x: x.move_to(
                axes.get_x_axis_label(Tex("$x$")).get_center())
        )
        y_label = axes.get_y_axis_label(Tex("$t$"))
        y_label.add_updater(
            lambda x: x.move_to(
                axes.get_y_axis_label(Tex("$t$")).get_center())
        )

        line = Line(start=axes.c2p(3 - np.interp(5, [0, 5], [0, 1]), 2, 0),
                    end=axes.c2p(3 + np.interp(5, [0, 5], [0, 1]), 2, 0))
        dot = Dot(radius=line.stroke_width / DOT_TO_LINE)
        dot.move_to(line.get_center())

        growth_cone_graph_left = axes.get_graph(
            lambda x: -5 * (x - 3), x_range=[2, 3]
        )
        growth_cone_graph_right = axes.get_graph(
            lambda x: 5 * (x - 3), x_range=[3, 4]
        )
        growth_cones_1d = VGroup(growth_cone_graph_left, growth_cone_graph_right)

        line_1 = Line(start=axes.c2p(3 - np.interp(1, [0, 5], [0, 1]), 1, 0),
                      end=axes.c2p(3 + np.interp(1, [0, 5], [0, 1]), 1, 0))
        line_2 = Line(start=axes.c2p(3 - np.interp(2, [0, 5], [0, 1]), 2, 0),
                      end=axes.c2p(3 + np.interp(2, [0, 5], [0, 1]), 2, 0))
        line_3 = Line(start=axes.c2p(3 - np.interp(3, [0, 5], [0, 1]), 3, 0),
                      end=axes.c2p(3 + np.interp(3, [0, 5], [0, 1]), 3, 0))
        line_4 = Line(start=axes.c2p(3 - np.interp(4, [0, 5], [0, 1]), 4, 0),
                      end=axes.c2p(3 + np.interp(4, [0, 5], [0, 1]), 4, 0))
        line_5 = Line(start=axes.c2p(3 - np.interp(5, [0, 5], [0, 1]), 5, 0),
                      end=axes.c2p(3 + np.interp(5, [0, 5], [0, 1]), 5, 0))
        lines_lst = [line_5, line_4, line_3, line_2, line_1]

        twod_group = VGroup(
            axes.get_x_axis(), axes.get_y_axis(), growth_cones_1d
        )

        self.play(FadeIn(dot))
        self.play(GrowFromCenter(line), rate_func=linear)
        self.remove(dot)
        self.wait()
        self.play(line.animate.move_to(line_5.get_center()), run_time=2)
        self.play(Create(axes.get_x_axis()))
        self.play(Write(x_label))
        self.play(Create(axes.get_y_axis()))
        self.play(Write(y_label))
        self.wait()

        self.play(LaggedStart(
            *[FadeIn(ind_line) for ind_line in lines_lst],
            lag_ratio=0.5
        ))
        self.remove(line)
        self.play(Create(growth_cones_1d))
        self.wait()
        self.play(LaggedStart(
            *[FadeOut(ind_line) for ind_line in lines_lst]
        ))
        self.wait()

        self.play(twod_group.animate.shift(2 * UP + 2 * RIGHT))
        self.play(twod_group.animate.scale(1 / axes_scale_factor))
        self.add(axes.get_z_axis().shift(2 * UP + 2 * RIGHT).scale(1 / axes_scale_factor),
                 axes.get_z_axis_label(Tex("$y$")))
        self.wait()

        self.move_camera(phi=60 * DEGREES, theta=-45 * DEGREES)
        self.wait()

        circle_1 = Circle(radius=1 / 10)
        circle_2 = Circle(radius=2 / 10)
        circle_3 = Circle(radius=3 / 10)
        circle_4 = Circle(radius=4 / 10)
        circle_5 = Circle(radius=5 / 10)
        circle_lst = [circle_1, circle_2, circle_3, circle_4, circle_5]
        circle_lst = [circle.shift(axes.c2p(3, num + 1, 0)).rotate(PI / 2, axis=X_AXIS)
                      for num, circle in enumerate(circle_lst)]

        cone = Cone(base_radius=1 / 2, height=5 / 2, direction=-Y_AXIS)
        cone.shift(axes.c2p(3, 0, 0))

        self.play(LaggedStart(
            *[GrowFromCenter(circle) for circle in circle_lst],
            lag_ratio=0.5
        ))
        self.play(Uncreate(growth_cones_1d))
        self.play(Create(cone), run_time=2)
        self.play(LaggedStart(
            *[FadeOut(circle) for circle in circle_lst]
        ))
        self.wait()


def create_single_growth_cone_with_shade(axes,
                                         growth_cone_tip_x, growth_cone_tip_y,
                                         growth_cone_base_radius, growth_cone_height,
                                         shade_opacity=0.5):
    growth_cone_graph_left = axes.get_graph(
        lambda x: -growth_cone_height / growth_cone_base_radius * (x - growth_cone_tip_x) + growth_cone_tip_y,
        x_range=[growth_cone_tip_x - growth_cone_base_radius, growth_cone_tip_x]
    )
    growth_cone_graph_right = axes.get_graph(
        lambda x: growth_cone_height / growth_cone_base_radius * (x - growth_cone_tip_x) + growth_cone_tip_y,
        x_range=[growth_cone_tip_x, growth_cone_tip_x + growth_cone_base_radius]
    )
    growth_cone = VGroup(growth_cone_graph_left, growth_cone_graph_right)

    growth_cone_shade = Polygon(growth_cone.get_bottom(),
                                axes.c2p(growth_cone_tip_x - growth_cone_base_radius,
                                         growth_cone_height),
                                axes.c2p(growth_cone_tip_x + growth_cone_base_radius,
                                         growth_cone_height))
    growth_cone_shade.set_fill(color=WHITE, opacity=shade_opacity)
    growth_cone_shade.set_stroke(color=WHITE, opacity=0)
    growth_cone_shade.add_updater(
        lambda x: x.move_to(growth_cone.get_center())
    )

    return growth_cone, growth_cone_shade


def create_triangle_path_points(top_point, left_point, right_point, num_segments=4):
    path_increment_vector_top_right = (right_point - top_point) / num_segments
    path_increment_vector_left_right = (right_point - left_point) / num_segments
    # path_points = [left_point, top_point,
    #                top_point + path_increment_vector_top_right,
    #                left_point + path_increment_vector_left_right,
    #                left_point + (2 * path_increment_vector_left_right),
    #                top_point + (2 * path_increment_vector_top_right),
    #                top_point + (3 * path_increment_vector_top_right),
    #                left_point + (3 * path_increment_vector_left_right),
    #                right_point]
    path_points = [left_point, top_point]
    for num in range(1, num_segments):
        if num % 2 != 0:
            path_points.append(top_point + num * path_increment_vector_top_right)
            path_points.append(left_point + num * path_increment_vector_left_right)
        else:
            path_points.append(left_point + num * path_increment_vector_left_right)
            path_points.append(top_point + num * path_increment_vector_top_right)
    path_points.append(right_point)

    return path_points


class TimeConeIntuition(Scene):
    def construct(self):
        slide_title = Tex(
            "How do we use the time cone?", font_size=40
        )
        slide_title.to_corner(UL).shift(UP * 0.3)

        page_4 = Tex("$4$", font_size=40)
        page_4.to_edge(DR)

        self.add(slide_title, page_4)
        self.wait()

        axes = Axes(
            x_range=[0, 7],
            y_range=[0, 7],
            x_length=8,
            y_length=6
        )

        x_label = axes.get_x_axis_label(Tex("$x$"))
        y_label = axes.get_y_axis_label(Tex("$t$"))

        growth_cone, growth_cone_shade = create_single_growth_cone_with_shade(axes, 3.5, 0, 1.5, 6)

        point_of_interest = Dot()
        point_of_interest_x, point_of_interest_y = 3.5, 3.5
        point_of_interest.move_to(axes.c2p(point_of_interest_x, point_of_interest_y))
        point_of_interest_label = Tex("$(x, t)$")
        point_of_interest_label.next_to(point_of_interest, RIGHT)

        self.play(Create(axes))
        self.play(Write(x_label), Write(y_label))
        self.wait()

        self.play(FadeIn(point_of_interest))
        self.play(Write(point_of_interest_label))
        self.wait()
        self.play(Unwrite(point_of_interest_label))
        self.wait()

        self.play(Create(growth_cone))
        self.wait()
        self.play(FadeIn(growth_cone_shade))
        self.wait()

        self.play(growth_cone.animate.shift(LEFT * 1.5))
        self.wait()
        self.play(growth_cone.animate.shift(RIGHT * 3))
        self.wait()
        self.play(growth_cone.animate.shift(LEFT * 0.5))
        self.wait()
        self.play(growth_cone.animate.shift(LEFT * 0.5))
        self.wait()
        self.play(FadeOut(growth_cone), FadeOut(growth_cone_shade))
        self.wait()

        growth_cone_lst = VGroup()
        growth_cone_tip_y = 0
        growth_cone_base_radius = 1.5
        growth_cone_height = 6
        for x_tip in np.arange(0, 7.5, 0.5):
            ind_growth_cone, ind_growth_cone_shade = create_single_growth_cone_with_shade(
                axes, x_tip, growth_cone_tip_y, growth_cone_base_radius, growth_cone_height, shade_opacity=0.2
            )
            growth_cone_lst.add(VGroup(ind_growth_cone, ind_growth_cone_shade))
        self.play(
            LaggedStart(
                *[FadeIn(growth_cone) for growth_cone in growth_cone_lst],
                lag_ratio=0.5
            )
        )
        self.wait()
        self.play(
            *[growth_cone_lst[_][0].animate.set_color(YELLOW) for _ in [6, 7, 8]],
            *[growth_cone_lst[_][1].animate.set_fill(YELLOW, opacity=0.2) for _ in [6, 7, 8]],
        )
        self.wait()

        nucleation_region_top = point_of_interest.get_center()
        nucleation_region_left = axes.c2p(point_of_interest_x -
                                          (point_of_interest_y / (growth_cone_height / growth_cone_base_radius)), 0)
        nucleation_region_right = axes.c2p(point_of_interest_x +
                                           (point_of_interest_y / (growth_cone_height / growth_cone_base_radius)), 0)
        nucleation_region_line_1 = DashedLine(start=nucleation_region_top,
                                              end=nucleation_region_left)
        nucleation_region_line_2 = DashedLine(start=nucleation_region_left,
                                              end=nucleation_region_right)
        nucleation_region_line_3 = DashedLine(start=nucleation_region_right,
                                              end=nucleation_region_top)
        nucleation_region = VGroup(nucleation_region_line_1, nucleation_region_line_2, nucleation_region_line_3)
        nucleation_region.set_stroke(color=BLUE)
        self.play(Create(nucleation_region))
        self.wait()

        self.play(FadeOut(growth_cone_lst))
        self.wait()

        path_correction_vector = growth_cone.get_center() - growth_cone.get_bottom()
        growth_cone.move_to(nucleation_region_left + path_correction_vector)
        path = VMobject()
        path_points = create_triangle_path_points(nucleation_region_top,
                                                  nucleation_region_left,
                                                  nucleation_region_right,
                                                  num_segments=6)
        path.set_points_as_corners([point + path_correction_vector for point in path_points])

        self.play(FadeIn(growth_cone), FadeIn(growth_cone_shade))
        self.wait()
        self.play(MoveAlongPath(growth_cone, path, rate_func=linear), run_time=10)
        self.wait()


class ConnectTimeConeWithPhaseFraction(Scene):
    def construct(self):
        slide_title = Text("Let's connect time cone to phase fraction transformed", font_size=30)
        slide_title.to_edge(UL)

        page_5 = Tex("$5$", font_size=40)
        page_5.to_edge(DR)

        statement_1 = Tex("$\\text{the phase fraction transformed}=1-\\text{the phase fraction untransformed}$",
                          font_size=40)
        statement_1.shift(UP * 2)

        statement_2 = Tex("$\\text{the phase fraction untransformed}=P(\\text{no nucleation has occurred})$",
                          font_size=40)
        statement_2.next_to(statement_1, DOWN, buff=0.5)

        statement_3 = Tex("$\\text{the phase fraction transformed}=1-P(\\text{no nucleation has occurred})$",
                          font_size=40)

        statement_4 = Tex("$P(\\text{no nucleation has occurred})$",
                          font_size=40)
        statement_4.align_to(statement_3, RIGHT)

        next_slide_title = Tex(
            "$\\text{How to calculate } P(\\text{no nucleation has occurred})\\text{ from time cone?}$",
            font_size=40)
        next_slide_title.to_corner(UL)

        self.add(slide_title, page_5)
        self.wait()

        self.play(Write(statement_1))
        self.wait()

        self.play(Write(statement_2), statement_1.animate.set_fill(color=WHITE, opacity=0.4))
        self.wait()

        self.play(Write(statement_3), statement_2.animate.set_fill(color=WHITE, opacity=0.4))
        self.add(statement_4)
        self.wait()

        self.play(
            LaggedStart(
                FadeOut(slide_title, statement_1, statement_2, statement_3),
                statement_4.animate.move_to(ORIGIN).scale(1.2),
                lag_ratio=0.5
            ),
            run_time=3
        )
        self.wait()

        page_6 = Tex("$6$", font_size=40)
        page_6.to_edge(DR)

        self.play(
            Transform(statement_4, next_slide_title),
            Transform(page_5, page_6),
        )
        self.wait()


class ProbToPhaseFraction(SpaceScene):
    def construct(self):
        sampling_box_axes = Axes(
            x_range=[0, 10],
            y_range=[0, 10],
            x_length=6,
            y_length=6,
            y_axis_config={
                "numbers_to_include": np.arange(0, 10 + 1, 1)
            },
            tips=False
        )

        grid_dimension_x = sampling_box_axes.x_range[1]
        grid_dimension_y = sampling_box_axes.y_range[1]

        bounding_box_left = Line(
            start=sampling_box_axes.c2p(0, grid_dimension_y),
            end=sampling_box_axes.c2p(0, 0)
        )
        bounding_box_bottom = Line(
            start=sampling_box_axes.c2p(0, 0),
            end=sampling_box_axes.c2p(grid_dimension_x, 0)
        )
        bounding_box_right = Line(
            start=sampling_box_axes.c2p(grid_dimension_x, 0),
            end=sampling_box_axes.c2p(grid_dimension_x, grid_dimension_y)
        )
        bounding_box = VGroup(bounding_box_left, bounding_box_bottom, bounding_box_right)

        self.add(sampling_box_axes, bounding_box)
        self.wait()

        # add a 10*10 grid of random numbers in the range of [0, 1)
        rng = np.random.default_rng(0)  # use a seeded random number generator to ensure reproducibility
        random_num_grid = rng.uniform(size=(grid_dimension_x, grid_dimension_y))
        # for each random number, determine if it is greater than the probability threshold
        # if higher than the threshold, desiganate this point as transformed
        prob_untransformed = 0.2
        transformed_grid = (random_num_grid >= prob_untransformed)

        # create a grid of white dots at those points which are transformed
        transformed_dot_lst = [Dot(point=sampling_box_axes.c2p(row_index + 0.5, col_index + 0.5),
                                   radius=(sampling_box_axes.c2p(0.5, 0) - sampling_box_axes.c2p(0, 0))[0])
                               for row_index in range(grid_dimension_x) for col_index in range(grid_dimension_y)
                               if transformed_grid[row_index][col_index]]

        count_tracker = ValueTracker(0)
        count_tracker_label = always_redraw(
            lambda: Integer(count_tracker.get_value()).next_to(bounding_box, RIGHT, buff=1)
        )
        num_of_dots = len(transformed_dot_lst)
        count_perc_eqn = MathTex(
            r"\frac{}{100}",
            r"={num_dots}\%".format(num_dots=num_of_dots),
            arg_separator=""
        )
        count_perc_eqn.align_to(count_tracker_label, UL).shift(
            np.array([(count_tracker_label.get_top() - count_perc_eqn.submobjects[0].get_top())[0]*0.7,
                      0, 0])
        ).shift(
            np.array([0,
                      (count_perc_eqn.submobjects[0].get_center() - count_tracker_label.get_center())[1]*0.8,
                      0])
        )

        self.play(LaggedStart(
            *[FadeIn(dot) for dot in transformed_dot_lst],
            lag_ratio=0.01
        ))
        self.wait()
        self.play(Write(count_tracker_label))

        self.play(
            LaggedStart(
                *[ShowPassingFlash(Circle().surround(dot, buffer_factor=1.1).set_color(YELLOW), time_width=2)
                  for dot in transformed_dot_lst],
                lag_ratio=0.05, run_time=4),
            count_tracker.animate(run_time=4, rate_func=linear).set_value(num_of_dots)
        )
        self.wait()
        self.play(Write(count_perc_eqn))
        self.wait()

        transformed_fraction_fill_line = DashedLine(
            start=sampling_box_axes.c2p(0, 10 * (1 - prob_untransformed)),
            end=sampling_box_axes.c2p(grid_dimension_x, 10 * (1 - prob_untransformed)),
            stroke_color=YELLOW
        )
        self.play(Create(transformed_fraction_fill_line))
        self.wait()

        self.make_static_body(bounding_box)
        self.make_rigid_body(*transformed_dot_lst)
        self.wait(8)


# noinspection DuplicatedCode
class ConnectProbToTimeCone(Scene):
    def construct(self):
        slide_title = Tex(
            "$\\text{How to calculate } P(\\text{no nucleation has occurred})\\text{ from time cone?}$",
            font_size=40)
        slide_title.to_corner(UL)

        page_6 = Tex("$6$", font_size=40)
        page_6.to_edge(DR)

        poisson_assumption = Text(
            "Assumption: nucleation events are random and independent", font_size=35
        )
        poisson_assumption.shift(UP * 2)
        poisson_assumption.set_color_by_t2c({"[:10]": YELLOW})

        poisson_equation = MathTex(
            r"P(X=k)=\frac{\lambda^ke^{-\lambda}}{k!}"
        )

        poisson_equation_zero_event_1 = MathTex(
            r"P(X=0)=\frac{\lambda^0e^{-\lambda}}{0!}"
        )

        poisson_equation_zero_event_2 = MathTex(
            r"P(X=0)=", r"\frac{\lambda^0e^{-\lambda}}{0!}"
        )

        poisson_equation_zero_event_simplified = MathTex(
            r"P(X=0)=", r"e^{-\lambda}"
        )
        poisson_equation_zero_event_simplified.align_to(poisson_equation_zero_event_2, LEFT)

        poisson_equation_zero_event_expected_number = MathTex(
            r"P(X=0)=", r"e^{-\langle X\rangle}"
        )
        poisson_equation_zero_event_expected_number.align_to(poisson_equation_zero_event_simplified, LEFT)

        expect_num_nuclei = Tex(
            "$\\langle X\\rangle$", "$:\\text{ expected number of nucleation events}$",
            arg_separator=" ")
        expect_num_nuclei.next_to(
            poisson_equation_zero_event_expected_number, DOWN)

        expect_num_nuclei_new_symbol = Tex(
            "$\\langle X\\rangle$", "$\\equiv$", "$\\langle N\\rangle$",
            arg_separator=" ")
        expect_num_nuclei_new_symbol.move_to(expect_num_nuclei.get_center())

        poisson_equation_zero_event_expected_number_new_symbol = MathTex(
            r"P(X=0)=", r"e^{-\langle N\rangle}"
        )

        expect_num_nuclei_unit_analysis = Tex(
            "$\\langle N\\rangle$", "$=$",
            "$\\frac{\\text{\\# of nucleation events}}{\\text{unit volume}\\cdot\\text{unit time}}$",
            "$\\cdot$", "$\\text{unit volume}$", "$\\cdot$", "$\\text{unit time}$",
            arg_separator=" "
        )
        expect_num_nuclei_unit_analysis.move_to(expect_num_nuclei_new_symbol.get_center())

        expect_num_nuclei_unit_analysis_symbols = Tex(
            "$\\langle N\\rangle$", "$=$",
            "$J$", "$\\cdot$", "$V$", "$\\cdot$", "$t$",
            arg_separator=" "
        )
        expect_num_nuclei_unit_analysis_symbols.move_to(expect_num_nuclei_unit_analysis.get_center())

        expect_num_nuclei_jvt = Tex(
            "$\\langle N\\rangle$", " $=$",
            " $J$", " $\\cdot$", " $($", "$V$", "$t$", "$)$",
            arg_separator=""
        )
        expect_num_nuclei_jvt.move_to(expect_num_nuclei_unit_analysis_symbols.get_center())

        expect_num_nuclei_jomega = Tex(
            "$\\langle N\\rangle$", " $=$",
            " $J$", " $\\cdot$", " $||\\Omega_N||$",
            arg_separator=""
        )
        expect_num_nuclei_jomega.move_to(expect_num_nuclei_jvt.get_center())

        prob_no_nucleation = MathTex(
            r"P(X=0)=", r"e^{-J||\Omega_N||}"
        )
        # prob_no_nucleation.align_to(poisson_equation_zero_event_expected_number_new_symbol, LEFT)

        nucleation_rate_sym = expect_num_nuclei_unit_analysis.submobjects[2].copy()
        nucleation_rate_brace = BraceText(nucleation_rate_sym, "$J$")

        volume_sym = expect_num_nuclei_unit_analysis.submobjects[4].copy()
        volume_brace = BraceText(volume_sym, "$V$")

        time_sym = expect_num_nuclei_unit_analysis.submobjects[6].copy()
        time_brace = BraceText(time_sym, "$t$")

        self.add(slide_title, page_6)
        self.wait()

        self.play(Write(poisson_assumption))
        self.wait()

        self.play(Write(poisson_equation))
        self.wait()

        self.play(ReplacementTransform(poisson_equation, poisson_equation_zero_event_1))
        self.add(poisson_equation_zero_event_2)
        self.remove(poisson_equation_zero_event_1)
        self.wait()

        self.play(TransformMatchingTex(poisson_equation_zero_event_2, poisson_equation_zero_event_simplified))
        self.wait()

        self.play(ReplacementTransform(poisson_equation_zero_event_simplified,
                                       poisson_equation_zero_event_expected_number))
        self.wait()

        self.play(Write(expect_num_nuclei))
        self.wait()

        self.play(TransformMatchingTex(expect_num_nuclei, expect_num_nuclei_new_symbol))
        self.wait()

        self.play(LaggedStart(
            ReplacementTransform(poisson_equation_zero_event_expected_number,
                                 poisson_equation_zero_event_expected_number_new_symbol),
            TransformMatchingTex(expect_num_nuclei_new_symbol, expect_num_nuclei_unit_analysis),
            lag_ratio=1
        ))
        self.wait()

        self.play(LaggedStart(
            Write(nucleation_rate_brace),
            Write(volume_brace),
            Write(time_brace),
            lag_ratio=1
        ), run_time=2)
        self.wait()

        self.play(TransformMatchingShapes(VGroup(expect_num_nuclei_unit_analysis,
                                                 nucleation_rate_brace,
                                                 volume_brace,
                                                 time_brace),
                                          VGroup(expect_num_nuclei_unit_analysis_symbols)))
        self.wait()

        self.play(TransformMatchingTex(expect_num_nuclei_unit_analysis_symbols,
                                       expect_num_nuclei_jvt))
        self.wait()

        self.play(TransformMatchingTex(expect_num_nuclei_jvt,
                                       expect_num_nuclei_jomega))
        self.wait()

        self.play(
            LaggedStart(
                ReplacementTransform(poisson_equation_zero_event_expected_number_new_symbol,
                                     prob_no_nucleation),
                FadeOut(expect_num_nuclei_jomega),
                lag_ratio=1
            ))
        self.wait()
        self.play(prob_no_nucleation.animate.next_to(poisson_assumption, DOWN))
        self.wait()

        axes = Axes(
            x_range=[0, 7],
            y_range=[0, 7],
            x_length=8,
            y_length=6
        )

        axes.move_to(prob_no_nucleation.get_bottom(), aligned_edge=UP)
        axes.scale(0.6)

        x_label = axes.get_x_axis_label(Tex("$x$"))
        y_label = axes.get_y_axis_label(Tex("$t$"))

        point_of_interest = Dot()
        point_of_interest_x, point_of_interest_y = 3.5, 3.5
        point_of_interest.move_to(axes.c2p(point_of_interest_x, point_of_interest_y))

        growth_cone_height = 6
        growth_cone_base_radius = 1.5

        nucleation_region_top = point_of_interest.get_center()
        nucleation_region_left = axes.c2p(point_of_interest_x -
                                          (point_of_interest_y / (growth_cone_height / growth_cone_base_radius)), 0)
        nucleation_region_right = axes.c2p(point_of_interest_x +
                                           (point_of_interest_y / (growth_cone_height / growth_cone_base_radius)), 0)
        nucleation_region_line_1 = DashedLine(start=nucleation_region_top,
                                              end=nucleation_region_left)
        nucleation_region_line_2 = DashedLine(start=nucleation_region_left,
                                              end=nucleation_region_right)
        nucleation_region_line_3 = DashedLine(start=nucleation_region_right,
                                              end=nucleation_region_top)
        nucleation_region = VGroup(nucleation_region_line_1, nucleation_region_line_2, nucleation_region_line_3)
        nucleation_region.set_stroke(color=BLUE)
        nucleation_region_polygon = Polygon(nucleation_region_top,
                                            nucleation_region_left,
                                            nucleation_region_right)
        nucleation_region_polygon.set_stroke(opacity=0)
        nucleation_region_polygon.set_fill(BLUE, opacity=0.2)

        self.play(Create(axes))
        self.play(Create(x_label), Create(y_label))
        self.play(FadeIn(point_of_interest))
        self.play(Create(nucleation_region), FadeIn(nucleation_region_polygon), run_time=2)
        self.wait()


class CalcTimeConeVolume(ThreeDScene):
    def construct(self):
        slide_title = Tex(
            "$\\text{Calculate } ||\\Omega_N|| \\text{ in 1D, 2D, 3D space-time}$",
            font_size=40)
        slide_title.to_corner(UL)

        page_7 = Tex("$7$", font_size=40)
        page_7.to_edge(DR)

        axes = ThreeDAxes(
            x_range=[0, 7],
            y_range=[0, 7],
            z_range=[0, 7],
            x_length=8,
            y_length=8,
            z_length=6
        )

        axes.shift(-Z_AXIS * 3.3)
        axes.scale(0.6)

        x_axis = axes.get_x_axis()
        y_axis = axes.get_y_axis()
        z_axis = axes.get_z_axis()

        axes.x_axis = x_axis.rotate(PI / 2, axis=X_AXIS)

        x_label = axes.get_x_axis_label(Tex("$x$")).rotate(PI / 2, axis=X_AXIS)
        y_label = axes.get_y_axis_label(Tex("$y$")).rotate(PI / 2, axis=X_AXIS)
        z_label = axes.get_z_axis_label(Tex("$t$"))

        point_of_interest_x, point_of_interest_y, point_of_interest_z = 3.5, 0, 3.5
        point_of_interest = Dot3D(point=axes.c2p(point_of_interest_x, point_of_interest_y, point_of_interest_z))

        growth_cone_height = 6
        growth_cone_base_radius = 1.5

        nucleation_region_top = point_of_interest.get_center()
        nucleation_region_left = axes.c2p(point_of_interest_x -
                                          (point_of_interest_z / (growth_cone_height / growth_cone_base_radius)), 0, 0)
        nucleation_region_right = axes.c2p(point_of_interest_x +
                                           (point_of_interest_z / (growth_cone_height / growth_cone_base_radius)), 0, 0)
        nucleation_region_bottom_mid = axes.c2p(point_of_interest_x, 0, 0)
        nucleation_region_line_1 = DashedLine(start=nucleation_region_top,
                                              end=nucleation_region_left)
        nucleation_region_line_2 = DashedLine(start=nucleation_region_left,
                                              end=nucleation_region_right)
        nucleation_region_line_3 = DashedLine(start=nucleation_region_right,
                                              end=nucleation_region_top)
        nucleation_region = VGroup(nucleation_region_line_1, nucleation_region_line_2, nucleation_region_line_3)
        nucleation_region.set_stroke(color=BLUE)
        # nucleation_region_polygon = Polygon(nucleation_region_top,
        #                                     nucleation_region_left,
        #                                     nucleation_region_right)
        # nucleation_region_polygon.set_stroke(opacity=0)
        # nucleation_region_polygon.set_fill(BLUE, opacity=0.2)

        height_2d = DashedLine(start=nucleation_region_top, end=nucleation_region_bottom_mid)
        # TODO: why is the height_2d_brace not rendered
        height_2d_brace = Brace(height_2d, direction=X_AXIS).rotate(
            PI / 2, axis=X_AXIS, about_point=(nucleation_region_bottom_mid + nucleation_region_top) / 2
        )
        height_2d_brace_label = height_2d_brace.get_tex("t").rotate(
            PI / 2, axis=X_AXIS, about_point=(nucleation_region_bottom_mid + nucleation_region_top) / 2
        )
        base_2d_brace = BraceText(Line(start=nucleation_region_bottom_mid,
                                       end=nucleation_region_right), "$\\dot{R}t$").rotate(PI / 2, axis=X_AXIS,
                                                                                           about_point=nucleation_region_bottom_mid)

        self.add_fixed_in_frame_mobjects(slide_title, page_7)
        self.wait()

        self.set_camera_orientation(phi=90 * DEGREES, theta=-90 * DEGREES)
        self.play(Create(x_axis), Create(z_axis))
        self.play(Create(x_label), Create(z_label))
        self.play(FadeIn(point_of_interest))
        self.play(Create(nucleation_region), run_time=2)
        self.wait()

        self.play(Create(height_2d))
        self.wait()
        self.play(Write(height_2d_brace))
        self.play(Write(height_2d_brace_label))
        self.wait()
        self.play(Write(base_2d_brace))
        self.wait()

        time_cone_1d_area_eqn = Tex(
            "$||\\Omega_N||_{\\text{1D}}$", " $=$",
            " $\\frac{1}{2}$", " $\\cdot$",
            " $2$", "$\\dot{R}$", "$t$",
            " $\\cdot$", " $t$"
        )
        time_cone_1d_area_eqn.shift(UP * 2.8 + RIGHT * 4)
        # time_cone_1d_area_eqn.rotate(PI/2, axis=X_AXIS)

        time_cone_1d_area_eqn_symbols = Tex(
            "$||\\Omega_N||_{\\text{1D}}$", " $=$",
            " $\\dot{R}$", "$t^2$"
        )
        time_cone_1d_area_eqn_symbols.align_to(time_cone_1d_area_eqn, LEFT)
        time_cone_1d_area_eqn_symbols.align_to(time_cone_1d_area_eqn, UP)

        self.add_fixed_in_frame_mobjects(time_cone_1d_area_eqn)
        self.play(Write(time_cone_1d_area_eqn))
        self.wait()

        self.play(ReplacementTransform(time_cone_1d_area_eqn, time_cone_1d_area_eqn_symbols))
        self.add_fixed_in_frame_mobjects(time_cone_1d_area_eqn_symbols)
        self.wait()

        # convert from quasi-2d to 3d
        self.move_camera(phi=85 * DEGREES, theta=-85 * DEGREES)
        self.wait()
        self.play(Create(y_axis))
        self.play(Create(y_label))
        self.wait()

        cone_2d_space_time = Cone(base_radius=nucleation_region_bottom_mid[0] - nucleation_region_left[0],
                                  height=nucleation_region_top[-1] - nucleation_region_bottom_mid[-1],
                                  direction=Z_AXIS)
        cone_2d_space_time.shift(axes.c2p(point_of_interest_x,
                                          point_of_interest_y,
                                          point_of_interest_z))

        self.play(
            LaggedStart(
                Rotate(
                    nucleation_region, angle=PI, axis=Z_AXIS
                ),
                FadeIn(cone_2d_space_time),
                lag_ratio=0.1
            ),
            run_time=2
        )
        self.wait()

        time_cone_2d_area_eqn = Tex(
            "$||\\Omega_N||_{\\text{2D}}$", " $=$",
            " $\\frac{1}{3}$", " $\\cdot$",
            " $\\pi$", "$(\\dot{R}t)^2$", " $\\cdot$", " $t$"
        )
        time_cone_2d_area_eqn.shift(UP * 2.8 + RIGHT * 4)
        time_cone_2d_area_eqn.next_to(time_cone_1d_area_eqn_symbols, DOWN)
        time_cone_2d_area_eqn.align_to(time_cone_1d_area_eqn_symbols, LEFT)

        self.add_fixed_in_frame_mobjects(time_cone_2d_area_eqn)
        self.play(Write(time_cone_2d_area_eqn))
        self.wait()

        time_cone_2d_area_eqn_symbols = Tex(
            "$||\\Omega_N||_{\\text{2D}}$", " $=$",
            " $\\frac{\\pi}{3}$",
            "$\\dot{R}^2$", "$t^3$"
        )
        time_cone_2d_area_eqn_symbols.rotate(
            PI / 2,
            axis=(axes.c2p(2, 0, 0) - axes.c2p(1, 0, 0)) / np.linalg.norm(axes.c2p(2, 0, 0) - axes.c2p(1, 0, 0)),
            about_point=time_cone_2d_area_eqn_symbols.get_center())
        time_cone_2d_area_eqn_symbols.align_to(time_cone_2d_area_eqn, LEFT)
        time_cone_2d_area_eqn_symbols.align_to(time_cone_2d_area_eqn, UP)
        time_cone_2d_area_eqn_symbols.shift(np.array((0.0, 0.0, 2.0)))

        self.play(ReplacementTransform(time_cone_2d_area_eqn, time_cone_2d_area_eqn_symbols))
        # self.add_fixed_in_frame_mobjects(time_cone_2d_area_eqn_symbols)
        self.wait()

        circle_lst = []
        for circle_height in np.linspace(0, point_of_interest_z, 4, endpoint=False):
            base_rad = np.interp(
                circle_height,
                [0, point_of_interest_z],
                [nucleation_region_bottom_mid[0] - nucleation_region_left[0], 0]
            )
            circle_to_add = Circle(radius=base_rad)
            circle_to_add.set_fill(color=RED, opacity=0.6)
            circle_lst.append(circle_to_add)

        circle_lst = [circle.shift(axes.c2p(point_of_interest_x, 0, num))
                      for num, circle in enumerate(circle_lst)]

        self.play(FadeIn(*circle_lst))
        self.wait()

        alter_time_cone_2d_area_eqn = MathTex(
            r"\int_{t'=0}^{t'=t}\pi(\dot{R}t')^2dt'",
            r"&=\pi\dot{R}^2\frac{{t'}^3}{3}\bigg|_{t'=0}^{t'=t}\\",
            r"&=\frac{\pi}{3}\dot{R}^2t^3",
            arg_separator=""
        )
        alter_time_cone_2d_area_eqn.shift(UP * 0.3 + RIGHT)

        self.add_fixed_in_frame_mobjects(alter_time_cone_2d_area_eqn)
        self.play(Write(alter_time_cone_2d_area_eqn))
        self.wait()

        alter_time_cone_3d_area_eqn = MathTex(
            r"\int_{t'=0}^{t'=t}\frac{4}{3}\pi(\dot{R}t')^3dt'",
            r"&=\frac{4\pi}{3}\dot{R}^3\frac{{t'}^4}{4}\bigg|_{t'=0}^{t'=t}\\",
            r"&=\frac{\pi}{3}\dot{R}^3t^4",
            arg_separator=""
        )
        alter_time_cone_3d_area_eqn.rotate(
            PI / 2, axis=X_AXIS, about_point=alter_time_cone_3d_area_eqn.get_center()
        )
        alter_time_cone_3d_area_eqn.align_to(alter_time_cone_2d_area_eqn, LEFT)
        alter_time_cone_3d_area_eqn.align_to(alter_time_cone_2d_area_eqn, UP)

        self.play(TransformMatchingShapes(alter_time_cone_2d_area_eqn,
                                          alter_time_cone_3d_area_eqn))
        self.wait()

        time_cone_3d_area_eqn = Tex(
            "$||\\Omega_N||_{\\text{3D}}=\\frac{\\pi}{3}\\dot{R}^3t^4$"
        )
        time_cone_3d_area_eqn.next_to(time_cone_2d_area_eqn_symbols, DOWN, buff=0.8)
        time_cone_3d_area_eqn.align_to(time_cone_2d_area_eqn_symbols, LEFT)
        # time_cone_3d_area_eqn.shift(UP * 2.2 + RIGHT * 4)

        self.add_fixed_in_frame_mobjects(time_cone_3d_area_eqn)
        self.play(
            LaggedStart(
                FadeOut(alter_time_cone_3d_area_eqn),
                Write(time_cone_3d_area_eqn),
                lag_ratio=0.8
            ),
            run_time=2
        )
        self.wait()


class SinglePhaseFraction(Scene):
    def construct(self):
        slide_title = Tex(
            "Phase fraction transformed as a function of time (single-phase)",
            font_size=40)
        slide_title.to_corner(UL)

        page_8 = Tex("$8$", font_size=40)
        page_8.to_edge(DR)

        statement_1 = Tex(
            "the phase fraction transformed", "$=$",
            "$1-$", "the phase fraction untransformed",
            font_size=40,
            arg_separator=" "
        )
        statement_1.shift(UP * 2)

        statement_2 = Tex(
            "the phase fraction transformed", "$=$",
            "$1-$", "$P(\\text{no nucleation has occurred})$",
            font_size=40,
            arg_separator=" "
        )
        statement_2.move_to(statement_1.get_center())

        statement_3 = Tex(
            "the phase fraction transformed", "$=$",
            "$1-$", "$e^{-J||\\Omega_N||}$",
            font_size=40,
            arg_separator=" "
        )
        statement_3.move_to(statement_2.get_center()).align_to(statement_2, LEFT)

        statement_4 = Tex(
            "the phase fraction transformed", "$=$",
            "$\\begin{cases} 1-\\exp(-J\\dot{R}t^2) & \\text{1D}\\\\ 1-\\exp(-\\frac{\\pi}{3}J\\dot{R}^2t^3) & \\text{2D}\\\\ 1-\\exp(-\\frac{\\pi}{3}J\\dot{R}^3t^4) & \\text{3D}\\end{cases}$",
        )

        jmak_box = Rectangle(
            width=6, height=0.8
        )
        jmak_box.set_stroke(color=YELLOW)
        jmak_box.move_to(statement_4.get_corner(direction=DR)).shift(
            (statement_4.get_top()[1] - statement_4.get_bottom()[1]) / 6 + LEFT * 3)
        jmak_box_brace = BraceText(jmak_box, "JMAK Equation")

        phase_fraction_all_cases = VGroup(
            statement_4, jmak_box, jmak_box_brace
        )

        axes = Axes(
            x_range=[0, 6, 1],
            y_range=[0, 1.2, 0.2],
            x_length=8,
            y_length=4,
            x_axis_config={
                "numbers_to_include": np.arange(0, 5 + 1, 1)
            },
            y_axis_config={
                "numbers_to_include": np.arange(0, 1 + 0.2, 0.2)
            }
        )

        axes.shift(DOWN * 1.4)
        x_label = axes.get_x_axis_label(Tex("$t$"))
        y_label = axes.get_y_axis_label(Tex("$f_{\\text{transformed}}$"))

        nucleation_rate = ValueTracker(0.24)
        growth_rate = ValueTracker(0.76)
        t_valuetracker = ValueTracker(0)
        fraction_curve = always_redraw(
            lambda:
            axes.get_graph(
                lambda t: 1 - np.exp(
                    -np.pi / 3 * nucleation_rate.get_value() * (growth_rate.get_value() ** 3) * (t ** 4)),
                x_range=[0, t_valuetracker.get_value()],
                color=YELLOW
            )
        )

        nucleation_rate_text = Tex(
            "$J$: "
        ).next_to(axes, RIGHT, buff=0.5).shift(UP)

        nucleation_rate_value = always_redraw(
            lambda: DecimalNumber(num_decimal_places=2).set_value(
                nucleation_rate.get_value()).next_to(nucleation_rate_text)
        )

        growth_rate_text = Tex(
            "$\\dot{R}$: "
        ).next_to(nucleation_rate_text, DOWN, buff=0.3)

        growth_rate_value = always_redraw(
            lambda: DecimalNumber(num_decimal_places=2).set_value(
                growth_rate.get_value()).next_to(growth_rate_text).align_to(
                growth_rate_text.get_bottom(), DOWN).align_to(nucleation_rate_value, LEFT)
        )

        j_decrease_arrow = Arrow(
            start=nucleation_rate_text.get_top(),
            end=nucleation_rate_text.get_bottom(),
            stroke_width=2
        ).next_to(nucleation_rate_text, LEFT, buff=0.1)

        r_increase_arrow = Arrow(
            start=nucleation_rate_text.get_bottom(),
            end=nucleation_rate_text.get_top(),
            stroke_width=2
        ).next_to(growth_rate_text, LEFT, buff=0.1).align_to(
            j_decrease_arrow, LEFT
        ).align_to(
            growth_rate_text, DOWN
        )

        self.add(slide_title, page_8)
        self.wait()

        self.play(Write(statement_1))
        self.wait()

        self.play(TransformMatchingTex(statement_1, statement_2))
        self.wait()

        self.play(ReplacementTransform(statement_2, statement_3))
        self.wait()

        self.play(TransformMatchingTex(statement_3, statement_4))
        self.wait()

        self.play(Create(jmak_box))
        self.wait()
        self.play(Write(jmak_box_brace))
        self.wait()
        self.play(phase_fraction_all_cases.animate.shift(UP * 2).scale(0.8), run_time=2)
        self.wait()
        self.play(FadeOut(jmak_box_brace))
        self.wait()

        self.play(Create(axes))
        self.play(Write(x_label), Write(y_label))
        self.wait()

        self.add(fraction_curve)
        self.play(Write(VGroup(nucleation_rate_text, nucleation_rate_value,
                               growth_rate_text, growth_rate_value)))
        self.play(t_valuetracker.animate.set_value(6), run_time=5, rate_func=linear)
        self.wait()

        self.play(growth_rate.animate(run_time=5, rate_func=linear).set_value(1.5),
                  FadeIn(r_increase_arrow, run_time=0.5))
        self.play(FadeOut(r_increase_arrow, run_time=0.5))
        self.wait()

        self.play(nucleation_rate.animate(run_time=5, rate_func=linear).set_value(0.01),
                  FadeIn(j_decrease_arrow, run_time=0.5))
        self.play(FadeOut(j_decrease_arrow, run_time=0.5))
        self.wait()


class PhaseFractionToTTT(Scene):
    def construct(self):
        animation_duration = 5

        phase_fraction_axes = Axes(
            x_range=[0, 1000, 100],
            y_range=[0, 1.2, 0.2],
            x_length=8,
            y_length=2,
            y_axis_config={
                "numbers_to_include": np.arange(0, 1 + 0.2, 0.2)
            }
        )

        phase_fraction_axes.shift(LEFT * 0.5)
        phase_fraction_axes.to_edge(UP, buff=0.5)

        phase_fraction_x_label = phase_fraction_axes.get_x_axis_label(Tex("$t$"))
        phase_fraction_y_label = phase_fraction_axes.get_y_axis_label(
            Tex("$f_{\\text{transformed}}$").scale(0.7).rotate(PI / 2),
            edge=LEFT, direction=LEFT, buff=0.4
        )

        ttt_axes = Axes(
            x_range=[0, 1000, 100],
            y_range=[0, animation_duration + 2, 1],
            x_length=8,
            y_length=4,
            x_axis_config={
                "numbers_to_include": np.arange(0, 1000 + 100, 100)
            }
        )

        ttt_axes.shift(LEFT * 0.5)
        ttt_axes.to_edge(DOWN, buff=0.1)

        ttt_x_label = ttt_axes.get_x_axis_label(Tex("$t$"))
        ttt_y_label = ttt_axes.get_y_axis_label(
            Tex("$T$"),
            edge=LEFT, direction=LEFT, buff=0.4
        )

        nucleation_start_value, nucleation_end_value = 0.01, 0.00001
        growth_rate_start_value, growth_rate_end_value = 0.0008, 0.01

        # a proxy for temperature
        temp_proxy_tracker = ValueTracker(0)
        nucleation_rate = ValueTracker(nucleation_start_value).add_updater(
            lambda num: num.set_value(
                nucleation_start_value + (
                        nucleation_end_value - nucleation_start_value) / animation_duration * temp_proxy_tracker.get_value()
            ))
        growth_rate = ValueTracker(growth_rate_start_value).add_updater(
            lambda num: num.set_value(
                growth_rate_start_value + (
                        growth_rate_end_value - growth_rate_start_value) / animation_duration * temp_proxy_tracker.get_value()
            ))

        self.add(nucleation_rate, growth_rate)

        fraction_curve = always_redraw(
            lambda:
            phase_fraction_axes.get_graph(
                lambda t: 1 - np.exp(
                    -np.pi / 3 * nucleation_rate.get_value() * (growth_rate.get_value() ** 3) * (t ** 4))
            )
        )

        nucleation_rate_text = Tex(
            "$J$: "
        ).next_to(phase_fraction_axes, RIGHT, buff=0.4).shift(UP)

        nucleation_rate_value = always_redraw(
            lambda: DecimalNumber(num_decimal_places=1).set_value(
                nucleation_rate.get_value() * 1e4).next_to(nucleation_rate_text)
        )

        scale_factor_nucleation = Tex(
            "$\\times 10^{-4}$"
        )
        scale_factor_nucleation.add_updater(lambda x: x.next_to(nucleation_rate_value, RIGHT, buff=0.1).align_to(
            nucleation_rate_value, DOWN
        ))

        growth_rate_text = Tex(
            "$\\dot{R}$: "
        ).next_to(nucleation_rate_text, DOWN, buff=0.3)

        growth_rate_value = always_redraw(
            lambda: DecimalNumber(num_decimal_places=1).set_value(
                growth_rate.get_value() * 1e4).next_to(growth_rate_text).align_to(
                growth_rate_text.get_bottom(), DOWN).align_to(nucleation_rate_value, LEFT)
        )

        scale_factor_growth_rate = Tex(
            "$\\times 10^{-4}$"
        )
        scale_factor_growth_rate.add_updater(lambda x: x.next_to(growth_rate_value, RIGHT, buff=0.1).align_to(
            growth_rate_value, DOWN
        ))

        self.play(Create(phase_fraction_axes))
        self.play(Write(phase_fraction_x_label), Write(phase_fraction_y_label))
        self.wait()

        self.play(Create(ttt_axes))
        self.play(Write(ttt_x_label), Write(ttt_y_label))
        self.wait()

        self.play(Write(VGroup(nucleation_rate_text, nucleation_rate_value, scale_factor_nucleation,
                               growth_rate_text, growth_rate_value, scale_factor_growth_rate)))
        self.play(Create(fraction_curve), run_time=3)
        self.wait()

        # figure out the time it takes to reach 1%, 50%, 99% transformed
        # AHA: need to be more careful about object reference
        # fraction_transformed_lst = [0.01, 0.5, 0.99]
        # time_it_takes_lst = []
        # fraction_dot_lst = []
        # for fraction_transformed, dot_color in zip(fraction_transformed_lst, [BLUE, YELLOW, PURPLE]):
        #     time_it_takes = ValueTracker(np.power(
        #         - 3 * np.log(1 - fraction_transformed) / (
        #                     np.pi * nucleation_rate.get_value() * (growth_rate.get_value() ** 3)),
        #         0.25
        #     )).add_updater(lambda num: num.set_value(np.power(
        #         - 3 * np.log(1 - fraction_transformed) / (
        #                     np.pi * nucleation_rate.get_value() * (growth_rate.get_value() ** 3)),
        #         0.25
        #     )))
        #     time_it_takes_lst.append(time_it_takes)
        #     # fraction_dot_x = time_it_takes.get_value()
        #     # fraction_dot_y = fraction_transformed
        #     fraction_dot = always_redraw(
        #         lambda: Dot(point=phase_fraction_axes.c2p(time_it_takes.get_value(), fraction_transformed),
        #                     color=dot_color)
        #     )
        #     fraction_dot_lst.append(fraction_dot)
        #
        # self.add(*time_it_takes_lst)
        # self.play(Create(fraction_dot_lst[0]), Create(fraction_dot_lst[1]), Create(fraction_dot_lst[2]))
        # self.wait()

        time_it_takes_10_perc = ValueTracker(np.power(
            - 3 * np.log(1 - 0.1) / (np.pi * nucleation_rate.get_value() * (growth_rate.get_value() ** 3)), 0.25
        )).add_updater(lambda num: num.set_value(np.power(
            - 3 * np.log(1 - 0.1) / (np.pi * nucleation_rate.get_value() * (growth_rate.get_value() ** 3)), 0.25)))
        fraction_line_horizontal_10_perc = always_redraw(
            lambda: DashedLine(
                start=phase_fraction_axes.c2p(0, 0.1),
                end=phase_fraction_axes.c2p(time_it_takes_10_perc.get_value(), 0.1),
                color=BLUE
            )
        )
        fraction_line_vertical_10_perc = always_redraw(
            lambda: DashedLine(
                start=phase_fraction_axes.c2p(time_it_takes_10_perc.get_value(), 0.1),
                end=ttt_axes.c2p(time_it_takes_10_perc.get_value(), temp_proxy_tracker.get_value() + 1),
                color=BLUE
            )
        )
        fraction_dot_10_perc = always_redraw(
            lambda: Dot(point=phase_fraction_axes.c2p(time_it_takes_10_perc.get_value(), 0.1),
                        color=BLUE)
        )
        ttt_dot_10_perc = always_redraw(
            lambda: Dot(point=ttt_axes.c2p(time_it_takes_10_perc.get_value(), temp_proxy_tracker.get_value() + 1),
                        color=BLUE)
        )
        ttt_dot_10_perc_path = TracedPath(ttt_dot_10_perc.get_center, stroke_color=BLUE, stroke_width=4)

        self.add(time_it_takes_10_perc)
        self.play(Create(fraction_line_horizontal_10_perc))
        self.play(Create(fraction_dot_10_perc))
        self.play(Create(fraction_line_vertical_10_perc))
        self.play(Create(ttt_dot_10_perc))
        self.add(ttt_dot_10_perc_path)

        time_it_takes_50_perc = ValueTracker(np.power(
            - 3 * np.log(1 - 0.5) / (np.pi * nucleation_rate.get_value() * (growth_rate.get_value() ** 3)), 0.25
        )).add_updater(lambda num: num.set_value(np.power(
            - 3 * np.log(1 - 0.5) / (np.pi * nucleation_rate.get_value() * (growth_rate.get_value() ** 3)), 0.25)))
        fraction_line_horizontal_50_perc = always_redraw(
            lambda: DashedLine(
                start=phase_fraction_axes.c2p(0, 0.5),
                end=phase_fraction_axes.c2p(time_it_takes_50_perc.get_value(), 0.5),
                color=YELLOW
            )
        )
        fraction_line_vertical_50_perc = always_redraw(
            lambda: DashedLine(
                start=phase_fraction_axes.c2p(time_it_takes_50_perc.get_value(), 0.5),
                end=ttt_axes.c2p(time_it_takes_50_perc.get_value(), temp_proxy_tracker.get_value() + 1),
                color=YELLOW
            )
        )
        fraction_dot_50_perc = always_redraw(
            lambda: Dot(point=phase_fraction_axes.c2p(time_it_takes_50_perc.get_value(), 0.5),
                        color=YELLOW)
        )
        ttt_dot_50_perc = always_redraw(
            lambda: Dot(point=ttt_axes.c2p(time_it_takes_50_perc.get_value(), temp_proxy_tracker.get_value() + 1),
                        color=YELLOW)
        )
        ttt_dot_50_perc_path = TracedPath(ttt_dot_50_perc.get_center, stroke_color=YELLOW, stroke_width=4)

        self.add(time_it_takes_50_perc)
        self.play(Create(fraction_line_horizontal_50_perc))
        self.play(Create(fraction_dot_50_perc))
        self.play(Create(fraction_line_vertical_50_perc))
        self.play(Create(ttt_dot_50_perc))
        self.add(ttt_dot_50_perc_path)

        time_it_takes_99_perc = ValueTracker(np.power(
            - 3 * np.log(1 - 0.99) / (np.pi * nucleation_rate.get_value() * (growth_rate.get_value() ** 3)), 0.25
        )).add_updater(lambda num: num.set_value(np.power(
            - 3 * np.log(1 - 0.99) / (np.pi * nucleation_rate.get_value() * (growth_rate.get_value() ** 3)), 0.25)))
        fraction_line_horizontal_99_perc = always_redraw(
            lambda: DashedLine(
                start=phase_fraction_axes.c2p(0, 0.99),
                end=phase_fraction_axes.c2p(time_it_takes_99_perc.get_value(), 0.99),
                color=RED
            )
        )
        fraction_line_vertical_99_perc = always_redraw(
            lambda: DashedLine(
                start=phase_fraction_axes.c2p(time_it_takes_99_perc.get_value(), 0.99),
                end=ttt_axes.c2p(time_it_takes_99_perc.get_value(), temp_proxy_tracker.get_value() + 1),
                color=RED
            )
        )
        fraction_dot_99_perc = always_redraw(
            lambda: Dot(point=phase_fraction_axes.c2p(time_it_takes_99_perc.get_value(), 0.99),
                        color=RED)
        )
        ttt_dot_99_perc = always_redraw(
            lambda: Dot(point=ttt_axes.c2p(time_it_takes_99_perc.get_value(), temp_proxy_tracker.get_value() + 1),
                        color=RED)
        )
        ttt_dot_99_perc_path = TracedPath(ttt_dot_99_perc.get_center, stroke_color=RED, stroke_width=4)

        self.add(time_it_takes_99_perc)
        self.play(Create(fraction_line_horizontal_99_perc))
        self.play(Create(fraction_dot_99_perc))
        self.play(Create(fraction_line_vertical_99_perc))
        self.play(Create(ttt_dot_99_perc))
        self.add(ttt_dot_99_perc_path)

        ttt_temperature_line = always_redraw(
            lambda: DashedLine(
                start=ttt_axes.c2p(0, temp_proxy_tracker.get_value() + 1),
                end=ttt_axes.c2p(ttt_axes.x_range[1], temp_proxy_tracker.get_value() + 1),
                stroke_opacity=0.6
            )
        )
        self.play(Create(ttt_temperature_line))
        self.wait()

        self.play(
            temp_proxy_tracker.animate(run_time=animation_duration, rate_func=smooth).set_value(animation_duration))
        self.wait()


class TwoPhaseTimeCone(Scene):
    def construct(self):
        slide_title = Tex(
            "What happens when we have two phases?",
            font_size=40
        )
        slide_title.to_corner(UL)

        page_9 = Tex("$9$", font_size=40)
        page_9.to_edge(DR)

        axes = Axes(
            x_range=[-3, 3, 0.5],
            y_range=[0, 3, 0.5],
            x_length=10,
            y_length=6
        )

        axes.shift(DOWN * 0.8)

        x_label = axes.get_x_axis_label(Tex("$x$"))
        y_label = axes.get_y_axis_label(Tex("$t$"))

        point_of_interest_x_tracker = ValueTracker(axes.c2p(1.5, 2.5)[0])
        point_of_interest_y_tracker = ValueTracker(axes.c2p(1.5, 2.5)[1])

        green_phase_growth_rate_tracker = ValueTracker(0.2)
        orange_phase_growth_rate_tracker = ValueTracker(1 / 1.1)

        point_of_interest = always_redraw(
            lambda: Dot(point=np.array([point_of_interest_x_tracker.get_value(),
                                        point_of_interest_y_tracker.get_value(),
                                        0]))
        )

        inverted_time_cone_green_poi_left = always_redraw(
            lambda: DashedVMobject(
                axes.get_graph(
                    lambda x: (1 / green_phase_growth_rate_tracker.get_value()) * (
                            x - axes.p2c(point_of_interest.get_center())[0]) + axes.p2c(point_of_interest.get_center())[
                                  1],
                    x_range=[
                        axes.p2c(point_of_interest.get_center())[0] - axes.p2c(point_of_interest.get_center())[1] / (
                                1 / green_phase_growth_rate_tracker.get_value()),
                        axes.p2c(point_of_interest.get_center())[0]],
                    color=GREEN
                )
            )
        )

        inverted_time_cone_green_poi_right = always_redraw(
            lambda: DashedVMobject(
                axes.get_graph(
                    lambda x: -(1 / green_phase_growth_rate_tracker.get_value()) * (
                            x - axes.p2c(point_of_interest.get_center())[0]) + axes.p2c(point_of_interest.get_center())[
                                  1],
                    x_range=[
                        axes.p2c(point_of_interest.get_center())[0],
                        axes.p2c(point_of_interest.get_center())[0] + axes.p2c(point_of_interest.get_center())[1] / (
                                1 / green_phase_growth_rate_tracker.get_value())],
                    color=GREEN
                )
            )
        )

        inverted_time_cone_green_poi = VGroup(inverted_time_cone_green_poi_left,
                                              inverted_time_cone_green_poi_right)

        inverted_time_cone_orange_poi_left = always_redraw(
            lambda: DashedVMobject(
                axes.get_graph(
                    lambda x: (1 / orange_phase_growth_rate_tracker.get_value()) * (
                            x - axes.p2c(point_of_interest.get_center())[0]) + axes.p2c(point_of_interest.get_center())[
                                  1],
                    x_range=[
                        axes.p2c(point_of_interest.get_center())[0] - axes.p2c(point_of_interest.get_center())[1] / (
                                1 / orange_phase_growth_rate_tracker.get_value()),
                        axes.p2c(point_of_interest.get_center())[0]],
                    color=ORANGE
                )
            )
        )

        inverted_time_cone_orange_poi_right = always_redraw(
            lambda: DashedVMobject(
                axes.get_graph(
                    lambda x: -(1 / orange_phase_growth_rate_tracker.get_value()) * (
                            x - axes.p2c(point_of_interest.get_center())[0]) + axes.p2c(point_of_interest.get_center())[
                                  1],
                    x_range=[
                        axes.p2c(point_of_interest.get_center())[0],
                        axes.p2c(point_of_interest.get_center())[0] + axes.p2c(point_of_interest.get_center())[1] / (
                                1 / orange_phase_growth_rate_tracker.get_value())],
                    color=ORANGE
                )
            )
        )

        inverted_time_cone_orange_poi = VGroup(inverted_time_cone_orange_poi_left,
                                               inverted_time_cone_orange_poi_right)

        grow_rate_relation_1 = Tex(
            "$\\dot{R}_O$", color=ORANGE
        )
        grow_rate_relation_1.shift(RIGHT * 3 + UP * 3)

        grow_rate_relation_2 = Tex(
            "$>$"
        )
        grow_rate_relation_2.next_to(grow_rate_relation_1, RIGHT)

        grow_rate_relation_3 = Tex(
            "$\\dot{R}_G$", color=GREEN
        )
        grow_rate_relation_3.next_to(grow_rate_relation_2, RIGHT)

        grow_rate_relation = VGroup(grow_rate_relation_1,
                                    grow_rate_relation_2,
                                    grow_rate_relation_3)

        self.add(slide_title, page_9)
        self.wait()

        self.play(Create(axes))
        self.play(Write(VGroup(x_label, y_label)))
        self.wait()

        self.play(FadeIn(point_of_interest), FadeIn(inverted_time_cone_green_poi))
        self.wait()
        self.play(FadeIn(inverted_time_cone_orange_poi))
        self.wait()
        self.play(
            point_of_interest_x_tracker.animate.set_value(axes.c2p(0, 2.5)[0]),
            run_time=2
        )
        self.wait()

        self.play(Write(grow_rate_relation))
        self.wait()
        self.play(FadeOut(grow_rate_relation))
        self.wait()

        orange_x_tracker = ValueTracker(1)
        orange_t_tracker = ValueTracker(0.3)
        orange_t_cone_opacity_tracker = ValueTracker(1)

        time_cone_orange_left = always_redraw(
            lambda:
            axes.get_graph(
                lambda x: -(1 / orange_phase_growth_rate_tracker.get_value()) * (
                        x - orange_x_tracker.get_value()) + orange_t_tracker.get_value(),
                color=ORANGE,
                x_range=[axes.x_range[0], orange_x_tracker.get_value()]
            ).set_stroke(opacity=orange_t_cone_opacity_tracker.get_value())
        )

        time_cone_orange_right = always_redraw(
            lambda:
            axes.get_graph(
                lambda x: (1 / orange_phase_growth_rate_tracker.get_value()) * (
                        x - orange_x_tracker.get_value()) + orange_t_tracker.get_value(),
                color=ORANGE,
                x_range=[orange_x_tracker.get_value(), axes.x_range[1]]
            ).set_stroke(opacity=orange_t_cone_opacity_tracker.get_value())
        )

        time_cone_orange = VGroup(time_cone_orange_left, time_cone_orange_right)

        # TODO: polygon vertices always in 3D
        orange_phase_allowed_region = always_redraw(
            lambda: Polygon(
                point_of_interest.get_center(),
                axes.c2p(axes.p2c(point_of_interest.get_center())[0] - axes.p2c(point_of_interest.get_center())[1] / (
                        1 / orange_phase_growth_rate_tracker.get_value()), 0, 0),
                axes.c2p(axes.p2c(point_of_interest.get_center())[0] + axes.p2c(point_of_interest.get_center())[1] / (
                        1 / orange_phase_growth_rate_tracker.get_value()), 0, 0)
            ).set_stroke(color=ORANGE, opacity=0.4).set_fill(color=ORANGE, opacity=0.3)
        )

        self.play(FadeIn(time_cone_orange))
        self.wait()

        self.play(DrawBorderThenFill(orange_phase_allowed_region))
        self.wait()
        self.play(orange_x_tracker.animate(run_time=3).set_value(0.5),
                  orange_t_tracker.animate(run_time=3).set_value(0.8))
        self.wait()
        self.play(orange_x_tracker.animate(run_time=3).set_value(-1.5),
                  orange_t_tracker.animate(run_time=3).set_value(0.2))
        self.wait()
        self.play(FadeOut(orange_phase_allowed_region),
                  orange_t_cone_opacity_tracker.animate.set_value(0))
        self.wait()

        green_x_tracker = ValueTracker(0)
        green_t_tracker = ValueTracker(0.1)

        time_cone_green_left = always_redraw(
            lambda:
            axes.get_graph(
                lambda x: -(1 / green_phase_growth_rate_tracker.get_value()) * (
                        x - green_x_tracker.get_value()) + green_t_tracker.get_value(),
                color=GREEN,
                x_range=[axes.x_range[0], green_x_tracker.get_value()]
            )
        )

        time_cone_green_right = always_redraw(
            lambda:
            axes.get_graph(
                lambda x: (1 / green_phase_growth_rate_tracker.get_value()) * (
                        x - green_x_tracker.get_value()) + green_t_tracker.get_value(),
                color=GREEN,
                x_range=[green_x_tracker.get_value(), axes.x_range[1]]
            )
        )

        time_cone_green = VGroup(time_cone_green_left, time_cone_green_right)

        green_phase_allowed_region = always_redraw(
            lambda: Polygon(
                point_of_interest.get_center(),
                axes.c2p(axes.p2c(point_of_interest.get_center())[0] - axes.p2c(point_of_interest.get_center())[1] / (
                        1 / green_phase_growth_rate_tracker.get_value()), 0, 0),
                axes.c2p(axes.p2c(point_of_interest.get_center())[0] + axes.p2c(point_of_interest.get_center())[1] / (
                        1 / green_phase_growth_rate_tracker.get_value()), 0, 0)
            ).set_stroke(color=GREEN, opacity=0.4).set_fill(color=GREEN, opacity=0.3)
        )

        self.play(FadeIn(time_cone_green))
        self.wait()

        self.play(DrawBorderThenFill(green_phase_allowed_region))
        self.wait()
        self.play(green_t_tracker.animate(run_time=3).set_value(0.8))
        self.wait()
        self.play(green_x_tracker.animate(run_time=3).set_value(0.4),
                  green_t_tracker.animate(run_time=3).set_value(0.3))
        self.wait()
        self.play(green_x_tracker.animate(run_time=3).set_value(-1),
                  green_t_tracker.animate(run_time=3).set_value(0.5))
        self.wait()
        self.play(FadeOut(green_phase_allowed_region))
        self.wait()

        self.play(
            orange_x_tracker.animate(run_time=2).set_value(1.5),
            orange_t_tracker.animate(run_time=2).set_value(0.1),
            orange_t_cone_opacity_tracker.animate.set_value(1)
        )
        self.play(
            green_x_tracker.animate(run_time=2).set_value(-0.26),
            green_t_tracker.animate(run_time=2).set_value(0.4)
        )

        # impingement_x = always_redraw(
        #     lambda: (
        #                     orange_phase_growth_rate_tracker.get_value() * green_phase_growth_rate_tracker.get_value()
        #             ) / (
        #                     orange_phase_growth_rate_tracker.get_value() + green_phase_growth_rate_tracker.get_value()
        #             ) * (
        #                     (orange_x_tracker.get_value() / orange_phase_growth_rate_tracker.get_value()) + (
        #                     green_x_tracker.get_value() / green_phase_growth_rate_tracker.get_value()) +
        #                     orange_t_tracker.get_value() - green_t_tracker.get_value()
        #             ))
        #
        # impingement_t = always_redraw(
        #     lambda: (
        #                     green_phase_growth_rate_tracker.get_value() * green_t_tracker.get_value() +
        #                     orange_phase_growth_rate_tracker.get_value() * orange_t_tracker.get_value() -
        #                     green_x_tracker.get_value() + orange_x_tracker.get_value()
        #             ) / (orange_phase_growth_rate_tracker.get_value() + green_phase_growth_rate_tracker.get_value()))

        impingement_line = always_redraw(
            lambda: DashedLine(
                start=axes.c2p(((
                                        orange_phase_growth_rate_tracker.get_value() * green_phase_growth_rate_tracker.get_value()
                                ) / (
                                        orange_phase_growth_rate_tracker.get_value() + green_phase_growth_rate_tracker.get_value()
                                ) * (
                                        (
                                                orange_x_tracker.get_value() / orange_phase_growth_rate_tracker.get_value()) + (
                                                green_x_tracker.get_value() / green_phase_growth_rate_tracker.get_value()) +
                                        orange_t_tracker.get_value() - green_t_tracker.get_value()
                                )),
                               ((
                                        green_phase_growth_rate_tracker.get_value() * green_t_tracker.get_value() +
                                        orange_phase_growth_rate_tracker.get_value() * orange_t_tracker.get_value() -
                                        green_x_tracker.get_value() + orange_x_tracker.get_value()
                                ) / (
                                        orange_phase_growth_rate_tracker.get_value() + green_phase_growth_rate_tracker.get_value())),
                               0),
                end=axes.c2p(((
                                      orange_phase_growth_rate_tracker.get_value() * green_phase_growth_rate_tracker.get_value()
                              ) / (
                                      orange_phase_growth_rate_tracker.get_value() + green_phase_growth_rate_tracker.get_value()
                              ) * (
                                      (
                                              orange_x_tracker.get_value() / orange_phase_growth_rate_tracker.get_value()) + (
                                              green_x_tracker.get_value() / green_phase_growth_rate_tracker.get_value()) +
                                      orange_t_tracker.get_value() - green_t_tracker.get_value()
                              )),
                             5, 0)
            ).set_stroke(color=(ORANGE if ((
                                                   orange_phase_growth_rate_tracker.get_value() * green_phase_growth_rate_tracker.get_value()
                                           ) / (
                                                   orange_phase_growth_rate_tracker.get_value() + green_phase_growth_rate_tracker.get_value()
                                           ) * (
                                                   (
                                                           orange_x_tracker.get_value() / orange_phase_growth_rate_tracker.get_value()) + (
                                                           green_x_tracker.get_value() / green_phase_growth_rate_tracker.get_value()) +
                                                   orange_t_tracker.get_value() - green_t_tracker.get_value()
                                           )) <= 0 else GREEN))
        )

        # self.play(FadeIn(
        #     impingement_line
        # ))
        # self.wait()
        self.add(impingement_line)
        self.play(
            orange_x_tracker.animate(run_time=2).set_value(1.0)
        )
        self.wait()
        self.play(
            orange_x_tracker.animate(run_time=2).set_value(1.2)
        )
        self.wait()
        self.play(
            orange_t_tracker.animate(run_time=2).set_value(1.0)
        )
        self.wait()

        orange_win_region_1 = always_redraw(
            lambda: Polygon(
                axes.c2p(0,
                         -green_x_tracker.get_value() / green_phase_growth_rate_tracker.get_value() + green_t_tracker.get_value()),
                axes.c2p(green_x_tracker.get_value(), green_t_tracker.get_value()),
                axes.c2p(
                    green_x_tracker.get_value() + green_t_tracker.get_value() * orange_phase_growth_rate_tracker.get_value(),
                    0),
                axes.c2p(orange_phase_growth_rate_tracker.get_value() * (
                        green_t_tracker.get_value() - green_x_tracker.get_value() / green_phase_growth_rate_tracker.get_value()),
                         0)
            ).set_stroke(color=BLUE, opacity=0.6).set_fill(color=BLUE, opacity=0.3)
        )

        self.play(FadeIn(orange_win_region_1))
        self.wait()
        self.play(
            orange_x_tracker.animate.set_value(orange_phase_growth_rate_tracker.get_value() * (
                    green_t_tracker.get_value() - green_x_tracker.get_value() / green_phase_growth_rate_tracker.get_value()) - 0.1),
            orange_t_tracker.animate.set_value(0), run_time=2
        )
        self.play(
            orange_x_tracker.animate.set_value(0),
            orange_t_tracker.animate.set_value(
                -green_x_tracker.get_value() / green_phase_growth_rate_tracker.get_value() + green_t_tracker.get_value() - 0.1),
            run_time=2
        )
        self.wait()
        self.play(
            orange_x_tracker.animate.set_value(0),
            orange_t_tracker.animate.set_value(0.5),
            run_time=2
        )
        self.wait()
        self.play(
            orange_x_tracker.animate.set_value(
                green_x_tracker.get_value() + green_t_tracker.get_value() * orange_phase_growth_rate_tracker.get_value()),
            orange_t_tracker.animate.set_value(0),
            run_time=2
        )
        self.wait()
        self.play(
            orange_x_tracker.animate.set_value(2),
            orange_t_tracker.animate.set_value(0.2),
            run_time=2
        )
        self.wait()

        self.play(FadeOut(impingement_line))
        self.wait()
        impingement_line.suspend_updating()
        self.play(
            green_t_tracker.animate(run_time=2).set_value(0.1)
        )
        self.wait()
        self.play(
            green_x_tracker.animate(run_time=2).set_value(-0.45),
            green_t_tracker.animate(run_time=2).set_value(0)
        )
        self.wait()
        self.play(
            green_x_tracker.animate(run_time=2).set_value(0)
        )
        self.wait()
        impingement_line.resume_updating()
        self.add(impingement_line)
        self.play(
            orange_x_tracker.animate(run_time=2).set_value(0.2)
        )
        self.wait()
        self.play(FadeOut(impingement_line))
        self.wait()
        impingement_line.suspend_updating()

        self.play(
            green_x_tracker.animate(run_time=2).set_value(-0.15),
            green_t_tracker.animate(run_time=2).set_value(0.7)
        )
        self.wait()

        orange_win_region_2 = always_redraw(
            lambda: Polygon(
                axes.c2p(green_x_tracker.get_value(), green_t_tracker.get_value()),
                axes.c2p(
                    green_x_tracker.get_value() - green_t_tracker.get_value() * orange_phase_growth_rate_tracker.get_value(),
                    0),
                axes.c2p(
                    green_x_tracker.get_value() + green_t_tracker.get_value() * orange_phase_growth_rate_tracker.get_value(),
                    0)
            ).set_stroke(color=BLUE, opacity=0.6).set_fill(color=BLUE, opacity=0.3)
        )
        self.play(DrawBorderThenFill(orange_win_region_2))
        self.wait()
        self.play(
            orange_x_tracker.animate(run_time=2).set_value(-0.4)
        )
        self.wait()
        self.play(
            orange_t_tracker.animate(run_time=2).set_value(0)
        )
        self.wait()
        self.play(FadeOut(time_cone_orange))
        self.wait()
        time_cone_orange.suspend_updating()

        self.play(
            green_t_tracker.animate(run_time=2).set_value(1.35)
        )
        self.wait()
        self.play(
            green_x_tracker.animate(run_time=2).set_value(-0.05),
            green_t_tracker.animate(run_time=2).set_value(0.05)
        )
        self.wait()

        question_statement = Tex(
            "Question: ", "what is the fraction of the\\\\total nucleation volume where\\\\the orange phase can win?",
            font_size=30
        )
        question_statement.to_corner(UR)

        self.play(Write(question_statement))
        self.wait()
        self.play(FadeIn(orange_phase_allowed_region))
        self.wait()

        orange_phase_allowed_region_area_eqn_1 = Tex(
            "$||\\Omega_N||_{\\text{1D}}$", " $=$", " $\\dot{R}$", "$t^2$"
        )
        orange_phase_allowed_region_area_eqn_1.shift(LEFT * 3)

        orange_phase_allowed_region_area_eqn_2 = Tex(
            "$||\\Omega_N||_{\\text{1D}}$", " $=$", " $\\dot{R}_O$", "$t^2$"
        )
        orange_phase_allowed_region_area_eqn_2.align_to(orange_phase_allowed_region_area_eqn_1,
                                                        LEFT)

        self.play(Write(orange_phase_allowed_region_area_eqn_1))
        self.wait()
        self.play(TransformMatchingTex(orange_phase_allowed_region_area_eqn_1,
                                       orange_phase_allowed_region_area_eqn_2))
        self.wait()

        fraction = Tex(
            "$\\iint$", "$($",
            "$||\\Omega||_{\\text{orange wins}}$", "$/$", "$||\\Omega_N||_{\\text{tot}}$",
            "$)$", "$dx_Gdt_G$",
            font_size=40
        )
        fraction.to_edge(RIGHT, buff=0.1)
        fraction.shift(UP)

        fraction_orange_wins_eqn = Tex(
            "$||\\Omega||_O$", " $=$ ",
            "$\\frac{1}{6}$", "$t^2$", "$($", "$5$", "$\\dot{R}_O$", " $-$ ", "2", "$\\dot{R}_G$", "$)$"
        )
        fraction_orange_wins_eqn.next_to(fraction, DOWN)

        self.play(
            ReplacementTransform(VGroup(orange_win_region_1, orange_win_region_2).copy(),
                                 VGroup(fraction.submobjects[2]))
        )
        self.wait()
        self.play(
            ReplacementTransform(VGroup(orange_phase_allowed_region).copy(),
                                 VGroup(fraction.submobjects[4]))
        )
        self.wait()
        self.play(Write(fraction.submobjects[1]),
                  Write(fraction.submobjects[3]),
                  Write(fraction.submobjects[5]))
        self.wait()
        self.play(Write(fraction.submobjects[0]),
                  Write(fraction.submobjects[-1]))
        self.wait()
        self.play(ShowCreationThenFadeOut(green_phase_allowed_region))
        self.wait()

        self.play(Write(fraction_orange_wins_eqn))
        self.wait()
        self.play(FadeOut(fraction))
        self.wait()
        self.play(fraction_orange_wins_eqn.animate.move_to(question_statement.get_bottom(), UP).shift(DOWN * 0.1))
        self.wait()

        fraction_green_wins_eqn = Tex(
            "$||\\Omega||_G$", " $=$ ",
            "$\\frac{1}{6}$", "$t^2$", "$($", "$\\dot{R}_O$", " $+$ ", "2", "$\\dot{R}_G$", "$)$"
        )
        fraction_green_wins_eqn.next_to(fraction_orange_wins_eqn, DOWN).shift(DOWN * 0.1)
        self.play(Write(fraction_green_wins_eqn))
        self.wait()

        fraction_combined = Tex(
            "$||\\Omega||_O$", " $+$ ", "$||\\Omega||_G$",
            " $=$ ", "$\\dot{R}_O$", "$t^2$"
        )
        fraction_combined.next_to(fraction_green_wins_eqn, DOWN).shift(DOWN * 0.1)
        self.play(TransformMatchingShapes(VGroup(fraction_orange_wins_eqn, fraction_green_wins_eqn).copy(),
                                          VGroup(fraction_combined)))
        self.wait()


class TestSlide2(Scene):
    def construct(self):
        dot = Dot(radius=0.06)
        line = Line()
        line.set_stroke(width=dot.radius * DOT_TO_LINE)
        self.play(FadeIn(dot), run_time=2)
        self.wait(1)
        self.play(GrowFromPoint(line, point=dot), run_time=2, rate_func=linear)
        self.remove(dot)
        self.wait(1)
        self.play(line.animate.shift(UP * 2), run_time=2)
        self.wait(0.5)


class TestAxes(Scene):
    def construct(self):
        ax = Axes(x_range=[0, 15], y_range=[0, 10], x_length=8, y_length=6)
        # ax.move_to(DL)

        x_ax = ax.get_x_axis()
        x_ax_label = ax.get_x_axis_label(Tex("$x$"))
        y_ax = ax.get_y_axis()
        y_ax_label = ax.get_y_axis_label(Tex("$t$"))

        self.play(Create(x_ax), run_time=2)
        self.play(Write(x_ax_label), run_time=2)
        self.wait()
        self.play(Create(y_ax), run_time=2)
        self.play(Write(y_ax_label), run_time=2)
        self.wait(1)


class TestTitle(Scene):
    def construct(self):
        title = Text("Using the Time Cone Method to\nModel Competitive Phase Transformation", font_size=50)

        self.add(title)
        # TODO: Text().set_color_by_t2c() is bugged when using substrings as keys?
        self.play(title.animate.set_color_by_t2c(t2c={"[8:22]": BLUE, "[29:-1]": RED}))
        self.wait()


class TestLineToCircle(Scene):
    def construct(self):
        circle = Circle()
        line = Rectangle(height=0.1, width=circle.radius * 2)
        line.set_fill(color=WHITE, opacity=1)

        # self.add(line)
        # self.play(ClockwiseTransform(line, circle))
        self.play(ReplacementTransform(line, circle))
        self.wait()


class ThreeDCameraIllusionRotation(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        circle = Circle()
        self.set_camera_orientation(phi=75 * DEGREES, theta=30 * DEGREES)
        self.add(circle, axes)
        self.begin_3dillusion_camera_rotation(rate=2)
        self.wait(PI / 2)
        self.stop_3dillusion_camera_rotation()
