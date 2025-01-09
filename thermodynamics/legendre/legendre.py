import numpy as np
import sympy as smp

from manim import *

# set manim configs
config.flush_cache = True
config.disable_caching = True

# Define some constants
a = 0.01
b = 0.005
c = 150
offset = 50

# Define thermodynamic variables and potentials
s, v = smp.symbols("s v", real=True)
# symbolic U(S,V)
u = a * s ** 2 + b * (v - c) ** 2 + offset
# use u_f internal as a function that takes entropy and v as arguments
u_f = smp.lambdify(args=[s, v], expr=u)
# print(u_f(0.1, 0.1))
# get the partial derivatives and their physical variables
duds = smp.diff(u, s)
duds_f = smp.lambdify(args=[s], expr=duds)
t = duds
dudv = smp.diff(u, v)
dudv_f = smp.lambdify(args=[v], expr=dudv)
p = - dudv
# use legendre transform to obtain G(T,P)
g = u - t * s + p * v


# cd .\thermodynamics\legendre\
# manim -pqm .\legendre.py Legendre2D
class Legendre2D(Scene):
    def construct(self):
        # define the U(S) axes
        smin, smax = 0, 15
        umin, umax = 0, 25
        us_axes = Axes(
            x_range=[smin, smax],
            y_range=[umin, umax],
            x_length=6,
            y_length=5,
            axis_config={"include_ticks": False}
        ).shift(LEFT * 3.6)

        u_label = us_axes.get_y_axis_label(MathTex("U"))
        s_label = us_axes.get_x_axis_label(MathTex("S")).shift(DOWN * 0.5)

        self.add(us_axes)
        self.wait(0.5)
        self.play(
            Write(u_label),
            Write(s_label)
        )
        self.wait()

        # show the dU formula
        du_formula = MathTex(
            "dU","=","T","dS", "-", "P","dV",
            arg_separator=""
        ).shift(RIGHT * 2 + UP * 1.5)
        t_formula = MathTex(
            r"\left(\frac{\partial U}{\partial S}\right)_V",
            "=", "T", ">", "0",
            arg_separator=""
        ).next_to(du_formula, DOWN).align_to(du_formula, LEFT)
        self.play(
            Write(du_formula)
        )
        self.wait()
        self.play(
            Write(t_formula[:2])
        )
        self.wait(0.2)
        self.play(
            TransformFromCopy(
                du_formula[2], t_formula[2]
            )
        )
        self.wait(0.2)
        self.play(
            Write(t_formula[3:])
        )
        self.wait()

        # draw the linear U(s)
        m = 0.04
        b = 10
        linear_us = us_axes.plot(
            lambda x:
            m * smax ** 2 / (smax - smin) * x + b
        )
        self.play(Create(linear_us))
        self.wait()
        quadratic_us = us_axes.plot(
            lambda x:
            m * x ** 2 + b
        )
        # show the second derivative
        sec_deriv_formula = MathTex(
            r"\left(\frac{\partial T}{\partial S}\right)_V",
            "&=", r"T\left(\frac{\partial T}{\delta Q}\right)_V\\",
            "&=", r"\frac{T}{c_V}", ">", "0",
            arg_separator=""
        ).next_to(t_formula, DOWN).align_to(t_formula, LEFT)
        # hide the T in sec_deriv_formula
        sec_deriv_formula[0][2].set_opacity(0)
        t_in_sec = sec_deriv_formula[0][2].copy().set_opacity(1)
        self.play(
            TransformFromCopy(
                t_formula[2], t_in_sec
            )
        )
        self.play(
            Write(sec_deriv_formula[:-4])
        )
        self.wait(0.5)
        self.play(
            Write(sec_deriv_formula[-4:])
        )
        self.wait(0.5)
        self.play(
            ReplacementTransform(
                linear_us, quadratic_us
            )
        )
        self.wait()

        # show the overall transformation
        overall_formula = MathTex(
            r"U(S,V)\leftrightarrow F(T,V)"
        ).shift(3.5 * UP)
        self.play(
            Write(overall_formula)
        )
        self.wait()
        # highlight U->F
        rightarrow = MathTex(r"\rightarrow", color=BLUE).move_to(overall_formula[0][6])
        self.play(
            Create(rightarrow),
            overall_formula[0][6].animate.set_opacity(0.3)
        )
        self.wait(0.5)
        self.play(
            FadeOut(t_formula, sec_deriv_formula, t_in_sec)
        )
        self.wait()
        # show that F is F(T,V)
        f_formula_old = MathTex(
            "F", "=", "U", "-", r"\left(\frac{\partial U}{\partial S}\right)_V", "S"
        ).next_to(du_formula, DOWN).align_to(du_formula, LEFT)
        f_formula_new = MathTex(
            "F", "=", "U", "-", "T", "S"
        ).next_to(du_formula, DOWN).align_to(du_formula, LEFT)
        self.play(Write(f_formula_old))
        self.wait()
        self.play(TransformMatchingTex(f_formula_old, f_formula_new))
        self.wait()
        df_formula_old = MathTex(
            "dF", "=", "dU", "-", "T", "dS", "-", "S", "dT"
        ).next_to(f_formula_new, DOWN).align_to(f_formula_new, LEFT)
        df_formula_new = MathTex(
            "dF", "=", "TdS", "-", "PdV", "-", "T", "dS", "-", "S", "dT",
        ).next_to(f_formula_new, DOWN).align_to(f_formula_new, LEFT)
        self.play(Write(df_formula_old))
        self.wait(0.5)
        self.play(TransformMatchingTex(df_formula_old, df_formula_new))
        self.wait()
        # cancel out terms
        cancel_line_1 = Line(
            start=df_formula_new[2].get_boundary_point(LEFT + DOWN),
            end=df_formula_new[2].get_boundary_point(RIGHT + UP),
            stroke_width=DEFAULT_STROKE_WIDTH * 0.7
        )
        cancel_line_2 = Line(
            start=df_formula_new[6:8].get_boundary_point(LEFT + DOWN),
            end=df_formula_new[6:8].get_boundary_point(RIGHT + UP),
            stroke_width=DEFAULT_STROKE_WIDTH * 0.7
        )
        self.play(
            Create(cancel_line_1),
            Create(cancel_line_2)
        )
        df_formula_final = MathTex(
            "dF", "=", "-", "S", "dT", "-", "PdV",
        ).next_to(f_formula_new, DOWN).align_to(f_formula_new, LEFT)
        self.play(
            FadeOut(cancel_line_1, cancel_line_2, run_time=0.5),
            TransformMatchingTex(df_formula_new, df_formula_final)
        )
        self.wait()

        # demonstrate the graphical intuition
        s_tracker = ValueTracker((smax - smin) / 2)
        su_vert = always_redraw(
            lambda:
            DashedLine(
                start=us_axes.c2p(s_tracker.get_value(), umin),
                end=us_axes.c2p(s_tracker.get_value(),
                                m * s_tracker.get_value() ** 2 + b),
            )
        )
        su_dot = Dot(
            point=us_axes.c2p(
                s_tracker.get_value(),
                m * s_tracker.get_value() ** 2 + b)
        ).add_updater(
            lambda mob:
            mob.move_to(
                us_axes.c2p(
                    s_tracker.get_value(),
                    m * s_tracker.get_value() ** 2 + b
                )
            )
        )

        su_hori = always_redraw(
            lambda:
            DashedLine(
                start=us_axes.c2p(s_tracker.get_value(),
                                  m * s_tracker.get_value() ** 2 + b),
                end=us_axes.c2p(smin,
                                m * s_tracker.get_value() ** 2 + b),
            )
        )
        self.play(
            LaggedStart(
                Create(su_vert),
                DrawBorderThenFill(su_dot),
                Create(su_hori),
                lag_ratio=1.1
            ),
            run_time=1.2
        )
        self.wait()

        u_label = MathTex("U", font_size=DEFAULT_FONT_SIZE * 0.7)
        minus_ts_label = MathTex("-TS", font_size=DEFAULT_FONT_SIZE * 0.7)
        minus_ts_label.add_updater(
            lambda mob:
            mob.next_to(
                us_axes.c2p(
                    smin, (
                    m * s_tracker.get_value() ** 2 + b
                    - 2 * m * s_tracker.get_value() ** 2
                    + m * s_tracker.get_value() ** 2 + b
                    ) / 2
                ),
                LEFT
            )
        )
        f_label = MathTex("F", font_size=DEFAULT_FONT_SIZE * 0.7)
        f_label.add_updater(
            lambda mob:
            mob.next_to(
                us_axes.c2p(
                    smin, (
                    - 2 * m * s_tracker.get_value() ** 2
                    + m * s_tracker.get_value() ** 2 + b
                    )
                ),
                LEFT
            )
        )
        # add the slope line
        slope = always_redraw(
            lambda:
            us_axes.plot(
                lambda x:
                2 * m * s_tracker.get_value() * (x-s_tracker.get_value())
                + m * s_tracker.get_value() ** 2 + b,
                x_range=[smax, smin, (smin - smax) / 100]
            ).set_opacity(0.5)
        )
        self.play(Create(slope))
        self.wait()
        # show the connection to f formula
        minus_ts_arrow = Arrow(
            start=us_axes.c2p(smin, m * s_tracker.get_value() ** 2 + b),
            end=us_axes.c2p(
                smin, (
                        - 2 * m * s_tracker.get_value() ** 2
                        + m * s_tracker.get_value() ** 2 + b
                )
            ),
            buff=0
        )
        minus_ts_arrow.sort(lambda mob: np.dot(mob, DOWN))
        self.play(DrawBorderThenFill(minus_ts_arrow))
        self.wait(0.2)
        f_dot = Dot(
            point=us_axes.c2p(
                smin,
                (
                        - 2 * m * s_tracker.get_value() ** 2
                        + m * s_tracker.get_value() ** 2 + b
                )
            ),
            z_index=100
        )
        f_dot.add_updater(
            lambda mob:
            mob.move_to(
                us_axes.c2p(
                    smin,
                    (
                            - 2 * m * s_tracker.get_value() ** 2
                            + m * s_tracker.get_value() ** 2 + b
                    )
                )
            )
        )
        self.play(DrawBorderThenFill(f_dot))
        self.wait()
        # indicate the graphical equivalent
        su_dot_final = su_dot.copy().set_color(GREEN)
        u_in_f = f_formula_new[2]
        u_in_f_final = u_in_f.copy().set_color(GREEN)
        self.play(
            Succession(
                Transform(su_dot, su_dot.copy().scale(1.4).set_color(GREEN), run_time=0.5),
                Wait(run_time=0.5),
                Transform(su_dot, su_dot_final, run_time=0.5)
            ),
            Succession(
                Transform(u_in_f, u_in_f.copy().scale(1.4).set_color(GREEN), run_time=0.5),
                Wait(run_time=0.5),
                Transform(u_in_f, u_in_f_final, run_time=0.5)
            ),
            rate_func=smooth
        )
        self.wait()
        minus_ts_in_f = f_formula_new[3:]
        minus_ts_in_f_final = minus_ts_in_f.copy().set_color(BLUE)
        minus_ts_arrow_final = minus_ts_arrow.copy().set_color(BLUE)
        self.play(
            Succession(
                Transform(minus_ts_arrow, minus_ts_arrow.copy().scale(1.4).set_color(BLUE), run_time=0.5),
                Wait(run_time=0.5),
                Transform(minus_ts_arrow, minus_ts_arrow_final, run_time=0.5)
            ),
            Succession(
                Transform(minus_ts_in_f, minus_ts_in_f.copy().scale(1.4).set_color(BLUE), run_time=0.5),
                Wait(run_time=0.5),
                Transform(minus_ts_in_f, minus_ts_in_f_final, run_time=0.5)
            ),
            rate_func=smooth
        )
        minus_ts_arrow.add_updater(
            lambda mob:
            mob.become(
                Arrow(
                    start=us_axes.c2p(smin, m * s_tracker.get_value() ** 2 + b),
                    end=us_axes.c2p(
                        smin, (
                                - 2 * m * s_tracker.get_value() ** 2
                                + m * s_tracker.get_value() ** 2 + b
                        )
                    ),
                    buff=0, color=BLUE
                )
            )
        )
        self.wait()
        f_in_f = f_formula_new[0]
        f_in_f_final = f_in_f.copy().set_color(ORANGE)
        f_dot_final = f_dot.copy().set_color(ORANGE)
        self.play(
            Succession(
                Transform(f_dot, f_dot.copy().scale(1.4).set_color(ORANGE), run_time=0.5),
                Wait(run_time=0.5),
                Transform(f_dot, f_dot_final, run_time=0.5)
            ),
            Succession(
                Transform(f_in_f, f_in_f.copy().scale(1.4).set_color(ORANGE), run_time=0.5),
                Wait(run_time=0.5),
                Transform(f_in_f, f_in_f_final, run_time=0.5)
            ),
            rate_func=smooth
        )
        self.wait()
        self.play(
            s_tracker.animate.set_value(smin)
        )
        self.wait()
        # fade away all the equations
        self.play(
            FadeOut(
                du_formula,
                f_formula_new,
                df_formula_final
            )
        )
        self.wait()

        # show the F(T) axes
        tmin, tmax = 2 * m * smin, 2 * m * smax
        fmin, fmax = umin, umax
        ft_axes = Axes(
            x_range=[tmin, tmax],
            y_range=[fmin, fmax],
            x_length=6,
            y_length=5,
            axis_config={"include_ticks": False}
        ).shift(RIGHT * 3.6)

        f_label = ft_axes.get_y_axis_label(MathTex("F"))
        t_label = ft_axes.get_x_axis_label(MathTex("T")).shift(DOWN * 0.5 + RIGHT * 0.2)

        self.play(FadeIn(ft_axes))
        self.play(
            Write(f_label),
            Write(t_label)
        )
        self.wait()
        # set a temperature tracker
        t_tracker = ValueTracker(2 * m * s_tracker.get_value()).add_updater(
            lambda tracker, dt: tracker.set_value(2 * m * s_tracker.get_value())
        )
        f_point = always_redraw(
            lambda:
            Dot(
                point=ft_axes.c2p(
                    t_tracker.get_value(),
                    m * s_tracker.get_value() ** 2 + b - t_tracker.get_value() * s_tracker.get_value()
                ),
                color=ORANGE
            )
        )
        # add horizontal lines connecting the f_intercept the f axis
        f_hori_connect = always_redraw(
            lambda:
            DashedLine(
                start=f_dot, end=f_point
            )
        )
        # trace out the helmholtz path
        f_plot = TracedPath(f_point.get_center, stroke_color=ORANGE, stroke_width=4)
        self.add(t_tracker)
        self.play(Create(f_hori_connect))
        self.play(FadeIn(f_point))
        self.add(f_plot)
        self.wait()
        self.play(
            s_tracker.animate.set_value(smax),
            run_time=2
        )
        f_plot.clear_updaters()
        self.wait()

        # show the equivalent of reverse transformation
        minus_ts_arrow.clear_updaters()
        self.play(
            FadeOut(minus_ts_arrow, su_hori, su_vert, su_dot),
            quadratic_us.animate.set_stroke(opacity=0.3)
        )
        self.wait(0.5)
        leftarrow = MathTex(r"\leftarrow", color=RED).move_to(overall_formula[0][6])
        self.play(
            ReplacementTransform(rightarrow, leftarrow),
            t_label.animate.set_color(RED)
        )
        self.wait()

        slope_red = always_redraw(
            lambda:
            us_axes.plot(
                lambda x:
                2 * m * s_tracker.get_value() * (x - s_tracker.get_value())
                + m * s_tracker.get_value() ** 2 + b
            ).set_color(RED).set_opacity(0.5)
        )
        # slope.sort(lambda mob: np.dot(mob, LEFT))
        self.play(
            Create(slope_red),
            Uncreate(slope)
        )
        self.wait()
        
        
        # add the outer envelopes
        def draw_u_envelopes():
            envelopes = [VMobject()]
            s_interval = (smax - smin) / 20
            s_upper = smax
            s_lower = smax - (smax - s_tracker.get_value()) // s_interval * s_interval
            s_sampled = np.arange(s_upper, s_lower, -s_interval)
            slopes = [
                us_axes.plot(
                    lambda x:
                    2 * m * s_val * (x - s_val)
                    + m * s_val ** 2 + b
                ).set_color(RED).set_opacity(0.5)
                for s_val in s_sampled
            ]
            envelopes.extend(slopes)
            return VGroup(
                *envelopes
            )
        
        u_envelopes = always_redraw(draw_u_envelopes)
        self.add(u_envelopes)
        
        self.play(
            s_tracker.animate.set_value(smin),
            run_time=3
        )
        u_envelopes.clear_updaters()
        self.wait()
        self.play(
            Create(quadratic_us.copy().set_stroke(opacity=1, color=GREEN)),
            run_time=2
        )
        self.wait()


# cd .\thermodynamics\legendre\
# manim -pqm -n 22 .\legendre.py Legendre3D
class Legendre3D(ThreeDScene):
    def construct(self):
        # Define the 3D axes
        xmin, xmax = 0, 100
        ymin, ymax = 0, 150
        zmin, zmax = 0, 300
        axes = ThreeDAxes(
            x_range=[xmin, xmax, (xmax - xmin) / 10],
            y_range=[ymin, ymax, (ymax - ymin) / 10],
            z_range=[zmin, zmax, (zmax - zmin) / 10],
            x_length=6,
            y_length=6,
            z_length=5
        )

        # get the individual axis
        x_axis = axes.get_x_axis()
        y_axis = axes.get_y_axis()
        z_axis = axes.get_z_axis()

        axes.x_axis = x_axis.rotate(PI / 2, axis=X_AXIS)
        unit_y_length = axes.c2p(0, 1, 0)[1] - axes.c2p(0, 0, 0)[1]
        unit_z_length = axes.c2p(0, 0, 1)[2] - axes.c2p(0, 0, 0)[2]
        unit_z_vector = axes.c2p(0, 0, 1) - axes.c2p(0, 0, 0)

        # add the axes labels
        x_label = axes.get_x_axis_label(Tex("$S$")).rotate(PI / 2, axis=X_AXIS)
        y_label = axes.get_y_axis_label(Tex("$V$")).rotate(PI / 2, axis=X_AXIS)
        z_label = axes.get_z_axis_label(Tex("$U$"))

        # show the quasi-2D view with U(S)
        self.set_camera_orientation(
            90 * DEGREES, -90 * DEGREES,
            frame_center=axes.get_center()
        )
        self.play(
            FadeIn(x_axis, z_axis)
        )
        self.play(
            Write(x_label),
            Write(z_label)
        )
        self.wait(0.1)

        # Create the U(S)|_V
        v_tracker = ValueTracker(0)
        u_s_curve = always_redraw(
            lambda:
            axes.plot_parametric_curve(
                lambda x:
                np.array([
                    x,
                    v_tracker.get_value(),
                    u_f(x, v_tracker.get_value())
                ]),
                color=RED, t_range=[xmin, xmax]
            )
        )
        self.play(
            Create(u_s_curve)
        )
        self.wait()

        # move camera angle to show the V axis
        self.move_camera(
            phi=70 * DEGREES
        )
        # display the hidden V axis
        self.play(
            Create(y_axis)
        )
        self.play(
            Write(y_label)
        )
        self.wait(0.5)
        # show different U(S)|_V cross-sections
        # first create a plane to indicate the value of V
        v_cross_sec = always_redraw(
            lambda:
            Polygon(
                axes.c2p(xmin, ymin, zmax),
                axes.c2p(xmax, ymin, zmax),
                axes.c2p(xmax, ymin, zmin),
                axes.c2p(xmin, ymin, zmin),
                color=WHITE, fill_opacity=0.2
            ).shift(
                v_tracker.get_value() * (axes.c2p(0, 1, 0) - axes.c2p(0, 0, 0))
            )
        )
        self.play(
            DrawBorderThenFill(v_cross_sec)
        )
        self.wait()
        self.play(
            v_tracker.animate.set_value(ymax),
            run_time=1.5
        )
        self.wait()

        # switch to U(V)|_S perspective
        self.move_camera(
            phi=90 * DEGREES, theta=0,
            added_anims=[
                # hide the entropy x-axis and the U(S)|_V cross-section
                FadeOut(u_s_curve), FadeOut(v_cross_sec),
                x_axis.animate.set_opacity(0),
                x_label.animate.set_opacity(0),
                # rotate the U and V axis and labels
                z_axis.animate.rotate(PI / 2, axis=Z_AXIS),
                z_label.animate.rotate(PI / 2, axis=Z_AXIS,
                               about_point=axes.c2p(xmin, ymin, zmax)),
                y_axis.animate.rotate(PI / 2, axis=Y_AXIS),
                y_label.animate.rotate(PI / 2, axis=Z_AXIS,
                               about_point=axes.c2p(xmin, ymax, zmin)),
            ]
        )

        # Create the U(V)|_S
        s_tracker = ValueTracker(0)
        u_v_curve = always_redraw(
            lambda:
            axes.plot_parametric_curve(
                lambda x:
                np.array([
                    s_tracker.get_value(),
                    x,
                    u_f(s_tracker.get_value(), x)
                ]),
                color=BLUE, t_range=[ymin, ymax]
            )
        )
        self.play(
            Create(u_v_curve)
        )
        self.wait()
        # move camera angle to show the S axis
        u_v_curve.suspend_updating()
        self.move_camera(
            phi=70 * DEGREES,
            added_anims=[
                # rotate the S axis
                x_axis.animate.set_opacity(1).rotate(-PI / 2, axis=X_AXIS),
                x_label.animate.set_opacity(1).rotate(-PI / 2, axis=X_AXIS)
            ]
        )
        u_v_curve.resume_updating()
        # show different U(V)|_S cross-sections
        # first create a plane to indicate the value of V
        s_cross_sec = always_redraw(
            lambda:
            Polygon(
                axes.c2p(xmin, ymin, zmax),
                axes.c2p(xmin, ymax, zmax),
                axes.c2p(xmin, ymax, zmin),
                axes.c2p(xmin, ymin, zmin),
                color=WHITE, fill_opacity=0.2
            ).shift(
                s_tracker.get_value() * (axes.c2p(1, 0, 0) - axes.c2p(0, 0, 0))
            )
        )
        self.play(
            DrawBorderThenFill(s_cross_sec)
        )
        self.play(
            s_tracker.animate.set_value(xmax),
            run_time=1.5
        )
        self.wait()

        # move back to the original 3D view
        self.move_camera(
            theta=-90 * DEGREES,
            added_anims=[
                # show the v cross-section again
                FadeIn(u_s_curve), FadeIn(v_cross_sec),
                # rotate the U and V axis and labels
                z_axis.animate.rotate(-PI / 2, axis=Z_AXIS),
                z_label.animate.rotate(-PI / 2, axis=Z_AXIS,
                                       about_point=axes.c2p(xmin, ymin, zmax)),
                y_axis.animate.rotate(-PI / 2, axis=Y_AXIS),
                y_label.animate.rotate(-PI / 2, axis=Z_AXIS,
                                       about_point=axes.c2p(xmin, ymax, zmin)),
                # rotate the S axis
                x_axis.animate.set_opacity(1).rotate(PI / 2, axis=X_AXIS),
                x_label.animate.set_opacity(1).rotate(PI / 2, axis=X_AXIS)
            ]
        )
        self.wait(0.5)

        # Create the U surface
        u_surf = Surface(
            lambda x, y: axes.c2p(
                x, y, u_f(x, y)
            ),
            u_range=[xmin, xmax],
            v_range=[ymin, ymax],
            fill_opacity=0.3,
            fill_color=WHITE,
            checkerboard_colors=[WHITE]
        ).set_shade_in_3d().sort(lambda num: np.dot(num, DOWN + LEFT))
        self.play(
            Create(u_surf),
            run_time=2
        )
        self.wait()
        # move to the middle point of the S, V plane
        self.play(
            FadeOut(s_cross_sec),
            FadeOut(v_cross_sec)
        )
        self.play(
            s_tracker.animate.set_value((xmax - xmin) / 2),
            v_tracker.animate.set_value((ymax - ymin) / 2)
        )
        self.wait()

        # Show the tangent plane
        # z_axis = Z_AXIS
        # normal_vect = np.array([-duds_f(s_tracker.get_value()), -dudv_f(v_tracker.get_value()), 1])
        # normal_vect /= np.linalg.norm(normal_vect)
        #
        # rotate_axis = np.cross(unit_z_vector, axes.c2p(*normal_vect)-axes.c2p(0,0,0))
        # rotate_axis /= np.linalg.norm(rotate_axis)
        # rotate_angle = angle_between_vectors(z_axis, normal_vect)
        # print(normal_vect)
        # print(rotate_angle / DEGREES)

        tangent_width = (xmax - xmin) / 5

        def get_tangent_z(x0, y0, x, y):
            return u_f(x0, y0) + duds_f(x0) * (x - x0) + dudv_f(y0) * (y - y0)

        def draw_tangent():
            x0, y0 = s_tracker.get_value(), v_tracker.get_value()
            # x1, y1 = x0 + tangent_width, y0 + tangent_width
            # x2, y2 = x0 + tangent_width, y0 - tangent_width
            # x3, y3 = x0 - tangent_width, y0 - tangent_width
            # x4, y4 = x0 - tangent_width, y0 + tangent_width
            x1, y1 = xmin, ymax
            x2, y2 = xmax, ymax
            x3, y3 = xmax, ymin
            x4, y4 = xmin, ymin
            plane = Polygon(
                axes.c2p(x1, y1, get_tangent_z(x0, y0, x1, y1)),
                axes.c2p(x2, y2, get_tangent_z(x0, y0, x2, y2)),
                axes.c2p(x3, y3, get_tangent_z(x0, y0, x3, y3)),
                axes.c2p(x4, y4, get_tangent_z(x0, y0, x4, y4)),
                color=YELLOW, stroke_opacity=0.8,
                fill_color=WHITE, fill_opacity=0.3,
                stroke_width=DEFAULT_STROKE_WIDTH * 0.5
            )
            return plane

        tangent_plane = always_redraw(
            draw_tangent
        )
        self.play(
            DrawBorderThenFill(tangent_plane)
        )
        self.wait()

        self.move_camera(
            theta=-110 * DEGREES,
            phi=80 * DEGREES
        )
        self.wait()
        self.move_camera(
            theta=-120 * DEGREES,
            phi=60 * DEGREES
        )
        self.wait()

        # show the F, H tangent lines
        u_point = always_redraw(
            lambda:
            Sphere(
                center=axes.c2p(s_tracker.get_value(), v_tracker.get_value(),
                                u_f(s_tracker.get_value(), v_tracker.get_value())),
                radius=0.1
            ).set_color(PURPLE_C)
        )
        f_point = always_redraw(
            lambda:
            Sphere(
                center=axes.c2p(0, v_tracker.get_value(),
                                get_tangent_z(s_tracker.get_value(), v_tracker.get_value(), 0, v_tracker.get_value())),
                radius=0.1
            ).set_color(RED)
        )
        h_point = always_redraw(
            lambda:
            Sphere(
                center=axes.c2p(s_tracker.get_value(), 0,
                                get_tangent_z(s_tracker.get_value(), v_tracker.get_value(), s_tracker.get_value(), 0)),
                radius=0.1
            ).set_color(BLUE)
        )
        f_line = always_redraw(
            lambda:
            DashedLine(
                start=u_point.get_center(),
                end=f_point.get_center(),
                color=RED,
                z_index=-10
            )
        )
        self.play(
            LaggedStart(
                FadeIn(u_point),
                Create(f_line),
                FadeIn(f_point),
                lag_ratio=1
            )
        )
        self.wait()
        h_line = always_redraw(
            lambda:
            DashedLine(
                start=u_point.get_center(),
                end=h_point.get_center(),
                color=BLUE,
                z_index=-10
            )
        )
        self.play(
            LaggedStart(
                Create(h_line),
                FadeIn(h_point),
                lag_ratio=1
            )
        )
        self.wait()

        # Show the connection to Gibbs
        # the following two are horizontal lines
        u_to_f = always_redraw(
            lambda:
            Line(
                start=u_point.get_center(),
                end=axes.c2p(0, v_tracker.get_value(),
                             u_f(s_tracker.get_value(), v_tracker.get_value())),
                color=WHITE
            )
        )
        u_to_h = always_redraw(
            lambda:
            Line(
                start=u_point.get_center(),
                end=axes.c2p(s_tracker.get_value(), 0,
                             u_f(s_tracker.get_value(), v_tracker.get_value())),
                color=WHITE
            )
        )
        minus_ts = Arrow3D(
            start=axes.c2p(0, v_tracker.get_value(),
                           u_f(s_tracker.get_value(), v_tracker.get_value())),
            end=f_point.get_center(),
            color=RED
        )
        plus_pv = Arrow3D(
            start=axes.c2p(s_tracker.get_value(), 0,
                           u_f(s_tracker.get_value(), v_tracker.get_value())),
            end=h_point.get_center(),
            color=BLUE
        )
        self.play(
            Create(u_to_f)
        )
        self.play(
            Create(minus_ts)
        )
        self.wait()
        self.play(
            Create(u_to_h)
        )
        self.play(
            Create(plus_pv)
        )
        self.wait()
        # show the u value on the z axis
        f_to_z_axis = always_redraw(
            lambda:
            Line(
                start=axes.c2p(0, v_tracker.get_value(),
                               u_f(s_tracker.get_value(), v_tracker.get_value())),
                end=axes.c2p(0, 0,
                             u_f(s_tracker.get_value(), v_tracker.get_value()))
            )
        )
        h_to_z_axis = always_redraw(
            lambda:
            Line(
                start=axes.c2p(s_tracker.get_value(), 0,
                               u_f(s_tracker.get_value(), v_tracker.get_value())),
                end=axes.c2p(0, 0, u_f(s_tracker.get_value(), v_tracker.get_value()))
            )
        )
        self.play(
            Create(f_to_z_axis),
            Create(h_to_z_axis)
        )

        # show to subtraction/addition process
        minus_ts_on_axis = minus_ts.copy()
        plus_pv_on_axis = plus_pv.copy()
        self.play(
            minus_ts_on_axis.animate.shift(axes.c2p(0, 0, 0) - axes.c2p(0, v_tracker.get_value(), 0))
        )
        self.play(
            plus_pv_on_axis.animate.shift(axes.c2p(0, 0, 0) - axes.c2p(s_tracker.get_value(), 0, 0))
        )
        self.play(
            plus_pv_on_axis.animate.shift(
                f_point.get_center() - axes.c2p(0, v_tracker.get_value(),
                                                u_f(s_tracker.get_value(), v_tracker.get_value())))
        )
        self.wait(0.5)
        self.move_camera(
            phi=80 * DEGREES
        )
        self.wait(0.5)

        g_point = always_redraw(
            lambda:
            Sphere(
                center=axes.c2p(0, 0,
                                get_tangent_z(s_tracker.get_value(), v_tracker.get_value(), 0, 0)),
                radius=0.1
            ).set_color(GREEN)
        )
        u_to_g = always_redraw(
            lambda:
            DashedLine(
                start=u_point.get_center(),
                end=g_point.get_center(),
                color=GREEN
            )
        )
        self.play(
            LaggedStart(
                FadeIn(g_point),
                Create(u_to_g),
                lag_ratio=1
            )
        )
        self.wait()
        self.move_camera(
            theta=-110 * DEGREES
        )
        self.wait()
        self.move_camera(
            theta=-120 * DEGREES, phi=75 * DEGREES
        )
        self.wait()

        # add the thermodynamic potential label
        F_label = MathTex(
            "F(T,V)", font_size=20, color=RED
        ).next_to(f_point, LEFT).rotate(PI/2, axis=X_AXIS)
        H_label = MathTex(
            "H(S,P)", font_size=20, color=BLUE
        ).next_to(h_point, -Y_AXIS).rotate(PI/2, axis=X_AXIS).rotate(-PI/2, axis=Z_AXIS)
        G_label = MathTex(
            "G(T,P)", font_size=20, color=GREEN
        ).next_to(g_point, LEFT).rotate(PI/2, axis=X_AXIS)
        self.play(
            LaggedStart(
                Write(F_label),
                Write(H_label),
                Write(G_label),
                lag_ratio=1
            )
        )
        self.wait(0.5)

        self.play(
            FadeOut(minus_ts), FadeOut(minus_ts_on_axis),
            FadeOut(plus_pv), FadeOut(plus_pv_on_axis)
        )
        self.wait(0.5)

        F_label.add_updater(lambda lab: lab.next_to(f_point, LEFT))
        H_label.add_updater(lambda lab: lab.next_to(h_point, -Y_AXIS))
        G_label.add_updater(lambda lab: lab.next_to(g_point, LEFT))

        # vary S and V
        self.play(
            s_tracker.animate.set_value(20),
            v_tracker.animate.set_value(30)
        )
        self.wait()
        self.play(
            s_tracker.animate.set_value(80),
            v_tracker.animate.set_value(130)
        )
        self.wait()
