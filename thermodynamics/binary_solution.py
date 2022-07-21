from abc import ABC

import numpy as np
from manim import *
from manim.opengl import *
from manim.mobject.opengl.opengl_compatibility import ConvertToOpenGL
from typing import List
from scipy.spatial import Delaunay

config.disable_caching = True
config.flush_cache = True


class BinarySolution(Scene):
    @staticmethod
    def get_axis_ticks_labels(values: List[tuple], axis: Axes, fontsize=20, label_dir=UP) -> VGroup:
        axis_labels = VGroup()
        for tick_val, tick_str in values:
            tex = Tex(tick_str, font_size=fontsize)
            tex.next_to(axis.n2p(tick_val), label_dir)
            axis_labels.add(tex)
        return axis_labels

    @staticmethod
    def add_axis_ticks(values, axis, inplace=True):
        """Add axis ticks inplace"""
        ticks = VGroup()
        for _ in values:
            ticks.add(axis.get_tick(_, axis.tick_size))
        if inplace:
            axis.add(ticks)
            axis.ticks = ticks
        return ticks

    @staticmethod
    def add_nl_ticks(values, nl, inplace=True):
        """Add NumberLine ticks inplace"""
        ticks = VGroup()
        for _ in values:
            ticks.add(nl.get_tick(_, nl.tick_size))
        if inplace:
            nl.add(ticks)
            nl.ticks = ticks
        return ticks

    def construct(self):
        x_min, x_max, x_step = 0, 1, 0.1
        y_min, y_max, y_step = -5000, 3000, 1000
        axes = Axes(
            x_range=[x_min, x_max, x_step],
            y_range=[y_min, y_max, y_step],
            x_length=6,
            y_length=6.5,
            axis_config={"include_ticks": False},
            tips=False
        ).shift(UP * 0.3 + LEFT * 2.5)

        # get the x, y axes separately
        x_axis = axes.get_x_axis().copy().shift(
            axes.c2p(axes.x_range[0], axes.y_range[0]) - axes.c2p(axes.x_range[0], 0)
        )
        y_axis = axes.get_y_axis()

        # add the ticks and tick labels
        values_x = [
            (tick_val, str(tick_val))
            for tick_val in np.arange(x_min * 10, (x_max + x_step) * 10, x_step * 10) / 10
        ]
        values_y = [
            (tick_val, str(int(tick_val / 1e3)))
            for tick_val in np.arange(y_min, y_max + y_step, y_step)
        ]

        x_tick_labels = self.get_axis_ticks_labels(values_x, x_axis, fontsize=25, label_dir=DOWN)
        y_tick_labels = self.get_axis_ticks_labels(values_y, y_axis, fontsize=25, label_dir=LEFT)
        x_ticks = self.add_axis_ticks([val[0] for val in values_x], x_axis, inplace=False)
        self.add_axis_ticks([val[0] for val in values_y], axes.y_axis)

        # get the up and right bounding line
        up_bound_line = Line(
            start=axes.c2p(x_min, y_max),
            end=axes.c2p(x_max, y_max),
            stroke_width=2
        )
        right_bound_line = Line(
            start=axes.c2p(x_max, y_min),
            end=axes.c2p(x_max, y_max),
            stroke_width=2
        )

        self.add(x_axis, y_axis, up_bound_line, right_bound_line, x_ticks)
        self.wait()

        # add the axis titles
        x_label = axes.get_x_axis_label(
            MathTex(r"X_\text{B}", font_size=35)
        ).next_to(x_tick_labels, DOWN)
        y_label = axes.get_y_axis_label(
            Tex("$G$ (kJ/mol)", font_size=35).rotate(PI / 2)
        ).next_to(y_tick_labels, LEFT)

        self.play(
            Write(x_label),
            Write(y_label)
        )
        self.wait()

        # add a dashed hline for the G=0
        zero_hline = DashedLine(
            start=axes.c2p(x_min, 0),
            end=axes.c2p(x_max, 0),
            stroke_width=2
        )

        # add the tick labels
        self.play(
            Write(x_tick_labels),
            Write(y_tick_labels),
            Create(zero_hline)
        )
        self.wait()

        # define relevant variables
        a_ref_tracker = ValueTracker(-2900)
        b_ref_tracker = ValueTracker(-1000)
        temp_tracker = ValueTracker(470)
        omega_tracker = ValueTracker(9500)

        # define a series of gibbs functions
        def get_ref_gibbs(x):
            gibbs = a_ref_tracker.get_value() * (1 - x) + b_ref_tracker.get_value() * x
            return gibbs

        def get_ideal_gibbs(x):
            if x in {0, 1}:
                return 0
            gibbs = 8.31 * temp_tracker.get_value() * (x * np.log(x) + (1 - x) * np.log(1 - x))
            return gibbs

        def get_excess_gibbs(x):
            gibbs = (1 - x) * x * omega_tracker.get_value()
            return gibbs

        gibbs_funcs = [get_ref_gibbs, get_ideal_gibbs, get_excess_gibbs]

        # define the formulas
        formula_fontsize = 30
        tot_gibbs_formula = MathTex(
            "G_m", "=", "G_m^0", "+", r"\Delta G_m^{\text{ideal}}", "+", r"\Delta G_m^{\text{xs}}",
            font_size=formula_fontsize
        ).set_opacity(0).next_to(axes, RIGHT, buff=2 * DEFAULT_MOBJECT_TO_MOBJECT_BUFFER).shift(UP * 1.5)
        tot_gibbs_formula[:3].set_opacity(1)

        ref_gibbs_formula = MathTex(
            "G_m^0", "=", r"x_\text{A}G_\text{A}^0+x_\text{B}G_\text{B}^0",
            font_size=formula_fontsize
        ).next_to(tot_gibbs_formula, DOWN).align_to(tot_gibbs_formula, LEFT)
        ref_gibbs_formula[2][2:5].set_color(GREEN)
        ref_gibbs_formula[2][-3:].set_color(RED)

        ideal_gibbs_formula = MathTex(
            r"\Delta G_m^{\text{ideal}}", "=", r"RT\left(x_\text{A}\ln{x_\text{A}}+x_\text{B}\ln{x_\text{B}}\right)",
            font_size=formula_fontsize
        ).next_to(ref_gibbs_formula, DOWN).align_to(ref_gibbs_formula, LEFT)

        excess_gibbs_formula = MathTex(
            r"\Delta G_m^{\text{xs}}", "=", r"x_\text{A}x_\text{B}L_{\text{A,B}}",
            font_size=formula_fontsize
        ).next_to(ideal_gibbs_formula, DOWN).align_to(ideal_gibbs_formula, LEFT)

        formulas = [ref_gibbs_formula, ideal_gibbs_formula, excess_gibbs_formula]

        self.play(
            Write(tot_gibbs_formula),
            Write(ref_gibbs_formula)
        )
        self.wait()

        # show the reference states
        a_ref_dot = Dot(
            point=axes.c2p(0, a_ref_tracker.get_value()),
            color=GREEN, z_index=10
        )
        b_ref_dot = Dot(
            point=axes.c2p(1, b_ref_tracker.get_value()),
            color=RED, z_index=10
        )
        self.play(
            DrawBorderThenFill(a_ref_dot),
            DrawBorderThenFill(b_ref_dot)
        )
        self.wait()

        gibbs_old = axes.plot(
            get_ref_gibbs
        )
        self.play(
            Create(gibbs_old),
            ref_gibbs_formula.animate.set_color(WHITE).set_opacity(0.6),
            a_ref_dot.animate.set_color(WHITE),
            b_ref_dot.animate.set_color(WHITE)
        )
        self.wait()

        num_vlines = 50

        for _ in range(1, len(gibbs_funcs)):
            # display the new formula to add
            formula_to_add = formulas[_]
            formula_to_add.set_color(YELLOW)

            self.play(
                Write(formula_to_add)
            )
            self.wait()

            # define the new gibbs term to add
            gibbs_to_add = axes.plot(
                gibbs_funcs[_],
                color=YELLOW, stroke_opacity=0.6
            )
            # get the vlines to show the addition/subtraction
            vlines = axes.get_vertical_lines_to_graph(
                gibbs_to_add, num_lines=num_vlines
            )
            for vline in vlines:
                vline.set_color_by_gradient(WHITE, YELLOW)
            # define the shift vectors for each vline
            shift_vecs = [
                axes.c2p(
                    x_coord, sum([func(x_coord) for func in gibbs_funcs[:_]])
                ) - axes.c2p(
                    x_coord, 0
                ) for x_coord in np.linspace(0, 1, num_vlines)
            ]
            # define the corresponding shift animations
            shift_anims = LaggedStart(
                *[vline.animate.shift(shift_vector)
                  for vline, shift_vector in zip(vlines, shift_vecs)],
                lag_ratio=0.3
            )

            # define the aggregated gibbs energy
            gibbs_new = axes.plot(
                lambda x: sum([func(x) for func in gibbs_funcs[:_ + 1]])
            )

            self.play(
                Create(gibbs_to_add)
            )
            self.wait()

            # add the formula_to_add to the total gibbs formula
            symbol_copy = formula_to_add[0].copy()
            self.play(
                # display the + sign
                Write(tot_gibbs_formula[_ * 2 + 1].copy().set_opacity(1)),
                # show the added symbol in the total gibbs formula
                symbol_copy.animate.move_to(tot_gibbs_formula[(_ + 1) * 2])
            )
            self.wait()

            self.play(
                Create(vlines)
            )
            self.wait()
            self.play(
                shift_anims,
                run_time=1.5
            )
            self.wait()

            self.play(
                Transform(gibbs_old, gibbs_new),
                formula_to_add.animate.set_color(WHITE).set_opacity(0.6),
                symbol_copy.animate.set_color(WHITE)
            )
            self.wait()

            self.play(
                FadeOut(gibbs_to_add, vlines)
            )
            self.wait()

        # change the variables
        gibbs_new_dynamic = always_redraw(
            lambda:
            axes.plot(
                lambda x: sum([func(x) for func in gibbs_funcs])
            )
        )
        self.play(
            ReplacementTransform(gibbs_old, gibbs_new_dynamic),
            run_time=0.1
        )

        def get_slide(center_val, offset) -> NumberLine:
            min_val, max_val = center_val - offset, center_val + offset
            print([min_val, max_val, offset])
            nl = NumberLine(
                x_range=[min_val, max_val, offset],
                length=2,
                include_ticks=False
            ).rotate(PI / 2)
            # manually add the ticks
            self.add_nl_ticks(
                values=np.arange(min_val, max_val + offset, offset),
                nl=nl
            )
            return nl

        def get_slider_content(nl: NumberLine, tracker: ValueTracker, div_num=1, **kwargs) -> VGroup:
            slider_dot = Dot(
                point=nl.n2p(tracker.get_value()),
                color=YELLOW
            )
            slider_val = DecimalNumber(
                number=tracker.get_value() / div_num,
                **kwargs
            ).next_to(slider_dot, RIGHT,
                      buff=0.5 * DEFAULT_MOBJECT_TO_MOBJECT_BUFFER
                      )
            return VGroup(slider_dot, slider_val)

        def get_slider(*args, **kwargs):
            slider = always_redraw(
                lambda:
                get_slider_content(*args, **kwargs)
            )
            return slider

        def get_slider_unit_label(mob: Mobject, label_str: str, alpha, **kwargs):
            return Tex(
                f"({label_str})",
                **kwargs
            ).set_opacity(alpha).next_to(
                mob, RIGHT, buff=0.3 * DEFAULT_MOBJECT_TO_MOBJECT_BUFFER
            )

        # show the variables
        sliders_buff = 4.5 * DEFAULT_MOBJECT_TO_MOBJECT_BUFFER
        unit_label_alpha = 0.6
        unit_label_fontsize = 20
        a_ref_label = ref_gibbs_formula[2][2:5].copy()
        self.play(
            a_ref_label.animate.next_to(
                excess_gibbs_formula, DOWN, buff=2 * DEFAULT_MOBJECT_TO_MOBJECT_BUFFER
            ).align_to(excess_gibbs_formula, LEFT).set_opacity(1)
        )
        a_ref_unit = get_slider_unit_label(a_ref_label, "kJ/mol", unit_label_alpha, font_size=unit_label_fontsize)
        self.play(
            Write(a_ref_unit)
        )
        a_ref_slide = get_slide(a_ref_tracker.get_value(), offset=2000)
        a_ref_slide.next_to(a_ref_label, DOWN)
        a_ref_slider = get_slider(a_ref_slide, a_ref_tracker,
                                  div_num=1000, num_decimal_places=1, font_size=20)
        self.play(
            LaggedStart(
                Create(a_ref_slide),
                DrawBorderThenFill(a_ref_slider),
                lag_ratio=0.5
            )
        )

        b_ref_label = ref_gibbs_formula[2][8:].copy()
        self.play(
            b_ref_label.animate.next_to(
                a_ref_label, RIGHT, buff=sliders_buff
            ).set_opacity(1)
        )
        b_ref_unit = get_slider_unit_label(b_ref_label, "kJ/mol", unit_label_alpha, font_size=unit_label_fontsize)
        self.play(
            Write(b_ref_unit)
        )
        b_ref_slide = get_slide(b_ref_tracker.get_value(), offset=2000)
        b_ref_slide.next_to(b_ref_label, DOWN).align_to(a_ref_slide, UP)
        b_ref_slider = get_slider(b_ref_slide, b_ref_tracker,
                                  div_num=1000, num_decimal_places=1, font_size=20)
        self.play(
            LaggedStart(
                Create(b_ref_slide),
                DrawBorderThenFill(b_ref_slider),
                lag_ratio=0.5
            )
        )

        temp_label = ideal_gibbs_formula[2][1].copy()
        self.play(
            temp_label.animate.next_to(
                b_ref_label, RIGHT, buff=sliders_buff
            ).set_opacity(1)
        )
        temp_unit = get_slider_unit_label(temp_label, "K", unit_label_alpha, font_size=unit_label_fontsize)
        self.play(
            Write(temp_unit)
        )
        temp_slide = get_slide(temp_tracker.get_value(), offset=300)
        temp_slide.next_to(temp_label, DOWN).align_to(a_ref_slide, UP)
        temp_slider = get_slider(temp_slide, temp_tracker,
                                 num_decimal_places=0, font_size=20)
        self.play(
            LaggedStart(
                Create(temp_slide),
                DrawBorderThenFill(temp_slider),
                lag_ratio=0.5
            )
        )

        omega_label = excess_gibbs_formula[2][-4:].copy()
        self.play(
            omega_label.animate.next_to(
                temp_label, RIGHT, buff=sliders_buff
            ).set_opacity(1)
        )
        omega_unit = get_slider_unit_label(omega_label, "kJ/mol", unit_label_alpha, font_size=unit_label_fontsize)
        self.play(
            Write(omega_unit)
        )
        omega_slide = get_slide(omega_tracker.get_value(), offset=10000)
        omega_slide.next_to(omega_label, DOWN).align_to(a_ref_slide, UP)
        omega_slider = get_slider(omega_slide, omega_tracker,
                                  div_num=1000, num_decimal_places=1, font_size=20)
        self.play(
            LaggedStart(
                Create(omega_slide),
                DrawBorderThenFill(omega_slider),
                lag_ratio=0.5
            )
        )
        self.wait()

        # add dot updaters
        a_ref_dot.add_updater(
            lambda mob: mob.move_to(axes.c2p(0, a_ref_tracker.get_value()))
        )
        b_ref_dot.add_updater(
            lambda mob: mob.move_to(axes.c2p(1, b_ref_tracker.get_value()))
        )

        # change the a_ref energy
        self.play(
            a_ref_tracker.animate.set_value(-4000)
        )
        self.wait()

        # change the b_ref energy
        self.play(
            b_ref_tracker.animate.set_value(-500)
        )
        self.wait()

        # change the temperature
        self.play(
            temp_tracker.animate.set_value(300)
        )
        self.wait()

        self.play(
            temp_tracker.animate.set_value(600)
        )
        self.wait()

        # change the interaction term
        self.play(
            omega_tracker.animate.set_value(12000)
        )
        self.wait()

        self.play(
            omega_tracker.animate.set_value(8000)
        )
        self.wait()


config.window_position = "0,100"


class Simplex:
    """Modified from https://github.com/materialsproject/pymatgen/blob/v2022.7.19/pymatgen/util/coord.py#L367
    Thank you, Pymatgen :)
    """

    def __init__(self, vertices: List[list]):
        self.vertices = np.array(vertices)
        self.dim = self.vertices.shape[0]
        self.augmented_matrix = np.concatenate([vertices, np.ones((self.dim, 1))], axis=-1)

    def bary_to_cart(self, bary_coords):
        return np.dot(bary_coords, self.augmented_matrix[:, :-1])


# simplex = Simplex(
#     [[0, 0],  # Component A
#      [1, 0],  # Component B
#      [0.5, np.sqrt(3)/2]]  # Component C
# )
# x, y, *_ = simplex.bary_to_cart(
#     [1/2, 1/2, 0]
# )
# print(np.sqrt(3)/2)
# print(x, y)


class TernarySurface(VGroup, metaclass=ConvertToOpenGL):
    def __init__(
            self,
            vertices,
            z_func: Callable,
            axes: ThreeDAxes,
            resolution=None,
            u_range=[0, 1],
            v_range=[0, 1],
            surface_piece_config: dict = {},
            fill_color: Color = BLUE_D,
            fill_opacity: float = 1.0,
            checkerboard_colors: Sequence[Color] = [BLUE_D, BLUE_E],
            stroke_color: Color = LIGHT_GREY,
            stroke_width: float = 0.5,
            should_make_jagged: bool = False,
            pre_function_handle_to_anchor_scale_factor: float = 0.00001,
            **kwargs,
    ):
        self.simplex = Simplex(vertices)
        self.z_func = z_func
        self.axes = axes
        self.u_range = u_range,
        self.v_range = v_range,
        super().__init__(**kwargs)
        self.surface_piece_config = surface_piece_config
        self.fill_color = fill_color
        self.fill_opacity = fill_opacity
        self.checkerboard_colors = checkerboard_colors
        self.stroke_color = stroke_color
        self.stroke_width = stroke_width
        self.should_make_jagged = should_make_jagged
        self.pre_function_handle_to_anchor_scale_factor = (
            pre_function_handle_to_anchor_scale_factor
        )

        res_value = (30, 30)

        self.resolution = resolution if resolution is not None else res_value
        self.delauney_tri = self.run_delauney()
        self.faces = None
        self.get_facets()
        # self.apply_function(lambda p: self.get_cart_coord(p[0], p[1]))
        # if self.should_make_jagged:
        #     self.make_jagged()

    def _get_u_values_and_v_values(self):
        res = tuplify(self.resolution)
        if len(res) == 1:
            u_res = v_res = res[0]
        else:
            u_res, v_res = res

        u_values = np.linspace(*self.u_range[0], u_res + 1)
        v_values = np.linspace(*self.v_range[0], v_res + 1)

        # create a mesh based on u, v values
        uu, vv = np.meshgrid(u_values, v_values)
        uv_sum_values = uu + vv
        # the sum should not be greater than 1
        mask = (uv_sum_values <= 1)
        u_values = uu[mask]
        v_values = vv[mask]
        return u_values, v_values

    def get_cart_coord(self, u, v):
        w = 1 - u - v
        x, y, *_ = self.simplex.bary_to_cart([u, v, w])
        z = self.z_func(u, v)
        return [x, y, z]

    def run_delauney(self):
        u_values, v_values = self._get_u_values_and_v_values()
        coords = [
            self.get_cart_coord(u, v)[:2]
            for u, v in zip(u_values, v_values)
        ]
        tri = Delaunay(
            points=np.array(coords)
        )
        return tri

    def get_facets(self):
        faces = VGroup()
        u_values, v_values = self._get_u_values_and_v_values()
        all_vertices = [
            self.get_cart_coord(u, v)
            for u, v in zip(u_values, v_values)
        ]
        simplices = self.delauney_tri.simplices
        for facet in simplices:
            face = ThreeDVMobject(shade_in_3d=False)
            vertices = [
                all_vertices[_]
                for _ in np.append(facet, [facet[0]])
            ]
            vertices = [
                self.axes.c2p(*coords)
                for coords in vertices
            ]
            face.set_points_as_corners(
                vertices
            )
            faces.add(face)
        # set the aesthetics
        faces.set_fill(color=self.fill_color, opacity=self.fill_opacity)
        faces.set_stroke(
            color=self.stroke_color,
            width=self.stroke_width,
            opacity=self.stroke_opacity,
        )
        self.faces = faces
        self.add(*faces)


class TernarySolution(ThreeDScene):
    @staticmethod
    def get_axis_ticks_labels(values: List[tuple], axis: Axes, fontsize=20, label_dir=UP) -> VGroup:
        axis_labels = VGroup()
        for tick_val, tick_str in values:
            tex = Tex(tick_str, font_size=fontsize)
            tex.next_to(axis.n2p(tick_val), label_dir)
            axis_labels.add(tex)
        return axis_labels

    @staticmethod
    def add_axis_ticks(values, axis: Axes, inplace=True):
        """Add axis ticks inplace"""
        # TODO: need to investigate add ticks in 3D
        ticks = VGroup()
        for _ in values:
            tick_line = axis.get_tick(_, axis.tick_size)
            tick_line.set_stroke(width=10)
            tick_line.set_angle(0)
            ticks.add(tick_line)
        if inplace:
            axis.add(ticks)
            axis.ticks = ticks
        return ticks

    def construct(self):
        # define the axes attributes
        x_min, x_max, x_step = 0, 1, 0.1
        y_min, y_max, y_step = 0, 1, 0.1
        z_min, z_max, z_step = 0, 8000, 1000
        z_offset = 5000
        axes_3d = ThreeDAxes(
            x_range=[x_min, x_max, x_step],
            y_range=[y_min, y_max, y_step],
            z_range=[z_min, z_max, z_step],
            x_length=6,
            y_length=6,
            z_length=6.5,
            tips=False
        )
        # set the ternary base
        default_vertices = [[0, 0, 0],
                            [1, 0, 0],
                            [1 / 2, np.sqrt(3) / 2, 0]]

        default_vertices_centroid = center_of_mass(default_vertices)
        default_vertices_centroid = axes_3d.c2p(*default_vertices_centroid)
        frame_center = [*default_vertices_centroid[:2], axes_3d.copy().get_center()[2]]

        self.set_camera_orientation(
            phi=90 * DEGREES, theta=-90 * DEGREES,
            frame_center=frame_center,
            focal_distance=100 * self.camera.get_focal_distance()
        )
        self.add(Dot3D(frame_center))  # for debug
        # shift vector for off-center view
        shift_vec_off_center = np.array([-2.5, 2, -3])
        # shift vector for center view
        shift_vec_center = np.array([0, 0, -1])
        axes_3d.shift(shift_vec_center)

        # get the individual axis
        x_axis = axes_3d.get_x_axis().rotate(PI / 2, axis=X_AXIS)
        z_axis = axes_3d.get_z_axis()

        # add the axis ticks and tick labels
        values_x = [
            (tick_val, str(tick_val))
            for tick_val in np.arange(x_min * 10, (x_max + x_step) * 10, x_step * 10) / 10
        ]
        values_z = [
            (tick_val, str(int((tick_val - z_offset) / 1000)))
            for tick_val in np.arange(z_min, z_max + z_step, z_step)
        ]

        x_tick_labels = self.get_axis_ticks_labels(values_x, x_axis, fontsize=25, label_dir=DOWN)
        for lab in x_tick_labels:
            lab.rotate(PI / 2, axis=X_AXIS, about_point=x_axis.get_center())
        z_tick_labels = self.get_axis_ticks_labels(values_z, z_axis, fontsize=25, label_dir=LEFT)
        for lab in z_tick_labels:
            lab.rotate(PI / 2, axis=X_AXIS)

        # define some aesthetic variable
        axis_title_fontsize = 35

        x_label = axes_3d.get_x_axis_label(
            MathTex(r"X_\text{B}", font_size=axis_title_fontsize)
        ).rotate(PI / 2, axis=X_AXIS).next_to(x_tick_labels, -Z_AXIS)
        z_label = Tex(
            "$G$ (kJ/mol)", font_size=axis_title_fontsize
        ).rotate(
            PI / 2, axis=X_AXIS
        ).rotate(
            -PI / 2, axis=Y_AXIS
        ).next_to(
            z_tick_labels, -X_AXIS
        )

        # add an extra tick at the left end
        old_ticks = x_axis.ticks
        new_ticks = VGroup()
        extra_tick = old_ticks[0].copy()
        extra_tick.move_to(x_axis.n2p(0))
        x_axis.add(extra_tick)
        new_ticks.add(extra_tick)
        for tick in old_ticks:
            new_ticks.add(tick)
        x_axis.ticks = new_ticks
        x_ticks = x_axis.ticks

        self.add(x_axis, z_axis, x_label, z_label, x_tick_labels, z_tick_labels)
        self.wait()

        def get_bisect_vec(point_a, point_b, point_c, ax):
            point_a = ax.c2p(*point_a)
            point_b = ax.c2p(*point_b)
            point_c = ax.c2p(*point_c)

            vec_1 = 0.1 * (point_a - point_b)
            vec_2 = 0.1 * (point_a - point_c)

            bisect_vec = (vec_1 + vec_2) / np.sqrt(2)
            return bisect_vec

        # move camera angle to show 3D view and the corner elements
        a_label = Tex(
            "A", font_size=axis_title_fontsize
        ).move_to(
            axes_3d.c2p(*default_vertices[0]) + get_bisect_vec(default_vertices[0], default_vertices[1],
                                                               default_vertices[2], axes_3d)
        )

        b_label = MathTex(
            r"\text{B}", font_size=axis_title_fontsize
        ).move_to(
            axes_3d.c2p(*default_vertices[1]) + get_bisect_vec(default_vertices[1], default_vertices[0],
                                                               default_vertices[2], axes_3d)
        )

        c_label = Tex(
            "C", font_size=axis_title_fontsize
        ).move_to(
            axes_3d.c2p(*default_vertices[2]) + get_bisect_vec(default_vertices[2], default_vertices[0],
                                                               default_vertices[1], axes_3d)
        )

        a_label_rotated = a_label.copy().rotate(PI / 2, axis=X_AXIS)
        b_label_rotated = b_label.copy().rotate(PI / 2, axis=X_AXIS)
        c_label_rotated = c_label.copy().rotate(PI / 2, axis=X_AXIS)
        self.move_camera(
            phi=80 * DEGREES, theta=-85 * DEGREES,
            added_anims=[
                ReplacementTransform(x_label[0][1].copy(), b_label_rotated),
                FadeOut(x_label)
            ]
        )
        # AHA: looks like need to add fixed orientation mobjects without prior rotation?
        self.add_fixed_orientation_mobjects(b_label)
        self.remove(b_label_rotated)

        # show the ternary base and move the axis
        ternary_base = ThreeDVMobject(
            fill_opacity=0.2, color=WHITE, stroke_width=1
        )
        ternary_base_vertices = [
            axes_3d.c2p(*vertex)
            # AHA: manim ThreeDVMobject need a closed loop for corner points
            for vertex in default_vertices + [default_vertices[0]]
        ]
        ternary_base.set_points_as_corners(
            ternary_base_vertices
        )

        z_axis_at_B = z_axis.copy()
        z_axis_at_C = z_axis.copy()

        self.play(
            LaggedStart(
                DrawBorderThenFill(ternary_base),
                LaggedStart(
                    AnimationGroup(
                        z_axis_at_B.animate.shift(ternary_base_vertices[1] - ternary_base_vertices[0]),
                        z_axis_at_C.animate.shift(ternary_base_vertices[2] - ternary_base_vertices[0]),
                    ),
                    AnimationGroup(
                        Write(a_label_rotated),
                        Write(c_label_rotated),
                    ),
                    lag_ratio=0.95
                ),
                lag_ratio=0.3
            )
        )
        self.add_fixed_orientation_mobjects(a_label, c_label)
        self.remove(a_label_rotated, c_label_rotated)
        self.wait()

        # move to top-down view
        self.move_camera(
            phi=0 * DEGREES, theta=-90 * DEGREES,
            # AHA: just the origin
            frame_center=Y_AXIS,
            added_anims=[
                # set the z axes opacity to 0
                z_axis.animate.set_opacity(0),
                z_axis_at_B.animate.set_opacity(0),
                z_axis_at_C.animate.set_opacity(0),
                z_tick_labels.animate.set_opacity(0),
                z_label.animate.set_opacity(0),
                # rotate the x-axis ticks and tick labels to x-y plane
                x_ticks.animate.rotate(-PI / 2, axis=X_AXIS),
                x_tick_labels.animate.rotate(-PI / 2, axis=X_AXIS, about_point=x_axis.get_center())
            ]
        )
        self.wait()

        # rotate each individual tick
        tick_rotation_anims = [
            tick.animate.rotate(-30 * DEGREES, axis=Z_AXIS, about_point=tick.get_center())
            for tick in x_ticks
        ]

        # shift each tick label to match the orientation of the tick
        tick_shift_vec = rotate_vector(
            vector=(x_axis.n2p(0.1) - x_tick_labels[1].get_top()) / np.sqrt(3),
            angle=PI / 2
        )
        print(tick_shift_vec)

        tick_label_shift_anims = [
            tick_label.animate.shift(tick_shift_vec)
            for tick_label in x_tick_labels
        ]

        tick_anim_combos = [
            AnimationGroup(
                anim_1, anim_2
            )
            for anim_1, anim_2 in zip(tick_rotation_anims, tick_label_shift_anims)
        ]

        self.play(
            LaggedStart(
                *tick_anim_combos,
                lag_ratio=0.5
            ),
            run_time=2
        )
        self.wait()

        # move back to the 3D view
        # restore the ticks and tick labels
        tick_label_restore_anims = AnimationGroup(
            *[
                tick_label.animate.shift(
                    -tick_shift_vec
                ).rotate(PI / 2, axis=X_AXIS, about_point=x_axis.get_center())
                for tick_label in x_tick_labels
            ]
        )

        tick_restore_anims = AnimationGroup(
            *[
                tick.animate.rotate(
                    30 * DEGREES, axis=Z_AXIS, about_point=tick.get_center()
                ).rotate(PI / 2, axis=X_AXIS)
                for tick in x_ticks
            ]
        )

        self.move_camera(
            phi=80 * DEGREES, theta=-90 * DEGREES,
            frame_center=frame_center,
            added_anims=[
                # set the z axes opacity to 1
                z_axis.animate.set_opacity(1),
                z_axis_at_B.animate.set_opacity(1),
                z_axis_at_C.animate.set_opacity(1),
                z_tick_labels.animate.set_opacity(1),
                z_label.animate.set_opacity(1),
                # rotate the x-axis ticks and tick labels back to x-z plane
                tick_label_restore_anims,
                tick_restore_anims,
            ]
        )
        self.wait()

        # define the thermodynamic variables
        a_ref = -3000
        b_ref = -1000
        c_ref = -2000
        temp = 400
        omega_ab = 9000
        omega_ac = 7000
        omega_bc = 4000
        omega_abc = 50000

        def get_gibbs(x_a, x_b):
            x_c = 1 - x_a - x_b
            gibbs = x_a * a_ref + x_b * b_ref + x_c * c_ref
            return gibbs + z_offset

        def get_entropy(x):
            if x <= 0:
                return 0
            return x * np.log(x)

        def get_gibbs_b(x_a, x_b):
            x_c = 1 - x_a - x_b

            gibbs_ref = x_a * a_ref + x_b * b_ref + x_c * c_ref
            gibbs_ideal = 8.31 * temp * (get_entropy(x_a) + get_entropy(x_b) + get_entropy(x_c))
            return gibbs_ref + gibbs_ideal + z_offset

        def get_gibbs_c(x_a, x_b):
            x_c = 1 - x_a - x_b

            gibbs_ref = x_a * a_ref + x_b * b_ref + x_c * c_ref
            gibbs_ideal = 8.31 * temp * (get_entropy(x_a) + get_entropy(x_b) + get_entropy(x_c))
            gibbs_excess = (
                    x_a * x_b * omega_ab +
                    x_a * x_c * omega_ac +
                    x_b * x_c * omega_bc +
                    x_a * x_b * x_c * omega_abc
            )
            return gibbs_ref + gibbs_ideal + gibbs_excess + z_offset

        gibbs_surface_1 = TernarySurface(vertices=default_vertices,
                                         z_func=get_gibbs,
                                         axes=axes_3d,
                                         fill_color=WHITE, fill_opacity=0.3,
                                         stroke_opacity=0.2)

        gibbs_surface_2 = TernarySurface(vertices=default_vertices,
                                         z_func=get_gibbs_b,
                                         axes=axes_3d,
                                         fill_color=PURPLE, fill_opacity=0.7,
                                         stroke_opacity=1)

        gibbs_surface_3 = TernarySurface(vertices=default_vertices,
                                         z_func=get_gibbs_c,
                                         axes=axes_3d,
                                         fill_color=BLUE, fill_opacity=0.7,
                                         stroke_opacity=1)

        self.play(FadeIn(gibbs_surface_1))
        self.wait()

        self.play(
            ReplacementTransform(gibbs_surface_1, gibbs_surface_2)
        )
        self.wait()

        self.play(
            ReplacementTransform(gibbs_surface_2, gibbs_surface_3)
        )
        self.wait()

        rotation_time = 5
        for _ in range(3):
            self.begin_ambient_camera_rotation(
                rate=120 / rotation_time * DEGREES
            )
            self.wait(rotation_time)
            self.stop_ambient_camera_rotation()
            self.wait()

        # rotation for off center view
        # axes_3d_copy = axes_3d.copy()
        # for _ in range(3):
        #     self.play(
        #         VGroup(
        #             axes_3d, gibbs_surface_2, x_label, y_label, z_label
        #         ).animate(
        #             run_time=5, rate_func=linear
        #         ).rotate(
        #             angle=-120 * DEGREES,
        #             axis=axes_3d_copy.c2p(0, 0, 1) - axes_3d_copy.c2p(0, 0, 0),
        #             about_point=axes_3d_copy.get_center()
        #         )
        #     )
        #     self.wait()

# debug
# TernarySolution().render()
