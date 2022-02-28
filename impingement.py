from manim import *


# noinspection DuplicatedCode
def get_circ_intersec(x1, y1, x2, y2, r1, r2):
    distance = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    delta = 0.25 * np.sqrt((distance + r1 + r2) * (distance + r1 - r2) * (distance - r1 + r2) * (-distance + r1 + r2))
    pm, mp = np.array([1, -1]), np.array([-1, 1])
    intersec_1_x, intersec_2_x = (x1 + x2) / 2 + (x2 - x1) * (r1 ** 2 - r2 ** 2) / (2 * distance ** 2) + 2 * pm * (
            y1 - y2) * delta / distance ** 2
    intersec_1_y, intersec_2_y = (y1 + y2) / 2 + (y2 - y1) * (r1 ** 2 - r2 ** 2) / (2 * distance ** 2) + 2 * mp * (
            x1 - x2) * delta / distance ** 2
    return VGroup(Dot(point=np.array([intersec_1_x, intersec_1_y, 0])),
                  Dot(point=np.array([intersec_2_x, intersec_2_y, 0])))


class TwoCircleImpingement(Scene):
    def construct(self):
        time_tracker = ValueTracker(0)

        circ_1_x, circ_1_y, circ_1_z = -1.5, 0, 0
        circ_1_start_radius = 1.5
        circ_1_growth_rate = 0.5
        circ_1_radius_tracker = ValueTracker(circ_1_start_radius).add_updater(
            lambda tracker: tracker.set_value(
                circ_1_start_radius + circ_1_growth_rate * time_tracker.get_value()
            )
        )

        circ_1 = always_redraw(
            lambda: Circle(radius=circ_1_radius_tracker.get_value(),
                           arc_center=np.array([circ_1_x, circ_1_y, circ_1_z]))
        )

        circ_2_x, circ_2_y, circ_2_z = 1.5, 0, 0
        circ_2_start_radius = 1.5
        circ_2_growth_rate = 1
        circ_2_radius_tracker = ValueTracker(circ_2_start_radius).add_updater(
            lambda tracker: tracker.set_value(
                circ_2_start_radius + circ_2_growth_rate * time_tracker.get_value()
            )
        )
        circ_2 = always_redraw(
            lambda: Circle(radius=circ_2_radius_tracker.get_value(),
                           arc_center=np.array([circ_2_x, circ_2_y, circ_2_z]))
        )

        self.add(circ_1_radius_tracker, circ_2_radius_tracker, circ_1, circ_2)
        self.wait()

        distance_between_circles = np.sqrt((circ_2_x - circ_1_x) ** 2 + (circ_2_y - circ_1_y) ** 2)
        # determine if the two circles have intersection point(s)
        # if ((circ_1_radius_tracker.get_value() +
        # circ_2_radius_tracker.get_value()) >= distance_between_circles) and ( distance_between_circles >= abs(
        # circ_1_radius_tracker.get_value() - circ_2_radius_tracker.get_value())):

        intersec_points = always_redraw(
            lambda: get_circ_intersec(circ_1_x, circ_1_y,
                                      circ_2_x, circ_2_y,
                                      circ_1_radius_tracker.get_value(), circ_2_radius_tracker.get_value())
        )

        intersec_point_1_path = TracedPath(intersec_points.submobjects[0].get_center, stroke_color=YELLOW,
                                           stroke_width=4)
        intersec_point_2_path = TracedPath(intersec_points.submobjects[1].get_center, stroke_color=YELLOW,
                                           stroke_width=4)
        self.add(intersec_points, intersec_point_1_path, intersec_point_2_path)

        # self.play(circ_1_radius_tracker.animate(run_time=2, rate_func=linear).set_value(4),
        #           circ_2_radius_tracker.animate(run_time=2, rate_func=linear).set_value(3))
        self.play(time_tracker.animate(run_time=2, rate_func=linear).set_value(2))
        self.wait()


class TestHorizontalSlice(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(
            x_range=[-5, 5],
            y_range=[-5, 5],
            z_range=[-5, 5],
            x_length=8,
            y_length=8,
            z_length=8
        )
        # make the axes completely transparent
        axes.set_opacity(opacity=0.3)

        self.set_camera_orientation(70 * DEGREES, -20 * DEGREES)
        self.add(axes)

        unit_y_length = axes.c2p(0, 1, 0)[1] - axes.c2p(0, 0, 0)[1]
        unit_z_length = axes.c2p(0, 0, 1)[2] - axes.c2p(0, 0, 0)[2]
        unit_y_vector = axes.c2p(0, 1, 0) - axes.c2p(0, 0, 0)
        unit_z_vector = axes.c2p(0, 0, 1) - axes.c2p(0, 0, 0)

        # add the time cone for the alpha phase (colored orange) and the beta phase (colored green)
        alpha_cone_radius = 2
        alpha_cone_height = 4
        alpha_rate_inverse = alpha_cone_radius / alpha_cone_height
        alpha_cone = Cone(
            base_radius=alpha_cone_radius * unit_y_length,
            height=alpha_cone_height * unit_z_length,
            stroke_color=ORANGE,
            fill_color=ORANGE,
            fill_opacity=0.3,
            stroke_opacity=0,
            resolution=50
        ).shift(alpha_cone_height * unit_z_vector)

        beta_cone_pos = 0.75
        beta_cone_time = 1
        beta_cone_radius = 0.75
        beta_cone_height = 3
        beta_rate_inverse = beta_cone_radius / beta_cone_height
        beta_cone = Cone(
            base_radius=beta_cone_radius * unit_y_length,
            height=beta_cone_height * unit_z_length,
            stroke_color=GREEN_E,
            fill_color=GREEN_E,
            fill_opacity=0.2,
            stroke_opacity=0,
            show_base=True,
            direction=-Z_AXIS,
            resolution=50
        ).shift(beta_cone_pos * unit_y_vector + beta_cone_time * unit_z_vector)

        self.add(alpha_cone, beta_cone)

        t_tracker = ValueTracker(0)
        plane = always_redraw(
            lambda:
            Square(side_length=10 * unit_y_length).shift(
                t_tracker.get_value() * unit_z_vector
            ).set_opacity(0.1).set_stroke(opacity=0)
        )
        circ_1 = always_redraw(
            lambda:
            Circle(radius=2 * unit_y_length - t_tracker.get_value() * alpha_rate_inverse * unit_z_length,
                   color=ORANGE, stroke_width=2).shift(
                t_tracker.get_value() * unit_z_vector
            )
        )
        circ_2 = always_redraw(
            lambda:
            Circle(radius=(t_tracker.get_value() - beta_cone_time) * beta_rate_inverse * unit_z_length,
                   color=GREEN, stroke_width=2).shift(
                beta_cone_pos * unit_y_vector + t_tracker.get_value() * unit_z_vector
            ) if t_tracker.get_value() > beta_cone_time else Circle(radius=0).shift(
                beta_cone_pos * unit_y_vector + t_tracker.get_value() * unit_z_vector
            )
        )
        circ_intersection = always_redraw(
            lambda:
            Intersection(circ_1, circ_2, color=YELLOW, fill_color=BLUE, fill_opacity=1, stroke_width=2).shift(
                t_tracker.get_value() * unit_z_vector
            )
        )

        self.add(circ_1, circ_2, circ_intersection, plane)
        self.play(
            t_tracker.animate.set_value(3),
            rate_func=linear,
            run_time=3
        )
        self.wait()
