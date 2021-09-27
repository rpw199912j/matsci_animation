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
