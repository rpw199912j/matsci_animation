from manim import *
from manim_editor import PresentationSectionType


def make_elements():  # only setting up the mobjects
    dots = VGroup(Dot(), Dot(), Dot(), Dot(), Dot(), Dot(), Dot(), z_index=0)
    dots.arrange(buff=0.7).scale(2).set_color(BLUE)
    dots[0].set_color(ORANGE)
    dots[-1].set_color(ORANGE)
    moving_dot = Dot(color=ORANGE, z_index=1).scale(2.5)
    moving_dot.move_to(dots[0])
    path = VGroup()
    path.add_updater(lambda x: x.become(Line(dots[0], moving_dot, stroke_width=10, z_index=1, color=ORANGE)))
    return dots, moving_dot, path


class MinimalPresentationExample(Scene):
    def construct(self):
        dots, moving_dot, path = make_elements()
        self.add(dots, moving_dot, path)

        self.next_section("A", PresentationSectionType.NORMAL)
        self.play(moving_dot.animate.move_to(dots[1]), rate_func=linear)

        self.next_section("A.1", PresentationSectionType.SUB_NORMAL)
        self.play(moving_dot.animate.move_to(dots[2]), rate_func=linear)

        self.next_section("B", PresentationSectionType.SKIP)
        self.play(moving_dot.animate.move_to(dots[3]), rate_func=linear)

        self.next_section("C", PresentationSectionType.LOOP)
        self.play(moving_dot.animate.move_to(dots[4]), rate_func=linear)

        self.next_section("D", PresentationSectionType.COMPLETE_LOOP)
        self.play(moving_dot.animate.move_to(dots[5]), rate_func=linear)

        self.next_section("E", PresentationSectionType.NORMAL)
        self.play(moving_dot.animate.move_to(dots[6]), rate_func=linear)


class CubeSymmetry(ThreeDScene):
    def construct(self):
        # define a 3D axes
        axes_3d = ThreeDAxes(
            tips=False
        )

        # define a cube with side length 2 placed at the origin
        cube = Cube(stroke_color=YELLOW, stroke_width=3)
        # define a line that aligns with one the edges
        line = Line3D(
            start=np.array([1, -1, 1]),
            end=np.array([1, 1, 1]),
            stroke_color=PURPLE
        )

        self.add(axes_3d)
        self.next_section("F", PresentationSectionType.NORMAL)

        self.play(
            FadeIn(cube)
        )
        self.next_section("G", PresentationSectionType.NORMAL)

        self.move_camera(phi=(90 - 35.26) * DEGREES, theta=-45 * DEGREES)
        self.wait()

        self.play(
            Create(line)
        )
        self.next_section("F", PresentationSectionType.NORMAL)

        for _ in range(3):
            self.play(
                Rotate(VGroup(cube, line), angle=120 * DEGREES, axis=np.array([1, -1, 1]))
            )
            self.next_section(f"{_}", PresentationSectionType.NORMAL)
