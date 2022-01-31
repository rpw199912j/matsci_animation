from manim import *


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
        self.wait()

        self.play(
            FadeIn(cube)
        )
        self.wait()

        self.move_camera(phi=(90 - 35.26) * DEGREES, theta=-45 * DEGREES)
        self.wait()

        self.play(
            Create(line)
        )
        self.wait()

        for _ in range(3):
            self.play(
                Rotate(VGroup(cube, line), angle=120 * DEGREES, axis=np.array([1, -1, 1]))
            )
            self.wait()
