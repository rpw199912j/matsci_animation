from typing import Tuple

import numpy as np
import pymunk
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist

from manim import *

# modified from https://github.com/Matheart/manim-physics/blob/main/src/manim_physics/rigid_mechanics.py
# For manim < 0.15.0
from manim.mobject.opengl_compatibility import ConvertToOpenGL


class Space(Mobject, metaclass=ConvertToOpenGL):
    def __init__(self, gravity: Tuple[float, float] = (0, -9.81), **kwargs):
        """An Abstract object for gravity.
        Parameters
        ----------
        gravity
            The direction and strength of gravity.
        """
        super().__init__(**kwargs)
        self.space = pymunk.Space()
        self.space.gravity = gravity
        self.space.sleep_time_threshold = 5


class CollisionScene(Scene):
    GRAVITY: Tuple[float, float] = 0, 0

    def __init__(self, crit_size=30, renderer=None, **kwargs):
        """A basis scene for all of rigid mechanics. The gravity vector
        can be adjusted with ``self.GRAVITY``.
        """
        self.space = Space(gravity=self.GRAVITY)
        self.joints = []
        self.nuc_cluster = {0}
        self.anchor_pos = None
        self.critical_size = crit_size
        super().__init__(renderer=renderer, **kwargs)

    def setup(self):
        """Used internally"""
        self.add(self.space)
        self.space.add_updater(_step)

    def add_body(self, body: Mobject):
        """Bodies refer to pymunk's object.
        This method ties Mobjects to their Bodies.
        """
        if body.body != self.space.space.static_body:
            self.space.space.add(body.body)
        self.space.space.add(body.shape)

    def make_rigid_body(
            self,
            *mobs: Mobject,
            elasticity: float = 0.8,
            density: float = 1,
            friction: float = 0.8,
            velocity: tuple = (0, 0),
            collision_type: int = 1
    ):
        """Make any mobject movable by gravity.
        Equivalent to ``Scene``'s ``add`` function.
        Parameters
        ----------
        mobs
            The mobs to be made rigid.
        elasticity
        density
        friction
            The attributes of the mobjects in regards to
            interacting with other rigid and static objects.
        velocity
        collision_type
        """
        for mob in mobs:
            if isinstance(mob, VGroup):
                return self.make_rigid_body(*mob)
            if not hasattr(mob, "body"):
                parts = mob.family_members_with_points()
                for p in parts:
                    self.add(p)
                    p.body = pymunk.Body()
                    p.body.position = p.get_x(), p.get_y()
                    p.body.velocity = velocity
                    # limit the velocity to prevent tunneling
                    p.body.velocity_func = self._limit_velocity
                    get_angle(p)
                    if not hasattr(p, "angle"):
                        p.angle = 0
                    p.body.angle = p.angle
                    get_shape(p)
                    p.shape.density = density
                    p.shape.elasticity = elasticity
                    p.shape.friction = friction
                    p.shape.collision_type = collision_type
                    p.spacescene = self

                    self.add_body(p)
                    p.add_updater(self._simulate)

            else:
                if mob.body.is_sleeping:
                    mob.body.activate()

    @staticmethod
    def _limit_velocity(body, gravity, damping, dt):
        """Set an upper limit for the speed to prevent tunneling"""
        # TODO: need to manually tweak max_speed (currently at v_max=3, it does not work)
        max_speed = 100
        pymunk.Body.update_velocity(body, gravity, damping, dt)
        current_speed = body.velocity.length
        if current_speed > max_speed:
            scale = max_speed / current_speed
            body.velocity = body.velocity * scale

    def _simulate(self, b):
        """Simulate the physics after each dt"""
        # get the space
        pymunk_space = self.space.space
        # get the default collision handler
        coll_handler = pymunk_space.add_default_collision_handler()
        # set custom collision behavior
        coll_handler.begin = self._begin
        x, y = b.body.position
        # set color based on speed
        vx, vy = b.body.velocity
        speed = np.sqrt(vx ** 2 + vy ** 2)
        b.set_color(interpolate_color(BLUE, RED,
                                      np.clip(speed / np.sqrt(v_max ** 2 + v_max ** 2), 0, 1)))
        b.move_to(x * RIGHT + y * UP)
        b.rotate(b.body.angle - b.angle)
        b.angle = b.body.angle
        # adjust the velocity based on the temperature
        velocity_vex = np.array([vx, vy])
        velocity_direction = velocity_vex / np.linalg.norm(velocity_vex)
        if temp_tracker.get_value() == HEAT:
            b.body.apply_force_at_local_point(force=tuple(0.03 * velocity_direction))
            # b.body.apply_force_at_local_point(force=tuple(rng.uniform(-1, 1, size=(2, 1))))
        elif temp_tracker.get_value() == COOL:
            b.body.apply_force_at_local_point(force=tuple(-0.03 * velocity_direction))

        # set detachment behavior
        current_size = len(self.nuc_cluster)

        # set detachment probability before reaching the critical size
        if current_size <= self.critical_size:
            detachment_prob = 1 / (90 * 7)
        # set detachment probability after reaching the critical size
        else:
            # linearly increase detachment prob based on nucleus size
            detachment_prob = np.interp(x=current_size,
                                        xp=[self.critical_size, 100],
                                        fp=[1 / (90 * 8),
                                            1 / (90 * 6)])
            # detachment_prob = -1
        if rng.uniform(0, 1) <= detachment_prob:
            existing_joints = self.joints
            if existing_joints:
                # randomly choose an arbitrary joint
                joint_selected_idx = rng.choice(a=len(existing_joints))
                joint_selected = existing_joints[joint_selected_idx]
                # get the two particles at the two ends of the selected joint
                b1, b2 = joint_selected.a, joint_selected.b
                # calculate each particle's distance to the initial nucleation site
                dist_1 = b1.position.get_distance(self.anchor_pos)
                dist_2 = b2.position.get_distance(self.anchor_pos)
                if dist_1 > dist_2:
                    cluster_remove_idx = list(b1.shapes)[0].collision_type
                    body_to_remove = b1
                else:
                    cluster_remove_idx = list(b2.shapes)[0].collision_type
                    body_to_remove = b2
                # remove the selected joint from the physics simulation
                pymunk_space.remove(joint_selected)
                existing_joints.remove(joint_selected)
                # apply force upon detachment
                if cluster_remove_idx != 0:
                    try:
                        self.nuc_cluster.remove(cluster_remove_idx)
                    except KeyError:
                        pass
                    # get the velocity direction
                    v_vec = np.array([body_to_remove.velocity[0], body_to_remove.velocity[1]])
                    v_mag = np.linalg.norm(v_vec)
                    # if non-zero velocity, apply a force in the direction of moving
                    if v_mag > 0:
                        v_dir = v_vec / v_mag
                        body_to_remove.apply_force_at_local_point(force=tuple(0.01 * v_dir))

    def _begin(self, arbiter, space, data):
        """Callback function to set attachment at a given probability upon collision"""
        b1, b2 = arbiter.shapes
        i1, i2 = b1.collision_type, b2.collision_type
        # set attachment behavior
        current_size = len(self.nuc_cluster)
        # attach only when one of the two particles colliding is in the existing nucleus
        if i1 in self.nuc_cluster or i2 in self.nuc_cluster:
            # set attachment probability before reaching the critical size
            if current_size <= self.critical_size:
                attachment_prob = 0.7
            # set attachment probability after reaching the critical size
            else:
                attachment_prob = 0.7

            if rng.uniform(0, 1) <= attachment_prob:
                joint = pymunk.PinJoint(b1.body, b2.body)
                space.add(joint)
                self.joints.append(joint)
                self.nuc_cluster = self.nuc_cluster.union({i1, i2})
                print(f" size of nucleus: {len(self.nuc_cluster)}")
        return True

    def make_static_body(
            self, *mobs: Mobject, elasticity: float = 1, friction: float = 0.8,
            collision_type: int = 10000
    ) -> None:
        """Make any mobject interactable by rigid objects.
        Parameters
        ----------
        mobs
            The mobs to be made static.
        elasticity
        friction
            The attributes of the mobjects in regards to
            interacting with rigid objects.
        collision_type
        """
        for mob in mobs:
            if isinstance(mob, VGroup or Group):
                return self.make_static_body(*mob)
            mob.body = self.space.space.static_body
            get_shape(mob)
            mob.shape.elasticity = elasticity
            mob.shape.friction = friction
            mob.shape.collision_type = collision_type
            if collision_type == 0:
                self.anchor_pos = mob.body.position
            self.add_body(mob)

    def stop_rigidity(self, *mobs: Mobject) -> None:
        """Stop the mobjects rigidity"""
        for mob in mobs:
            if isinstance(mob, VGroup or Group):
                self.stop_rigidity(*mob)
            if hasattr(mob, "body"):
                mob.body.sleep()


def _step(space, dt):
    space.space.step(dt)


def get_shape(mob: VMobject) -> None:
    """Obtains the shape of the body from the mobject"""
    if isinstance(mob, Circle):
        mob.shape = pymunk.Circle(body=mob.body, radius=mob.radius)
    elif isinstance(mob, Line):
        mob.shape = pymunk.Segment(
            mob.body,
            (mob.get_start()[0], mob.get_start()[1]),
            (mob.get_end()[0], mob.get_end()[1]),
            mob.stroke_width - 3.95,
        )
    elif issubclass(type(mob), Rectangle):
        width = np.linalg.norm(mob.get_vertices()[1] - mob.get_vertices()[0])
        height = np.linalg.norm(mob.get_vertices()[2] - mob.get_vertices()[1])
        mob.shape = pymunk.Poly.create_box(mob.body, (width, height))
    elif issubclass(type(mob), Polygram):
        vertices = [(a, b) for a, b, c in mob.get_vertices() - mob.get_center()]
        mob.shape = pymunk.Poly(mob.body, vertices)
    else:
        mob.shape = pymunk.Poly.create_box(mob.body, (mob.width, mob.height))


def get_angle(mob: VMobject) -> None:
    """Obtains the angle of the body from the mobject.
    Used internally for updaters.
    """
    if issubclass(type(mob), Polygon):
        vec1 = mob.get_vertices()[0] - mob.get_vertices()[1]
        vec2 = type(mob)().get_vertices()[0] - type(mob)().get_vertices()[1]
        mob.angle = angle_between_vectors(vec1, vec2)
    elif isinstance(mob, Line):
        mob.angle = mob.get_angle()


v_max = 3  # set the maximum velocity magnitude in either horizontal or vertical direction
NUM_PARTICLES = 100
MAINTAIN = 10
HEAT = 20
COOL = 5
temp_tracker = ValueTracker(MAINTAIN)
# set seed to ensure reproducible initial particle positions
seed = 31415926
rng = np.random.default_rng(seed)


# run the following command in the terminal to start simulation
# manim -p -r 1920,1080 --fps 90 --disable_caching --flush_cache <path to collision_fixed.py> CollisionFixed
class CollisionFixed(CollisionScene):
    def construct(self):
        axes = Axes(
            x_range=[0, 10],
            y_range=[0, 10],
            x_length=6,
            y_length=6,
            tips=False,
            axis_config={"include_ticks": False}
        )

        # set up the bounding box
        bounding_right = Line(
            start=axes.c2p(axes.x_range[1], axes.y_range[0]),
            end=axes.c2p(axes.x_range[1], axes.y_range[1])
        )
        bounding_top = Line(
            start=axes.c2p(axes.x_range[0], axes.y_range[1]),
            end=axes.c2p(axes.x_range[1], axes.y_range[1])
        )
        bounding_left = Line(
            start=axes.c2p(axes.x_range[0], axes.y_range[1]),
            end=axes.c2p(axes.x_range[0], axes.y_range[0])
        )
        bounding_bottom = Line(
            start=axes.c2p(axes.x_range[0], axes.y_range[0]),
            end=axes.c2p(axes.x_range[1], axes.y_range[0])
        )
        bounding_box = VGroup(bounding_right, bounding_top, bounding_left, bounding_bottom)
        self.add(bounding_box)
        self.wait()

        # add the particles
        num_particles = NUM_PARTICLES
        rad = 0.2
        particle_r = rad * axes.get_x_unit_size()
        particle_pos_lst = np.array([[5, 5]])
        # make sure the circles do not overlap
        while len(particle_pos_lst) < num_particles:
            # generate a random coordinate
            pos = np.array([rng.uniform(0 + rad, 10 - rad, size=2)])
            # calculate the distance to all the existing circles
            dist_arr = cdist(pos, particle_pos_lst).flatten()
            # check if any distance is smaller than the defined circle diameter
            is_overlapping = (dist_arr < rad * 2).sum()
            if is_overlapping > 0:
                continue
            else:
                particle_pos_lst = np.append(particle_pos_lst, pos, axis=0)

        particles = [
            Dot(point=axes.c2p(*pos), radius=particle_r)
            for pos in particle_pos_lst
        ]
        particles[0].set_color(YELLOW)
        self.play(
            *[DrawBorderThenFill(part) for part in particles]
        )
        self.wait()

        # start the collision simulation
        self.make_static_body(
            bounding_box,
            friction=0
        )
        self.make_static_body(
            particles[0],
            friction=0,
            collision_type=0
        )

        for _, part in enumerate(particles[1:]):
            self.make_rigid_body(
                part,
                elasticity=1,
                velocity=tuple(rng.uniform(-v_max, v_max, size=(2, 1))),
                friction=0,
                collision_type=_ + 1
            )

        self.wait(8)

        for _ in range(2):
            temp_tracker.set_value(HEAT)
            self.wait(1)

            temp_tracker.set_value(MAINTAIN)
            self.wait(4)

        self.wait(4)

        temp_tracker.set_value(COOL)
        self.wait(4)

        # temp_tracker.set_value(COOL)
        # self.wait(15)
