from typing import Tuple

import time
import numpy as np
import pymunk
import networkx as nx
from scipy.spatial import ConvexHull, KDTree
from scipy.spatial.distance import cdist, pdist

from manim import *
from manim.opengl import *

# modified from https://github.com/Matheart/manim-physics/blob/main/src/manim_physics/rigid_mechanics.py
# For manim < 0.15.0
# from manim.mobject.opengl_compatibility import ConvertToOpenGL
# For manim >= 0.15.0
from manim.mobject.opengl.opengl_compatibility import ConvertToOpenGL

config.window_position = "0,100"


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
        self.nuc_cluster = set()
        self.nuc_cluster_bodies = set()
        self.nuc_cluster_graph = nx.empty_graph()
        self.particle_mobs = []
        self.convex_vertices = []
        self.convex_mob_indices = []
        self.anchor_pos = None
        self.critical_size = crit_size
        super().__init__(renderer=renderer, **kwargs)

    def setup(self):
        """Used internally"""
        self.add(self.space)
        # step the physics simulation in pymunk
        self.space.add_updater(_step)
        # update the convex hull
        self.add_updater(self._update_convex_hull)
        # update the graph network
        self.add_updater(self._update_cluster_graph)

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
                # store each mobject
                self.particle_mobs.append(mob)
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
        velocity_vec = b.body.velocity
        speed = velocity_vec.length
        b.set_color(interpolate_color(BLUE, RED,
                                      np.clip(speed / np.sqrt(v_max ** 2 + v_max ** 2), 0, 1)))
        b.move_to(x * RIGHT + y * UP)
        b.rotate(b.body.angle - b.angle)
        b.angle = b.body.angle
        # adjust the velocity based on the temperature
        velocity_direction = velocity_vec.normalized()
        if temp_tracker.get_value() == HEAT:
            b.body.apply_force_at_local_point(force=0.03 * velocity_direction)
            # b.body.apply_force_at_local_point(force=tuple(rng.uniform(-1, 1, size=(2, 1))))
        elif temp_tracker.get_value() == COOL:
            b.body.apply_force_at_local_point(force=-0.03 * velocity_direction)

        # set detachment behavior
        current_size = len(self.nuc_cluster)

        # set detachment probability before reaching the critical size
        if current_size < self.critical_size:
            detachment_prob = 1 / (90 * 7)
        # set detachment probability after reaching the critical size
        else:
            # after reaching the critical size, turn off probabilistic detachment
            detachment_prob = -1

        # probabilistic detachment
        if rng.uniform(0, 1) <= detachment_prob:
            self._detach(pymunk_space)

    def _get_convex_hull(self):
        """Update the convex hull of the nucleus"""
        nuc_cluster_bodies_lst = list(self.nuc_cluster_bodies)
        # only construct the convex hull when there are at least 3 particles in the nucleus
        if len(nuc_cluster_bodies_lst) >= 3:
            nuc_cluster_pos = [bod.position for bod in nuc_cluster_bodies_lst]
            convex_hull = ConvexHull(nuc_cluster_pos)
            self.convex_vertices = convex_hull.vertices
            # convert from cluster index to the mobject index
            self.convex_mob_indices = [list(nuc_cluster_bodies_lst[_].shapes)[0].collision_type
                                       if list(nuc_cluster_bodies_lst[_].shapes)[0].collision_type != 10000 else 0
                                       for _ in self.convex_vertices]
        return nuc_cluster_bodies_lst

    def _update_convex_hull(self, dt):
        self._get_convex_hull()

    def _update_cluster_graph(self, dt):
        # create an ordered list of the nucleation cluster indices
        nuc_cluster_indices = list(self.nuc_cluster)
        # get the body positions in the nucleation cluster
        nuc_cluster_pos = np.array([self.particle_mobs[mob_idx].get_center()[:2]
                                    for mob_idx in nuc_cluster_indices])

        if len(nuc_cluster_pos) > 0:
            # create a KDTree for nearest neighbor lookup
            kd_tree = KDTree(nuc_cluster_pos)
            # get all pairs of particles that are within max_node_dist
            # AHA: need to be careful about extra argument for updater function
            neighbor_pairs_raw = kd_tree.query_pairs(r=0.25)
            neighbor_pairs_converted = [(nuc_cluster_indices[i1], nuc_cluster_indices[i2])
                                        for i1, i2 in neighbor_pairs_raw]

            # initialize an empty graph
            cluster_graph = nx.Graph()
            # add the mobject index as nodes
            cluster_graph.add_nodes_from(nuc_cluster_indices)
            # add edges from neighbor_pairs
            cluster_graph.add_edges_from(neighbor_pairs_converted)
            self.nuc_cluster_graph = cluster_graph

    def _detach(self, pymunk_space):
        # get the bodies in the nucleus cluster and right before detachment selection
        nuc_cluster_bodies_at_detach = self._get_convex_hull()
        if isinstance(self.convex_vertices, np.ndarray):
            # choose a random convex hull vertex
            vertex_selected_idx = rng.choice(a=self.convex_vertices)
            body_to_remove = nuc_cluster_bodies_at_detach[vertex_selected_idx]
            cluster_remove_idx = list(body_to_remove.shapes)[0].collision_type
            dist_to_core = body_to_remove.position.get_distance(self.anchor_pos)
            # get all the joints connected to the selected vertex
            body_selected_joints = body_to_remove.constraints

            # detach only when the particle to be removed is not the initial nucleation site,
            # and if it isn't, it should not be too close to the initial nucleation site
            if (cluster_remove_idx not in {0, 10000}) and (dist_to_core > 0.24 * 2):
                self.particle_mobs[cluster_remove_idx].set_color(PURPLE)
                # remove the selected joint from the physics simulation
                for joint in body_selected_joints:
                    try:
                        pymunk_space.remove(joint)
                    # TODO: temporary fix, need to look into why the same joint is removed more than once
                    except AssertionError:
                        continue

                try:
                    self.nuc_cluster.remove(cluster_remove_idx)
                    self.nuc_cluster_bodies.remove(body_to_remove)
                except KeyError:
                    pass
                # get the outward radial direction
                v_vec = body_to_remove.position - self.anchor_pos
                # apply a force in the outward radial direction of moving
                v_dir = v_vec.normalized()
                body_to_remove.apply_force_at_local_point(force=0.01 * v_dir)
                # # update the convex hull after each detachment
                # self._get_convex_hull()

    def _begin(self, arbiter, space, data):
        """Callback function to set attachment at a given probability upon collision"""
        s1, s2 = arbiter.shapes
        i1, i2 = s1.collision_type, s2.collision_type
        # set attachment behavior
        current_size = len(self.nuc_cluster)
        # attach only when one of the two particles colliding is in the existing nucleus
        if (i1 in self.nuc_cluster or i2 in self.nuc_cluster) and 10000 not in {i1, i2}:
            # set attachment probability before reaching the critical size
            if current_size <= self.critical_size:
                attachment_prob = 0.7
            # set attachment probability after reaching the critical size
            else:
                attachment_prob = 0.7

            if rng.uniform(0, 1) <= attachment_prob:
                joint = pymunk.PinJoint(s1.body, s2.body)
                space.add(joint)
                self.nuc_cluster = self.nuc_cluster.union({i1, i2})
                self.nuc_cluster_bodies = self.nuc_cluster_bodies.union({s1.body, s2.body})
                print(f" size of nucleus: {len(self.nuc_cluster)}")

                # when above the critical size,
                # maintain the size by detaching particles for each attachment to the nuc_cluster
                if len(self.nuc_cluster) >= self.critical_size:
                    # set up a breakout time for the while loop
                    timeout = time.time() + 2
                    while (len(self.nuc_cluster) > self.critical_size - rng.choice(a=3)) and (
                            growth_tracker.get_value() == MAINTAIN_CRIT_SIZE):
                        self._detach(space)
                        if time.time() > timeout:
                            break
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
                self.particle_mobs.append(mob)
                mob.body.position = mob.get_x(), mob.get_y()
                self.anchor_pos = mob.body.position
                self.nuc_cluster = {collision_type}
                self.nuc_cluster_bodies = {mob.body}
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
# set the growth behavior
MAINTAIN_CRIT_SIZE = 1
GROW_BEYOND_CRIT_SIZE = 0
growth_tracker = ValueTracker(MAINTAIN_CRIT_SIZE)
# set seed to ensure reproducible initial particle positions
seed = 31415926
rng = np.random.default_rng(seed)


# run the following command in the terminal to start simulation
# manim -p -r 3840,2160 --fps 120 --renderer=opengl --use_projection_stroke_shaders --use_projection_fill_shaders --disable_caching --flush_cache --write_to_movie .\thermodynamics\collision_fixed_opengl.py CollisionFixedGL
# manim -p -r 1920,1080 --fps 120 --renderer=opengl --use_projection_stroke_shaders --use_projection_fill_shaders --disable_caching --flush_cache .\thermodynamics\collision_fixed_opengl.py CollisionFixedGL
class CollisionFixedGL(CollisionScene):
    def __init__(self, crit_size=30, renderer=None, **kwargs):
        super().__init__(crit_size, renderer, **kwargs)
        self.critical_size = 30
        self.nucleus_size_axis = None
        self.size_count_arr = np.array([[0, 1]])

    def construct(self):
        axes = Axes(
            x_range=[0, 10],
            y_range=[0, 10],
            x_length=6,
            y_length=6,
            tips=False,
            axis_config={"include_ticks": False}
        ).shift(LEFT)

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
            Dot(point=axes.c2p(*pos), radius=particle_r, stroke_width=1)
            for pos in particle_pos_lst
        ]
        particles[0].set_color(YELLOW)
        self.play(
            *[DrawBorderThenFill(part) for part in particles]
        )
        self.wait()

        # show the nucleus size counter
        nucleus_size_axis_x_max = 30
        nucleus_size_axis = Axes(
            x_range=[0, nucleus_size_axis_x_max, nucleus_size_axis_x_max],
            y_range=[0, NUM_PARTICLES, 10],
            x_length=2.5,
            y_length=2,
            tips=False,
            y_axis_config={
                "numbers_to_include": [0, self.critical_size, NUM_PARTICLES],
                "font_size": 25,
                "include_ticks": False
            }
        ).to_edge(RIGHT).shift(0.5 * LEFT)

        self.nucleus_size_axis = nucleus_size_axis

        # add y_ticks inplace
        y_axis = nucleus_size_axis.y_axis
        ticks = VGroup()
        for _ in [self.critical_size, NUM_PARTICLES]:
            ticks.add(y_axis.get_tick(_, y_axis.tick_size))
        y_axis.add(ticks)
        y_axis.ticks = ticks

        time_tracker = ValueTracker(0)

        def update_axes(ax):
            ax_to_become = ax
            current_time = time_tracker.get_value()
            if current_time >= nucleus_size_axis_x_max:
                ax_to_become = Axes(
                    x_range=[0, current_time, current_time],
                    y_range=[0, NUM_PARTICLES, 10],
                    x_length=2.5,
                    y_length=2,
                    tips=False,
                    y_axis_config={"numbers_to_include": [0, self.critical_size, NUM_PARTICLES],
                                   "font_size": 25, "include_ticks": False}
                ).to_edge(RIGHT).shift(0.5 * LEFT)

                y_axis_updated = ax_to_become.y_axis
                ticks_updated = VGroup()
                for _ in [self.critical_size, NUM_PARTICLES]:
                    ticks_updated.add(y_axis_updated.get_tick(_, y_axis_updated.tick_size))
                y_axis_updated.add(ticks_updated)
                y_axis_updated.ticks = ticks_updated
            ax.become(ax_to_become)
            self.nucleus_size_axis = ax_to_become

        nucleus_size_axis.add_updater(update_axes)

        # add horizontal dashed line at critical size
        crit_hline = OpenGLDashedVMobject(
            vmobject=Line(
                start=nucleus_size_axis.c2p(nucleus_size_axis.x_range[0], self.critical_size),
                end=nucleus_size_axis.c2p(nucleus_size_axis.x_range[1], self.critical_size),
                color=YELLOW, stroke_width=1)
        )

        # get the axis labels
        nucleus_size_x_label = nucleus_size_axis.get_x_axis_label(
            "t"
        ).scale(0.7)
        nucleus_size_y_label = nucleus_size_axis.get_y_axis_label(
            r"\text{Nucleus Size}", edge=LEFT
        ).rotate(PI / 2).scale(0.7).shift(1.8 * LEFT + 0.3 * DOWN)

        def update_time(tracker, dt):
            tracker.increment_value(dt)
            self.size_count_arr = np.append(
                self.size_count_arr,
                np.array([[time_tracker.get_value(), len(self.nuc_cluster)]]),
                axis=0
            )

        time_tracker.add_updater(update_time)

        nucleus_size_counter_dot = always_redraw(
            lambda:
            Dot(
                point=self.nucleus_size_axis.c2p(time_tracker.get_value(), len(self.nuc_cluster)),
                radius=0.06, stroke_width=1
            )
        )
        counter_dot_path = always_redraw(
            lambda:
            OpenGLVMobject(
                stroke_width=2, color=WHITE
            ).set_points_as_corners(
                [self.nucleus_size_axis.c2p(*coord) for coord in self.size_count_arr]
            )
        )
        self.add(time_tracker)
        time_tracker.suspend_updating()
        self.play(
            FadeIn(nucleus_size_axis),
            Write(nucleus_size_x_label),
            Write(nucleus_size_y_label),
            DrawBorderThenFill(nucleus_size_counter_dot),
            FadeIn(counter_dot_path)
        )
        self.play(Create(crit_hline))
        self.wait()

        # prepare the collision simulation
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

        time_tracker.resume_updating()

        def draw_convex_hull(poly, dt):
            """Helper function to draw the convex hull"""
            _ = dt + 1  # TODO: just a dummy line to test the effect of dt
            if isinstance(self.convex_vertices, np.ndarray):
                if len(self.nuc_cluster) < self.critical_size:
                    poly.become(Polygon(
                        *[self.particle_mobs[idx].get_center() for idx in self.convex_mob_indices],
                        color=PURPLE
                    ))
                else:
                    poly.become(Polygon(
                        *[self.particle_mobs[idx].get_center() for idx in self.convex_mob_indices],
                        color=YELLOW
                    ))
            else:
                poly.become(Polygon(particles[0].get_center(),
                                    particles[0].get_center(),
                                    particles[0].get_center(),
                                    stroke_opacity=0, fill_opacity=0))

        # draw the convex hull
        convex_hull = Polygon(particles[0].get_center(),
                              particles[0].get_center(),
                              particles[0].get_center(), stroke_opacity=0, fill_opacity=0)
        convex_hull.add_updater(draw_convex_hull)
        self.add(convex_hull)

        def draw_graph_network(graph, dt):
            """Helper function to draw the convex hull"""
            _ = dt
            updated_networkx_graph = Graph.from_networkx(self.nuc_cluster_graph,
                                                         layout={
                                                             node: self.particle_mobs[node].get_center()
                                                             for node in self.nuc_cluster_graph.nodes}
                                                         ).get_group_class()
            graph.become(updated_networkx_graph)

        # draw the nucleation cluster graph
        print("Before graph")
        cluster_graph = always_redraw(
            lambda:
            Graph.from_networkx(self.nuc_cluster_graph,
                                layout={
                                    node: self.particle_mobs[node].get_center()
                                    for node in self.nuc_cluster_graph.nodes}
                                )
        )
        # cluster_graph.add_updater(draw_graph_network)
        print("Updating graph")
        self.add(cluster_graph)

        # start the collision simulation
        self.wait(8)

        for _ in range(5):
            temp_tracker.set_value(HEAT)
            self.wait(1)

            temp_tracker.set_value(MAINTAIN)
            self.wait(4)

        # let the nucleus grow beyond the critical size
        growth_tracker.set_value(GROW_BEYOND_CRIT_SIZE)

        for _ in range(5):
            temp_tracker.set_value(HEAT)
            self.wait(1)

            temp_tracker.set_value(MAINTAIN)
            self.wait(4)

        self.wait(4)

        temp_tracker.set_value(COOL)
        self.wait(4)


class OpenGLTest(Scene):
    def construct(self):
        circ = Circle(point=LEFT, radius=1, stroke_width=10, fill_color=YELLOW, fill_opacity=1)
        self.play(DrawBorderThenFill(circ))
        self.wait(10)
