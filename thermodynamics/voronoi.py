import numpy as np
import shapely.geometry as sg
from scipy.spatial import Voronoi

from manim import *


# voronoi_finite_polygons_2d source:
# https://stackoverflow.com/questions/20515554/colorize-voronoi-diagram/20678647#20678647
def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:, 1] - c[1], vs[:, 0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)


# set seed for reproducibility
seed = 31415926
rng = np.random.default_rng(seed)


class BackgroundVoronoi(Scene):
    def construct(self):
        # generate 10 random points as initial nucleation site
        points = 9 * rng.random(size=(10, 2))

        # compute Voronoi tesselation
        vor = Voronoi(points)

        # get the finite voronoi regions
        regions, vertices = voronoi_finite_polygons_2d(vor, radius=20)

        # construct the shapely polygons for each region
        shapely_polys = [sg.Polygon(vertices[region]) for region in regions]
        # construct the bounding box
        x_min, x_max, y_min, y_max = 0, 10, 0, 10
        bounding_box = sg.Polygon(np.array([[x_min, y_min],
                                            [x_max, y_min],
                                            [x_max, y_max],
                                            [x_min, y_max]]))
        # check the intersection with the bounding box
        truncated_polys = []
        for poly in shapely_polys:
            # check if the bounding box fully contains each region
            is_contained = bounding_box.contains(poly)
            new_poly = poly
            if not is_contained:
                # get the truncated polygons with the bounding box
                new_poly = poly.intersection(bounding_box)
            truncated_polys.append(new_poly)

        # get the vertices for each polygon
        truncated_polys_vertices = [np.asarray(poly.exterior.coords)
                                    for poly in truncated_polys]

        # define the axes to help position the mobjects
        axes = Axes(
            x_range=[x_min, x_max],
            y_range=[y_min, y_max],
            x_length=6,
            y_length=6,
            tips=False,
            axis_config={"include_ticks": False}
        ).shift(LEFT * 3)

        # set up the bounding box
        bounding_right = Line(
            start=axes.c2p(axes.x_range[1], axes.y_range[0]),
            end=axes.c2p(axes.x_range[1], axes.y_range[1]),
            z_index=10
        )
        bounding_top = Line(
            start=axes.c2p(axes.x_range[0], axes.y_range[1]),
            end=axes.c2p(axes.x_range[1], axes.y_range[1]),
            z_index=10
        )
        bounding_left = Line(
            start=axes.c2p(axes.x_range[0], axes.y_range[1]),
            end=axes.c2p(axes.x_range[0], axes.y_range[0]),
            z_index=10
        )
        bounding_bottom = Line(
            start=axes.c2p(axes.x_range[0], axes.y_range[0]),
            end=axes.c2p(axes.x_range[1], axes.y_range[0]),
            z_index=10
        )
        bounding_box = VGroup(bounding_right, bounding_top, bounding_left, bounding_bottom)
        self.add(bounding_box)
        # define a mask around the bounding box
        whole_screen_box = Rectangle(
            width=config.frame_width,
            height=config.frame_height,
            color=YELLOW,
            fill_opacity=0.2
        )
        poly_to_subtract = Polygon(
            *[bounding_box.get_corner(corner) for corner in [UR, UL, DL, DR]]
        )
        mask = Difference(whole_screen_box, poly_to_subtract, fill_opacity=1)
        mask.set_color(BLACK)
        self.add(mask)
        self.wait()
        # add the text description
        msg_1 = Tex(
            "Phase transformation", font_size=1.2 * DEFAULT_FONT_SIZE,
        ).next_to(bounding_box, RIGHT, buff=3 * DEFAULT_MOBJECT_TO_MOBJECT_BUFFER).align_to(bounding_box, UP)

        self.play(
            Write(msg_1),
            run_time=1
        )
        self.wait()

        msg_2 = Tex(
            "Nucleation", "+", "Growth", arg_separator=" ",
        ).next_to(msg_1, DOWN, buff=4 * DEFAULT_MOBJECT_TO_MOBJECT_BUFFER).align_to(msg_1, LEFT)
        msg_2[0].set_color(YELLOW)
        msg_2[2].set_color(ORANGE)

        # convert the polygon vertices into the axes basis
        converted_polys_vertices = [
            [axes.c2p(*vertex).tolist() for vertex in poly_vertices]
            for poly_vertices in truncated_polys_vertices
        ]
        # convert to manim polygons
        manim_polys = [Polygon(*vertex_arr, color=ORANGE, fill_opacity=0.2, stroke_width=4)
                       for vertex_arr in converted_polys_vertices]

        # convert the initial sites coordinates into the axes basis
        converted_sites_coords = [
            axes.c2p(*site_coord) for site_coord in points
        ]
        nucleation_sites = [
            Dot(point=coord, color=YELLOW, z_index=-10) for coord in converted_sites_coords
        ]

        growth_rate = 0.2
        time_tracker = ValueTracker(0)

        # show the initial nucleation sites
        self.play(
            Write(msg_2[0]),
            *[GrowFromCenter(mob) for mob in nucleation_sites]
        )
        self.wait()

        nucleation_sites_copy = [mob.copy() for mob in nucleation_sites]
        for _, dot in enumerate(nucleation_sites_copy):
            # AHA: this is how you add updaters in a for-loop
            def grow_particle(mob, dt, idx=_):
                new_dot = Dot(
                    point=converted_sites_coords[idx],
                    color=ORANGE,
                    radius=DEFAULT_DOT_RADIUS + time_tracker.get_value() * growth_rate,
                    stroke_width=4, fill_opacity=0.2
                )
                truncated_dot = Intersection(
                    new_dot, manim_polys[idx],
                    color=ORANGE, stroke_width=4, fill_opacity=0.2,
                    z_index=-10
                )
                mob.become(truncated_dot)

            dot.add_updater(grow_particle)

        self.add(*nucleation_sites_copy)
        self.play(
            Write(msg_2[1:]),
            time_tracker.animate(rate_func=linear, run_time=5).set_value(20),
        )
        self.wait()
        for voronoi_region in nucleation_sites_copy:
            voronoi_region.clear_updaters()

        # focus on the nucleation part
        self.play(
            Circumscribe(msg_2[0], run_time=1.5)
        )
        self.wait()

        # zoom in on one of the nuclei
        combined_vor_group = VGroup(*(nucleation_sites + nucleation_sites_copy), z_index=-10)
        center_idx = 8
        scale_factor = 10
        combined_vor_group_copy = combined_vor_group.copy().scale(scale_factor)
        shift_vec = axes.get_center() - combined_vor_group_copy[center_idx].get_center()
        combined_vor_group_without_center_mob = VGroup(
            *[mob for _, mob in enumerate(combined_vor_group) if _ != center_idx]
        )
        combined_vor_group_copy_shifted = combined_vor_group_copy.copy().shift(shift_vec)
        self.play(
            Transform(combined_vor_group[center_idx], combined_vor_group_copy_shifted[center_idx]),
            FadeOut(msg_2[1:], shift=RIGHT),
            FadeOut(combined_vor_group_without_center_mob,
                    scale=scale_factor,
                    shift=shift_vec)
        )
        self.wait()

        # focus on the energetic description of a stable nucleus
        target_nucleus = combined_vor_group[center_idx]
        target_nucleus_initial_r = (target_nucleus.copy().get_right() - target_nucleus.copy().get_center())[0]
        radius_line = always_redraw(
            lambda:
            DashedLine(
                start=target_nucleus.get_center(),
                end=target_nucleus.get_right(),
                color=ORANGE
            )
        )
        radius_line_label = always_redraw(
            lambda:
            MathTex("r", color=ORANGE).scale(
                (target_nucleus.get_right() - target_nucleus.get_center())[0] / target_nucleus_initial_r
            ).next_to(radius_line, UP, buff=0.3 * DEFAULT_MOBJECT_TO_MOBJECT_BUFFER)
        )
        self.play(
            LaggedStart(
                Create(radius_line, rate_func=smooth),
                Write(radius_line_label),
                lag_ratio=0.9, run_time=1
            )
        )
        self.wait()

        # emphasize the competition between the volume and surface gibbs energy terms
        msg_3 = MathTex(
            r"&\text{Gibbs energy: }", r"G\\",
            r"&\Delta G_{\text{tot}}", "(r)", "=", r"\Delta G_{\text{vol}}", "+", r"\Delta G_{\text{sur}}"
        ).next_to(msg_2, DOWN, buff=4 * DEFAULT_MOBJECT_TO_MOBJECT_BUFFER).align_to(msg_2, LEFT)
        self.play(
            Write(msg_3)
        )
        self.wait()

        self.play(
            Indicate(target_nucleus, color=PURPLE, run_time=2),
            Indicate(msg_3[5], color=PURPLE, run_time=2)
        )
        self.wait()

        self.play(
            Circumscribe(target_nucleus, shape=Circle,
                         buff=-0.3, color=PURPLE, time_width=1,
                         run_time=2),
            Indicate(msg_3[7], color=PURPLE, run_time=2)
        )
        self.wait()

        self.play(
            target_nucleus.animate.scale(2)
        )
        self.wait()
        self.play(
            target_nucleus.animate.scale(0.001)
        )
        self.wait()
        self.play(
            target_nucleus.animate.scale(1 / 0.002)
        )
        self.wait()

        # show the density assumption
        msg_4 = MathTex(
            r"&\text{Assumption: }", r"\text{no density change}\\",
            r"&G/\text{mol} = G/\text{vol}", font_size=40
        ).next_to(msg_3, DOWN, buff=4 * DEFAULT_MOBJECT_TO_MOBJECT_BUFFER).align_to(msg_3, LEFT)
        self.play(
            Write(msg_4)
        )
        self.wait()

        # transition into the next scene
        all_mobs = [mob for mob in self.mobjects]
        all_mobs.remove(msg_3)
        msg_3_partial = MathTex(
            r"\Delta G_{\text{tot}}", "&=", r"\Delta G_{\text{vol}}", "+", r"\Delta G_{\text{sur}}\\",
            "&=", r"\frac{4}{3}\pi r^3\left(G_\beta-G_\alpha\right)\\",
            "&+",
            r"4\pi r^2\gamma"
        ).to_corner(UR)
        self.play(
            *[FadeOut(mob) for mob in all_mobs + [msg_3[:2]]],
            ReplacementTransform(msg_3[2:], msg_3_partial[:5])
        )
        self.wait()
