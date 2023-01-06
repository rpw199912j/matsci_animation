import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from manim import *
from tqdm import tqdm

# set seed of the random number generator for reproducibility
seed = 31415926
rng = np.random.default_rng(seed)

# set manim configs
config.flush_cache = True
config.disable_caching = True
config.frame_rate = 60
config.frame_size = (720, 720)

# set the matplotlib color bar
cmap = mpl.cm.RdYlBu_r
norm = mpl.colors.Normalize(vmin=0, vmax=1)


class CahnHilliard2D(Scene):
    def __init__(self, nx=200, ny=200, steps=5000):
        super().__init__()
        # define the sampling resolution
        self.n_rows = ny
        self.n_cols = nx
        self.n_steps = steps
        # set the initial uniform composition
        self.c0 = 0.5
        # initialize the normalized current, the next and the Laplacian concentration grid
        self.c = np.zeros((self.n_rows, self.n_cols))  # holds the 2D concentration array at t
        self.c_new = np.zeros((self.n_rows, self.n_cols))  # holds the 2D concentration array at t+dt
        self.c_lap = np.zeros((self.n_rows, self.n_cols))
        # define the physical properties
        # self.mob = 1e-9  # mobility [m^2/sec]
        # self.sigma = 0.1  # interfacial energy [J/m^2]
        # self.wid = 1e-8  # interface width [m]
        # self.h = self.wid / 5
        # self.A = 3 * self.sigma / self.wid
        # self.K = 3 * self.sigma * self.wid
        # self.alpha = 1e15
        # self.dt = self.alpha * self.h ** 4
        self.mob = 0.01  # mobility [m^2/sec]
        self.h = 1e-2
        self.A = 100
        self.K = 1e-2
        self.dt = 1e-6

    def _init_c(self, beta=0.01):
        """Initialize the concentration array at t=0"""
        self.c = self.c0 + beta * rng.uniform(low=-1, high=1, size=(self.n_rows, self.n_cols))

    def _update_c_lap(self):
        """Update the Laplacian of the current concentration grid"""
        for i in range(self.n_rows):
            for j in range(self.n_cols):
                c_center = self.c[i, j]
                # get the 4 surrounding neighbors' concentration with periodicity
                c_right = self.c[i, (j + 1) % self.n_cols]
                c_left = self.c[i, (j + self.n_cols - 1) % self.n_cols]
                c_lower = self.c[(i + 1) % self.n_rows, j]
                c_upper = self.c[(i + self.n_rows - 1) % self.n_rows, j]
                self.c_lap[i, j] = (c_lower + c_upper + c_right + c_left - 4 * c_center) / (self.h ** 2)

    def _compute_node_u(self, node_i, node_j):
        """Compute each u_i,j value of the node (inside Laplacian)"""
        # get the concentration of the current node
        c_node = self.c[node_i, node_j]
        # get the Laplacian of the current node's concentration
        c_lap_node = self.c_lap[node_i, node_j]

        # compute df/dphi terms
        df_dphi = self.A * (2 * c_node - 6 * c_node ** 2 + 4 * c_node ** 3)
        # compute laplacian terms
        lap = 2 * self.K * c_lap_node
        
        return df_dphi - lap
    
    def _compute_node_lap(self, node_i, node_j):
        """Compute the laplacian term at each node i, j"""
        # get the current node's u_ij
        u_node = self._compute_node_u(node_i, node_j)
        # get the 4 surrounding neighbors' of u_ij with periodicity
        u_right = self._compute_node_u(node_i, (node_j + 1) % self.n_cols)
        u_left = self._compute_node_u(node_i, (node_j + self.n_cols - 1) % self.n_cols)
        u_lower = self._compute_node_u((node_i + 1) % self.n_rows, node_j)
        u_upper = self._compute_node_u((node_i + self.n_rows - 1) % self.n_rows, node_j)

        # compute the laplacian
        u_lap = (u_lower + u_upper + u_right + u_left - 4 * u_node) / (self.h ** 2)
        return u_lap

    def _update_c(self):
        """Update the concentration array after a single time step"""
        # update the laplacian
        self._update_c_lap()
        # update the normalized concentration grid
        for i in range(self.n_rows):
            for j in range(self.n_cols):
                self.c_new[i, j] = self.mob * self._compute_node_lap(i, j)

        self.c_new = self.c + self.dt * self.c_new
        # update c
        self.c[:, :] = self.c_new[:, :]

    def conc_to_color(self):
        """Convert the concentration arrays to a uint8 array in the rgba mode"""
        # get the image from matplotlib
        img = plt.imshow(self.c, cmap=cmap, norm=norm)
        # get the uint8 array
        color_arr = cmap(img.get_array(), bytes=True)
        return color_arr

    def construct(self):
        print(config.frame_rate)
        # initialize the concentration grid
        self._init_c(beta=0.1)
        # get the color array
        color_arr = self.conc_to_color()
        # read the color array as an ImageMobject
        img = ImageMobject(
            color_arr
        )
        # set the size of the ImageMobject
        img.height = 8

        self.play(
            FadeIn(img)
        )
        self.wait()

        # carry out the simulation
        tot_sec = 10
        # determine the second per frame
        spf_simu = tot_sec / self.n_steps  # second per frame required by the simulation
        spf_time = 1 / config.frame_rate  # second per frame allowed by the frame rate
        # determine every n simulation steps to update the image
        every_nsteps = int(np.ceil(spf_time / spf_simu))
        pbar = tqdm(range(1, self.n_steps + 1))
        for _ in pbar:
            pbar.set_description(f"Time step: {_}")
            self._update_c()
            if _ % every_nsteps == 0:
                new_color_arr = self.conc_to_color()
                new_img = ImageMobject(
                    new_color_arr
                )
                new_img.height = 8
                img.become(new_img)
                self.wait(spf_time)

# CahnHilliard2D().render()
