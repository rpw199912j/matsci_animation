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

# set the matplotlib color bar
cmap = mpl.cm.RdYlBu_r
norm = mpl.colors.Normalize(vmin=0, vmax=1)

# define physical variables and constants
R = 8.314  # gas constant
temp = 1000 + 273.15  # temperature [K]
dx, dy = 2.0e-9, 2.0e-9  # spacing of computational grids [m]
La = 39000 - 9. * temp  # Atom interaction constant [J/mol]
ac = 3.0e-14  # gradient coefficient [Jm2/mol]
Da = 1.0e-04 * np.exp(-300000 / R / temp)  # diffusion coefficient of A atom [m2/s]
Db = 2.0e-05 * np.exp(-300000 / R / temp)  # diffusion coefficient of B atom [m2/s]
dt = (dx * dx / Da) * 0.1  # time increment [s]
n_steps = 1200  # total number of time-steps


class CahnHilliard2DTest(Scene):
    def __init__(self):
        super().__init__()
        self.res = 100  # define the sampling resolution for the concentration grid
        self.c0 = 0.5  # set the initial uniform composition
        self.c = np.zeros((self.res, self.res))  # holds the 2D concentration array at t
        self.c_new = np.zeros((self.res, self.res))  # holds the 2D concentration array at t+dt

    def _init_c(self):
        """Initialize the concentration array at t=0"""
        self.c = self.c0 + rng.uniform(size=(self.res, self.res)) * 0.01

    def _update_c(self):
        """Update the concentration array after each time step"""
        # TODO: vectorize this part with matrix operation?
        for j in range(self.res):
            for i in range(self.res):

                ip = i + 1
                im = i - 1
                jp = j + 1
                jm = j - 1
                ipp = i + 2
                imm = i - 2
                jpp = j + 2
                jmm = j - 2

                if ip > self.res - 1:  # periodic boundary condition
                    ip = ip - self.res
                if im < 0:
                    im = im + self.res
                if jp > self.res - 1:
                    jp = jp - self.res
                if jm < 0:
                    jm = jm + self.res
                if ipp > self.res - 1:
                    ipp = ipp - self.res
                if imm < 0:
                    imm = imm + self.res
                if jpp > self.res - 1:
                    jpp = jpp - self.res
                if jmm < 0:
                    jmm = jmm + self.res

                cc = self.c[i, j]  # at (i,j) "centeral point"
                ce = self.c[ip, j]  # at (i+1.j) "eastern point"
                cw = self.c[im, j]  # at (i-1,j) "western point"
                cs = self.c[i, jm]  # at (i,j-1) "southern point"
                cn = self.c[i, jp]  # at (i,j+1) "northern point"
                cse = self.c[ip, jm]  # at (i+1, j-1)
                cne = self.c[ip, jp]
                csw = self.c[im, jm]
                cnw = self.c[im, jp]
                cee = self.c[ipp, j]  # at (i+2, j+1)
                cww = self.c[imm, j]
                css = self.c[i, jmm]
                cnn = self.c[i, jpp]

                mu_chem_c = R * temp * (np.log(cc) - np.log(1.0 - cc)) + La * (
                        1.0 - 2.0 * cc)  # chemical term of the diffusion potential
                mu_chem_w = R * temp * (np.log(cw) - np.log(1.0 - cw)) + La * (1.0 - 2.0 * cw)
                mu_chem_e = R * temp * (np.log(ce) - np.log(1.0 - ce)) + La * (1.0 - 2.0 * ce)
                mu_chem_n = R * temp * (np.log(cn) - np.log(1.0 - cn)) + La * (1.0 - 2.0 * cn)
                mu_chem_s = R * temp * (np.log(cs) - np.log(1.0 - cs)) + La * (1.0 - 2.0 * cs)

                mu_grad_c = -ac * ((ce - 2.0 * cc + cw) / dx / dx + (
                        cn - 2.0 * cc + cs) / dy / dy)  # gradient term of the diffusion potential
                mu_grad_w = -ac * ((cc - 2.0 * cw + cww) / dx / dx + (cnw - 2.0 * cw + csw) / dy / dy)
                mu_grad_e = -ac * ((cee - 2.0 * ce + cc) / dx / dx + (cne - 2.0 * ce + cse) / dy / dy)
                mu_grad_n = -ac * ((cne - 2.0 * cn + cnw) / dx / dx + (cnn - 2.0 * cn + cc) / dy / dy)
                mu_grad_s = -ac * ((cse - 2.0 * cs + csw) / dx / dx + (cc - 2.0 * cs + css) / dy / dy)

                mu_c = mu_chem_c + mu_grad_c  # total diffusion potental
                mu_w = mu_chem_w + mu_grad_w
                mu_e = mu_chem_e + mu_grad_e
                mu_n = mu_chem_n + mu_grad_n
                mu_s = mu_chem_s + mu_grad_s

                nabla_mu = (mu_w - 2.0 * mu_c + mu_e) / dx / dx + (mu_n - 2.0 * mu_c + mu_s) / dy / dy
                dc2dx2 = ((ce - cw) * (mu_e - mu_w)) / (4.0 * dx * dx)
                dc2dy2 = ((cn - cs) * (mu_n - mu_s)) / (4.0 * dy * dy)

                DbDa = Db / Da
                mob = (Da / R / temp) * (cc + DbDa * (1.0 - cc)) * cc * (1.0 - cc)
                dmdc = (Da / R / temp) * ((1.0 - DbDa) * cc * (1.0 - cc) + (cc + DbDa * (1.0 - cc)) * (1.0 - 2.0 * cc))

                dcdt = mob * nabla_mu + dmdc * (dc2dx2 + dc2dy2)  # right-hand side of Cahn-Hilliard equation
                self.c_new[i, j] = self.c[i, j] + dcdt * dt  # update order parameter self.c
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
        # initialize the concentration grid
        self._init_c()
        # get the color array
        color_arr = self.conc_to_color()
        # read the color array as an ImageMobject
        img = ImageMobject(
            color_arr
        )
        # set the size of the ImageMobject
        img.height = 4

        self.play(
            FadeIn(img)
        )
        self.wait()

        # carry out the simulation
        tot_sec = 10
        sec_per_step = tot_sec / n_steps
        pbar = tqdm(range(1, n_steps + 1))
        for _ in pbar:
            pbar.set_description(f"Time step: {_}")
            self._update_c()
            new_color_arr = self.conc_to_color()
            new_img = ImageMobject(
                new_color_arr
            )
            new_img.height = 4
            print()
            # self.wait(0.05)
            img.become(new_img)
            self.wait(0.1)
            # self.play(
            #     img.animate.become(new_img),
            #     run_time=sec_per_step,
            #     rate_func=linear
            # )


class CahnHilliard2D(Scene):
    def __init__(self):
        super().__init__()
        # define the sampling resolution
        self.n_rows = 200
        self.n_cols = 200
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
        # self.sigma = 0.1  # interfacial energy [J/m^2]
        # self.wid = 1e-8  # interface width [m]
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
        img.height = 4

        self.play(
            FadeIn(img)
        )
        self.wait()

        # carry out the simulation
        tot_sec = 10
        sec_per_step = tot_sec / n_steps
        pbar = tqdm(range(1, 5000 + 1))
        for _ in pbar:
            pbar.set_description(f"Time step: {_}")
            self._update_c()
            new_color_arr = self.conc_to_color()
            new_img = ImageMobject(
                new_color_arr
            )
            new_img.height = 4
            print()
            # self.wait(0.05)
            img.become(new_img)
            self.wait(0.1)
            # self.play(
            #     img.animate.become(new_img),
            #     run_time=sec_per_step,
            #     rate_func=linear
            # )

# CahnHilliard2D().render()
