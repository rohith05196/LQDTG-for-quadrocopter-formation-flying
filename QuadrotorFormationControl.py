import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm

class QuadrotorFormationControl:
    def __init__(self, section=1):
        self.n = 3
        self.N = 3 if section in [1, 3, 5, 6] else 4
        self.M = 2 if section in [1, 3, 5, 6] else 4
        self.ns = 12
        self.ni = 4
        self.tf = 20
        self.dt = 0.1
        self.eta = 4
        self.section = section
        self.step = 100
        self.setup_matrices()

    def setup_matrices(self):
        n, N = self.n, self.N
        A11i = np.zeros((2 * n, 2 * n))
        A12i = np.eye(2 * n)
        A21i = np.block([
            [np.zeros((n, n)), np.array([[0, -9.81, 0], [9.81, 0, 0], [0, 0, 0]])],
            [np.zeros((n, 2 * n))]
        ])
        A22i = np.zeros((2 * n, 2 * n))

        self.A = np.block([
            [np.kron(np.eye(N), A11i), np.zeros((2 * n * N, 1)), np.kron(np.eye(N), A12i)],
            [np.zeros((1, 2 * n * N)), 0, np.zeros((1, 2 * n * N))],
            [np.kron(np.eye(N), A21i), np.zeros((2 * n * N, 1)), np.kron(np.eye(N), A22i)]
        ])

        self.F = (np.eye(self.ns * N + 1) + self.dt * self.A + (self.dt ** 2 / 2) * (self.A @ self.A) +
                   (self.dt ** 3 / 6) * (self.A @ self.A @ self.A))

        self.m = 0.64
        B1i = np.vstack([
            np.zeros((self.ni, self.ni)),
            np.array([[0, 0, 1, 0], [0, 0, 0, 1 / self.m]])
        ])
        B3i = np.vstack([
            np.zeros((n, self.ni)),
            np.hstack([np.eye(n), np.zeros((n, 1))])
        ])
        self.B = np.vstack([
            np.kron(np.eye(N), B1i),
            np.zeros((1, N * self.ni)),
            np.kron(np.eye(N), B3i)
        ])

        self.Bb = np.zeros((self.B.shape[0], self.ni, N))
        for i in range(N):
            self.Bb[:, :, i] = self.B[:, i * self.ni: (i + 1) * self.ni]

        self.Gb = np.zeros_like(self.Bb)
        for i in range(N):
            self.Gb[:, :, i] = (self.dt * self.Bb[:, :, i] + (self.dt ** 2) / 2 * self.A @ self.Bb[:, :, i] +
                                (self.dt ** 3) / 6 * (self.A @ self.A) @ (self.Bb[:, :, i] ** 2) +
                                (self.dt ** 4) / 24 * (self.A @ self.A @ self.A) @ (self.Bb[:, :, i] ** 3))

    def run_section(self):
        if self.section == 1:
            self.run_section_1()

    def run_section_1(self):
        print("Running Section 1: Helical Trajectory Tracking for 3 agents")
      
    def plot_trajectories(self, xhat):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for i in range(self.N):
            idx = i * 6
            ax.plot(xhat[idx, :], xhat[idx + 1, :], xhat[idx + 2, :], label=f'Agent {i+1}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.legend()
        ax.set_title('3D Formation Trajectories')
        plt.show()
      
if __name__ == "__main__":
    formation_control = QuadrotorFormationControl(section=1)
    formation_control.run_section()
