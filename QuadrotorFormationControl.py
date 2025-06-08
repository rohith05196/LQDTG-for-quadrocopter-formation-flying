import numpy as np
import matplotlib.pyplot as plt

class QuadrotorFormationControl:
    def __init__(self, section=1):
        self.section = section
        self.initialize_parameters()

    def initialize_parameters(self):
        self.n = 3
        self.N = 3 if self.section in [1, 3, 5, 6] else 4
        self.M = 2 if self.section in [1, 3, 5, 6] else 4
        self.ns = 12
        self.ni = 4
        self.tf = 20
        self.dt = 0.1
        self.eta = 4
        self.step = 100

    def setup_matrices(self):
        self.A11i = np.zeros((2 * self.n, 2 * self.n))
        self.A12i = np.eye(2 * self.n)
        self.A21i = np.block([
            [np.zeros((self.n, self.n)), np.array([[0, -9.81, 0], [9.81, 0, 0], [0, 0, 0]])],
            [np.zeros((self.n, 2 * self.n))]
        ])
        self.A22i = np.zeros((2 * self.n, 2 * self.n))
        self.A = np.block([
            [np.kron(np.eye(self.N), self.A11i), np.zeros((2 * self.n * self.N, 1)), np.kron(np.eye(self.N), self.A12i)],
            [np.zeros((1, 2 * self.n * self.N)), 0, np.zeros((1, 2 * self.n * self.N))],
            [np.kron(np.eye(self.N), self.A21i), np.zeros((2 * self.n * self.N, 1)), np.kron(np.eye(self.N), self.A22i)]
        ])
        self.F = np.eye(self.ns * self.N + 1) + self.dt * self.A + (self.dt ** 2 / 2) * (self.A @ self.A) + (self.dt ** 3 / 6) * (self.A @ self.A @ self.A)

    def run_section(self):
        self.setup_matrices()
        t = np.linspace(0, self.tf, int(self.tf / self.dt))
        if self.section in [1, 2]:
            # Helical Trajectory
            x1 = np.sin(t)
            y1 = np.cos(t)
            z1 = t * 0.1
        else:
            # Infinity Trajectory
            x1 = np.sin(t)
            y1 = np.sin(2 * t)
            z1 = t * 0.1

        if self.section in [1, 3, 5, 6]:
            agents = 3
        else:
            agents = 4

        trajectories = []
        for i in range(agents):
            phase = i * 2 * np.pi / agents
            trajectories.append([
                np.sin(t + phase),
                np.cos(t + phase) if self.section in [1, 2] else np.sin(2 * (t + phase)),
                z1
            ])

        # Plot the 3D trajectories for all the cases 
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for i, traj in enumerate(trajectories):
            ax.plot(traj[0], traj[1], traj[2], label=f'Agent {i+1}')
        ax.legend()
        ax.set_title(f'Section {self.section}: Trajectory Tracking for {agents} agents')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()

if __name__ == '__main__':
    for section in range(1, 7):
        controller = QuadrotorFormationControl(section=section)
        controller.run_section()
