"""Quadrotor formation control trajectories and visualization."""

from dataclasses import dataclass
from typing import Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class SimulationParameters:
    """Container for common simulation parameters."""

    n: int = 3
    ns: int = 12
    ni: int = 4
    tf: float = 20.0
    dt: float = 0.1
    eta: float = 4.0
    step: int = 100


class QuadrotorFormationControl:
    """Generate reference trajectories and plot formation tracking."""

    def __init__(self, section: int = 1) -> None:
        self.section = section
        self.params = SimulationParameters()
        self.agent_count = 3 if self.section in (1, 3, 5, 6) else 4
        self.F = None

    def setup_matrices(self) -> None:
        """Create discrete time system matrix used by the controller."""

        n = self.params.n
        N = self.agent_count
        A11i = np.zeros((2 * n, 2 * n))
        A12i = np.eye(2 * n)
        A21i = np.block(
            [
                [np.zeros((n, n)), np.array([[0, -9.81, 0], [9.81, 0, 0], [0, 0, 0]])],
                [np.zeros((n, 2 * n))],
            ]
        )
        A22i = np.zeros((2 * n, 2 * n))

        A = np.block(
            [
                [np.kron(np.eye(N), A11i), np.zeros((2 * n * N, 1)), np.kron(np.eye(N), A12i)],
                [np.zeros((1, 2 * n * N)), 0, np.zeros((1, 2 * n * N))],
                [np.kron(np.eye(N), A21i), np.zeros((2 * n * N, 1)), np.kron(np.eye(N), A22i)],
            ]
        )

        dt = self.params.dt
        eye = np.eye(self.params.ns * N + 1)
        self.F = eye + dt * A + (dt**2 / 2) * (A @ A) + (dt**3 / 6) * (A @ A @ A)

    def _base_trajectory(self, t: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return reference (x, y, z) for the leader trajectory."""

        x = np.sin(t)
        y = np.cos(t) if self.section in (1, 2) else np.sin(2 * t)
        z = t * 0.1
        return x, y, z

    def _agent_phase_offsets(self) -> Iterable[float]:
        """Phases evenly distributed for all agents."""

        for idx in range(self.agent_count):
            yield idx * 2 * np.pi / self.agent_count

    def _generate_agent_trajectories(self, t: np.ndarray) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Create trajectories for each agent with phase shifts."""

        _, _, z_lead = self._base_trajectory(t)
        trajectories: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []

        for phase in self._agent_phase_offsets():
            x = np.sin(t + phase)
            y = np.cos(t + phase) if self.section in (1, 2) else np.sin(2 * (t + phase))
            trajectories.append((x, y, z_lead))

        return trajectories

    def _plot_trajectories(self, trajectories: List[Tuple[np.ndarray, np.ndarray, np.ndarray]]) -> None:
        """Plot trajectories for all agents with start and end markers."""

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        lines = []

        for idx, (x, y, z) in enumerate(trajectories):
            (line,) = ax.plot(x, y, z)
            lines.append(line)
            ax.scatter(x[0], y[0], z[0], c="black", marker="o", s=50)
            ax.scatter(x[-1], y[-1], z[-1], c="red", marker="^", s=50)

        ax.legend(lines, [f"Agent {idx + 1}" for idx in range(self.agent_count)])
        ax.set_title(f"Section {self.section}: Trajectory Tracking for {self.agent_count} agents")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")

        plt.show()

    def run_section(self) -> None:
        """Compute matrices, generate trajectories, and plot results."""

        self.setup_matrices()
        time_vector = np.linspace(0, self.params.tf, int(self.params.tf / self.params.dt))
        trajectories = self._generate_agent_trajectories(time_vector)
        self._plot_trajectories(trajectories)


def run_all_sections(sections: Iterable[int]) -> None:
    """Run the visualization for the provided section numbers."""

    for section in sections:
        controller = QuadrotorFormationControl(section=section)
        controller.run_section()


if __name__ == "__main__":
    run_all_sections(range(1, 7))
