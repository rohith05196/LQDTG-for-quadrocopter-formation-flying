# LQDTG
**LQDTG stands for Linear Quadratic Discrete time Games for Quadro copter formation flying.**

The Problem is modeled as a game where each vertices of the communication graph represents an agent/player
which minimizes its individual finite horizon-cost function in turn reaching the Nash
Equilibrium. This resulted in a coupled solution due to coupled cost function which cannot
be solved in a distributed manner. To handle this problem, the coupling is relocated to
the dynamics and this resulted in a decoupled cost that allowed solving the problem in a
distributed manner on the edges of the graph in contrary to the previous problem which
was solved by considering the agents on the vertices. Based on the results obtained and to
validate and verify the solution, the problem was further solved using receding horizon
and also, extended to trajectory tracking with formation control

**The Project invoved these steps.**

1. Familiarization with the concept of a Nash equilibrium and its Riccati solution [1]
2. Design and implementation of quadrocopter formation flying as an LQDTG
  • A linearized quadrocopter model as the dynamics
  • A group of quadrocopters maintains a geometric formation shape with the desired distance as the scenario
3. Implementation of a distributed framework for LQDTG
4. Extend the problem to a trajectory tracking as a scenario
5. Design an LQDTG framework for multi-agent collision avoidance

**Brief description of the workflow is depicted by the following flowchart.**
<img width="979" alt="image" src="https://github.com/RohithKamathMijar/LQDTGNew/assets/94147428/edb43c6a-6933-4514-b3da-d227d415920d">

