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

1. Familiarization with the concept of a Nash equilibrium and its Riccati solution.
2. Design and implementation of quadrocopter formation flying as an LQDTG.
  • A linearized quadrocopter model as the dynamics.
  • A group of quadrocopters maintains a geometric formation shape with the desired distance as the scenario.
3. Implementation of a distributed framework for LQDTG.
4. Extend the problem to a trajectory tracking as a scenario.
5. Design an LQDTG framework for multi-agent collision avoidance.

**Brief description of the workflow is depicted by the following flowchart.**
<img width="979" alt="image" src="https://github.com/RohithKamathMijar/LQDTGNew/assets/94147428/edb43c6a-6933-4514-b3da-d227d415920d">

Two Scenarios where conidered 
1. Three agest with the following communication pattern
<img width="684" alt="Screenshot 2024-03-22 at 10 20 48" src="https://github.com/RohithKamathMijar/LQDTG_Rohith/assets/94147428/07c68138-ee9f-4b53-8e72-66037561ff1a">


2. Four agents with the comminication pattern  shown in the picture.
<img width="684" alt="Screenshot 2024-03-22 at 10 21 31" src="https://github.com/RohithKamathMijar/LQDTG_Rohith/assets/94147428/c729b127-7135-40b4-8b8f-61f2328d9a02">


**Results Obtained**
1. Formation control of 3 agents starting at random initial position and following helical trajectory
<img width="726" alt="Screenshot 2024-03-22 at 10 22 13" src="https://github.com/RohithKamathMijar/LQDTG_Rohith/assets/94147428/b9ec46a8-22e8-483e-8bed-83bd8f15397c">

2.Formation control of 3 agents starting at random initial position and following helical trajectory
<img width="726" alt="Screenshot 2024-03-22 at 10 22 42" src="https://github.com/RohithKamathMijar/LQDTG_Rohith/assets/94147428/c51198fe-8114-4e51-8f49-b0d78a5d1069">
<img width="726" alt="Screenshot 2024-03-22 at 10 23 01" src="https://github.com/RohithKamathMijar/LQDTG_Rohith/assets/94147428/2c0181aa-f99f-4646-8a17-d0780df07b18">
