# Driven-dissipative free fermions
Codes for numerical simulation of driven-dissiparive free fermions via Lindblad master equations. 

## General setup
By saying "free" fermions, we restrict ourselves to:
- The Hamiltonian contians only quadratic terms in fermion operators.
- The jump operators describing dissipations are linear in fermion operators.

Under these requirements, the time derivatives of two-point correlation functions are determined by two-point functions only, i.e., the evolution of two-point functions is closed, resembling the situation in a closed system of non-interacting fermions. Thus, for simulating the dynamics, we numerically implement the equation of motion of the two-point functions.

This is the foundation of the the "third quantization" method: https://arxiv.org/pdf/0801.1257
