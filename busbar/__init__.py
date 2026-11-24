"""Fallback Project 4 topic (handout menu: "Busbar or resistive heater").

Transient, two-way coupled electro-thermal analysis of a copper EV battery-pack busbar:
    div(sigma(T) grad phi) = 0
    rho c dT/dt = div(k(T) grad T) + sigma(T) |grad phi|^2
with bolted lap-joint contacts (contact resistance, Joule heating in the joint, heat path into
the battery terminals), natural convection h ~ dT^(1/4) and radiation on the exposed surfaces.

This package is self-contained and does not modify the brake-disc code (`brakeopt/`, `optim/`).
It reuses `brakeopt.fem` (hex/quad elements) and `optim.adjoint` (generic backward adjoint sweep),
`optim.functionals` (KS) and `optim.gradcheck` (FD table, Taylor test).
"""
