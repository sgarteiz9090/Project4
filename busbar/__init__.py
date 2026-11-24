"""Fallback Project 4 topic (handout menu: "Busbar or resistive heater").

Transient, two-way coupled electro-thermal analysis of a copper EV battery-pack busbar:
    div(sigma(T) grad phi) = 0
    rho c dT/dt = div(k(T) grad T) + sigma(T) |grad phi|^2
with bolted lap-joint contacts (contact resistance, Joule heating in the joint, heat path into
the battery terminals), natural convection h ~ dT^(1/4) and radiation on the exposed surfaces.

The package is self-contained. The generic pieces were first written for our brake-disc study
(kept in `alternate_topic/`) and copied here unchanged: `fem.py` (hex/quad elements), `adjoint.py`
(backward adjoint sweep), `ks.py` (KS aggregation), `gradcheck.py` (FD table, Taylor test, cost)
and `optimize.py` (SLSQP driver with KKT read-out).
"""
