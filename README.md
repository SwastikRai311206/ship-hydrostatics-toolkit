# Ship Hydrostatics Toolkit

A Python-based tool for calculating ship hydrostatics and stability
parameters, built using Simpson's 1/3rd Rule and standard naval
architecture formulas.


## What This Project Does
This project calculates key hydrostatic and stability values for a
ship's hull, including:

- Waterplane Area (WPA)
- Longitudinal Center of Flotation (LCF)
- Moment of Inertia
- Tonnes per centimeter (TPC)
- KB, BM, KM, and GM (metacentric height)

## How It Works
Simpson's 1/3rd Rule is used to numerically integrate a table of
half-ordinates (half-breadths) taken at equally spaced stations
along the ship's waterline. This method works for any hull shape,
unlike simple box-hull formulas.

The calculated waterplane area and moment of inertia are then used
to compute the ship's stability parameters (KB, BM, KM, GM).

## Files
- `simpsons_13_rule.py` – Core functions for Simpson's Rule (waterplane
  area, moment of inertia, LCF).
- `demo_simpsons_hydrostatics.py` – Applies Simpson's Rule to a
  sample hull and calculates full hydrostatics, with a plot of the
  waterplane shape.

## Assumptions
- Displaced volume is calculated as `V = A_w × T`, assuming a
  wall-sided hull (vertical sides down to the keel).
- Ship stations are assumed to be equally spaced.

## Validation
The Simpson's Rule calculations were checked against a textbook
example and produced matching results, confirming the method is
implemented correctly.

## Tools Used
Python, NumPy, Matplotlib

## Author
[Swastik Rai] – Naval Architecture and Ocean Engineering,
Indian Maritime University, Visakhapatnam

