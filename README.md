# ROSIF-HEX

### Shell-and-Tube Heat Exchanger Simulation & Optimization Software

ROSIF-HEX is a desktop engineering software tool for the **simulation, thermal-hydraulic analysis, and optimization of shell-and-tube heat exchangers**.

The software was developed from chemical engineering heat-transfer and fluid-flow principles and extended into a general-purpose computational tool for engineers, students, researchers, operators, and process design applications.

---

## Overview

Shell-and-tube heat exchanger design requires the simultaneous consideration of heat-transfer performance, fluid-flow behavior, geometry, and pumping requirements.

ROSIF-HEX brings these calculations into a single interactive environment, allowing users to define heat exchanger parameters, perform thermal and hydraulic calculations, evaluate design alternatives, and investigate the relationship between **heat-transfer rate and pumping power**.

The software is not designed around a single academic case study. The underlying engineering principles are implemented as a **general-purpose modelling and analysis framework** that can be adapted to different heat exchanger configurations and operating conditions.

---

## Key Capabilities

* Shell-and-tube heat exchanger modelling
* Thermal performance calculations
* Hydraulic and pressure-drop calculations
* Pumping-power estimation
* Heat-transfer coefficient calculations
* Geometry-based design calculations
* Operating-condition analysis
* Temperature iteration
* Design-variable evaluation
* Performance comparison
* Multi-objective design analysis
* Q vs. pumping-power analysis
* KTHA-I optimization framework
* Interactive engineering interface
* Graphical presentation of simulation results

---

## KTHA-I

ROSIF-HEX incorporates the **Khokmah Thermal Hydraulic Algorithm (KTHA-I)** as an optimization framework for evaluating the relationship between heat-transfer performance and hydraulic energy requirements.

The framework was developed from the observation that increasing fluid velocity can improve heat transfer while simultaneously increasing pressure drop and pumping requirements.

Khokmah Thermal-Hydraulic Algorithm (KTHA)
1. Thermal-Hydraulic Leverage

Continuous form:
L₍Q,P₎ = d ln(Q) / d ln(Pₜₒₜₐₗ)
Discrete form:
Lᵢ = ln(Qᵢ / Qᵢ₋₁) / ln(Pᵢ / Pᵢ₋₁)

Where:

Q = heat-transfer rate
Pₜₒₜₐₗ = total pumping power
L₍Q,P₎ = thermal-hydraulic leverage
2. Normalized KTHA-I

Iᵦ = (Q / Q₀) / (Pₜₒₜₐₗ / P₀)ᵝ

Where:

Iᵦ = normalized KTHA-I performance index
Q = heat-transfer rate
Q₀ = reference heat-transfer rate
Pₜₒₜₐₗ = total pumping power
P₀ = reference total pumping power
β = hydraulic penalty exponent
3. Optimum Condition

At a smooth interior optimum:

d ln(Iᵦ) / d ln(Pₜₒₜₐₗ) = L₍Q,P₎ − β = 0

Therefore:

L₍Q,P₎ = β

This is the fundamental optimum condition of the KTHA-I framework.

4. Physics-Informed β = 1/3

Under simplified turbulent, fixed-geometry scaling:

ΔP ∝ v²

Pumping power is:

Pₜₒₜₐₗ = ṼΔP

Since:

Ṽ ∝ v

then:

Pₜₒₜₐₗ ∝ v³

Therefore:

Pₜₒₜₐₗ¹ᐟ³ ∝ v

This motivates:

β = 1/3

as a physics-informed candidate, not a universal law.

5. KTHOP

The Khokmah Thermal-Hydraulic Optimum Principle (KTHOP) states that the useful operating region is the region in which additional hydraulic expenditure continues to provide sufficient additional thermal benefit.

The KTHA peak is the maximum of the selected KTHA performance criterion.

So the framework can be summarized as:

Hydraulic expenditure → Thermal benefit → Leverage → Optimum

And the key relationship is:

L₍Q,P₎ = β

with β = 1/3 used as the physics-informed candidate under the stated scaling assumptions.
---

## Engineering Model

The software integrates the major components of shell-and-tube heat exchanger analysis, including:

### Thermal Analysis

* Heat-transfer rate
* Log mean temperature difference
* Overall heat-transfer coefficient
* Convective heat-transfer coefficients
* Heat-transfer area
* Temperature relationships

### Hydraulic Analysis

* Reynolds number
* Flow velocity
* Friction factor
* Pressure drop
* Pumping power
* Shell-side and tube-side hydraulic behaviour

### Geometry

The model can account for key shell-and-tube heat exchanger geometric parameters such as:

* Tube outside diameter
* Tube inside diameter
* Tube length
* Number of tubes
* Tube pitch
* Tube arrangement
* Shell diameter
* Baffle spacing
* Number of baffles
* Clearance
* Flow area

The specific variables available depend on the current implementation of the software.

---

## Optimization Concept

The design of a heat exchanger is inherently a trade-off.

Increasing flow velocity can increase:

**Heat transfer**

but can also increase:

**Pressure drop → Pumping power → Operating energy requirement**

ROSIF-HEX therefore provides a framework for examining the design space rather than simply searching for the highest possible heat-transfer rate.

A conceptual representation is:

```text
                 Higher velocity
                       │
                       ▼
                ┌─────────────┐
                │ Heat Transfer│
                │      ↑       │
                └──────┬──────┘
                       │
                       │
                       ▼
                ┌─────────────┐
                │ Pressure Drop│
                │      ↑       │
                └──────┬──────┘
                       │
                       ▼
                ┌─────────────┐
                │Pumping Power │
                │      ↑       │
                └─────────────┘
```

The objective is therefore not simply:

> Maximize heat transfer.

Instead, the software investigates how thermal performance can be achieved while accounting for the hydraulic energy required to operate the exchanger.

---

## Software Architecture

ROSIF-HEX is structured as a desktop engineering application with separate components for:

```text
User Input
    │
    ▼
Geometry Definition
    │
    ▼
Thermal Calculations
    │
    ▼
Hydraulic Calculations
    │
    ▼
Performance Evaluation
    │
    ▼
Optimization / KTHA-I
    │
    ▼
Results & Visualization
```

This separation allows the computational models to be extended without requiring the user interface to be redesigned.

---

## Technology Stack

ROSIF-HEX is built primarily with:

* **Python**
* **PyQt5**
* Numerical and engineering computation libraries
* Data-processing tools
* Graphical plotting/visualization libraries

The exact dependencies are listed in `requirements.txt`.

---

## Interface

### Main Interface

![ROSIF-HEX Interface](screenshots/dashboard.png)

The main interface provides access to the heat exchanger inputs, simulation controls, calculations, and results.

### Geometry & Design Inputs

![Geometry Inputs](screenshots/geometry.png)

Users can define the relevant heat exchanger geometry and operating parameters.

### Simulation Results

![Simulation Results](screenshots/simulation.png)

The simulation environment presents calculated thermal and hydraulic performance.

### Optimization

![Optimization Results](screenshots/optimization.png)

The optimization module evaluates different operating/design conditions and their corresponding thermal and hydraulic performance.

### Results

![Results](screenshots/results.png)

Results can be examined through numerical outputs and graphical representations.

---

## General-Purpose Design

ROSIF-HEX was developed using the principles and computational requirements of shell-and-tube heat exchanger analysis, but the software is **not locked to a single academic project, exchanger, or set of experimental values**.

The goal is to provide a reusable engineering tool that can support:

* Chemical engineers
* Process engineers
* Mechanical engineers
* Plant engineers
* Researchers
* Students
* Operators
* Heat exchanger designers
* Engineering educators

The software architecture is intended to allow future expansion into broader process-plant and thermal-system applications.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/IfiokUdoma/rosif-hex.git
```

Move into the project directory:

```bash
cd rosif-hex
```

Create a virtual environment:

```bash
python -m venv venv
```

Activate the environment.

### Windows

```bash
venv\Scripts\activate
```

### macOS / Linux

```bash
source venv/bin/activate
```

Install the required dependencies:

```bash
pip install -r requirements.txt
```

Run the application:

```bash
python main.py
```

> Replace `main.py` with the actual entry-point file used by your application.

---

## Project Structure

A simplified representation of the project is:

```text
rosif-hex/
│
├── src/                  # Application source code
├── screenshots/          # Interface and result screenshots
├── assets/               # Images and project assets
├── docs/                 # Technical documentation
├── tests/                # Tests
├── requirements.txt      # Python dependencies
├── README.md             # Project documentation
├── LICENSE               # Project license
└── .gitignore            # Ignored files
```

Update this section to match the actual repository structure.

---

## Development Status

🟢 **Active Development**

ROSIF-HEX is being developed as a general-purpose thermal-hydraulic engineering software tool.

Current development areas include:

* Model validation
* Expanded geometry support
* Improved temperature iteration
* Enhanced optimization workflows
* Additional visualization
* Engineering validation cases
* Improved user guidance
* Expanded documentation

---

## Roadmap

Future development may include:

* [ ] Expanded shell-side calculation methods
* [ ] Additional heat exchanger configurations
* [ ] Automated design comparison
* [ ] Advanced optimization algorithms
* [ ] Improved validation and benchmarking
* [ ] Exportable engineering reports
* [ ] Enhanced iteration and convergence controls
* [ ] More comprehensive design constraints
* [ ] Process-plant integration
* [ ] Additional equipment models

---

## Engineering Context

ROSIF-HEX originated from a chemical engineering heat exchanger modelling and optimization study and was subsequently developed into a broader software concept.

The project explores how computational engineering tools can move beyond static calculations toward interactive systems that help engineers evaluate **design trade-offs, operating conditions, and energy requirements**.

The KTHA-I framework represents one of the project's original contributions to this direction.

---

## Author

**Ifiok Essienubong Udoma**

Chemical Engineer • Product Builder • Software Developer

Interested in building software, engineering systems, AI applications, energy technologies, and products that solve real-world problems.

---

## License

This project is licensed under the **MIT License**.

See the `LICENSE` file for details.

---

## Disclaimer

ROSIF-HEX is an engineering software project intended for educational, research, development, and engineering analysis purposes.

Results should be independently validated against appropriate engineering standards, correlations, experimental data, and operating constraints before being used for safety-critical or commercial process-design decisions.
# ROSIF-HEX
Desktop shell-and-tube heat exchanger simulation and optimization software for thermal-hydraulic analysis and engineering design.
