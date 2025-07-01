# 🌐 CMDPSO

[ICACI2025] Official code for paper **"Coevolutionary Multi-objective Discrete Particle Swarm Optimization for Gateway Placement Optimization Problem"**

## 📋 Project Overview
This repository implements algorithms for solving the Gateway Placement Optimization (GPO) problem. It features a unified discrete PSO framework that supports different velocity update mechanisms to solve binary decision problems.

![Potential gateway positions around the world](fig/potential%20gateway%20position.png)
*The Gateway Placement Optimization problem aims to select the minimum number of gateways from these potential positions to satisfy traffic forwarding requirements and improve service quality.*

## ⚙️ Algorithms

### CMDPSO.m
A multi-objective discrete PSO framework that integrates two velocity update mechanisms:
- **Dual-velocity update**: Uses 2D vectors for position probabilities
- **Sigmoid mapping**: Maps 1D velocity through sigmoid function

Select mechanism by setting `velocity_type = 'bi-velocity'` or `'sigmoid'`.

### BVDPSO_IXP.m
Implementation of Bi-Velocity Discrete PSO that solves the GPO problem as a weighted single-objective problem. 
Unlike CMDPSO, it cannot directly handle multiple objectives but converts them into a single fitness value through weighted aggregation.

## 🗺️ Data & Visualization Files

### Data Processing
- **test_worldmap.m**: Generates gateway data used by optimization algorithms
- **population_print.m**: Processes population data for traffic demand calculation

### Visualization
- **Plot_gatewaymap.m**: Visualizes gateway locations on a world map
- **Plot_population.m**: Displays population distribution data

## 🚀 Usage

### Parameter Settings
```matlab
% General Parameters
NA = 100;        % Archive size limit
Tmax = 100;      % Maximum iterations
DIMENSION = 250; % Variable dimension
PARTICLE_NUM = 50; % Number of particles

% Algorithm Control Parameters
velocity_type = 'bi-velocity'; % Options: 'bi-velocity' or 'sigmoid'
```

### Running the Algorithms
- CMDPSO: Call the `single_time()` or `run_multiple_times()` function
- BVDPSO: Run the script directly

## 📊 Fitness Calculation
The algorithms optimize for a gateway layout problem with objectives:
1. Minimize traffic distribution imbalance
2. Minimize the number of enabled gateways

## 📚 Citation
If you use this code in your research, please cite:
```
@inproceedings{li2025cmdpso,
  title={Coevolutionary Multi-objective Discrete Particle Swarm Optimization for Gateway Placement Optimization Problem},
  author={Li, Jianyu and He, Zhida et al.},
  booktitle={International Conference on Advanced Computational Intelligence},
  year={2025}
}
``` 