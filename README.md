# Particle Filter for Time Series Model

## Overview
This MATLAB function implements a particle filter for estimating state variables in a time-series model. It incorporates different types of Gaussian and non-Gaussian noise models and performs state estimation using a particle-based Bayesian filtering approach.

## Features
- Implements a particle filter with residual resampling.
- Supports multiple noise models:
  - Gaussian noise
  - Composite periodic wave noise
  - Bernoulli-Gaussian impulse noise
  - Rician noise
- Models signal propagation over a uniform planar array.
- Includes state transition and measurement models.

## Function Signature
```matlab
[Xoutput, Xhatoutput, noise] = LEO_3_Particle(std_n, noise_switch)
```

## Inputs
- `std_n`: Standard deviation of the measurement noise.
- `noise_switch`: Integer parameter to select the noise model:
  - `0`: Gaussian noise
  - `1`: Composite periodic wave noise
  - `2`: Bernoulli-Gaussian impulse noise
  - `3`: Rician noise

## Outputs
- `Xoutput`: True state trajectory (theta, phi components).
- `Xhatoutput`: Estimated state trajectory using the particle filter.
- `noise`: Generated measurement noise based on the selected model.

## Usage Example
```matlab
std_n = 0.1;
noise_switch = 0; % Gaussian noise
[Xoutput, Xhatoutput, noise] = LEO_3_Particle(std_n, noise_switch);
```

## Dependencies
- MATLAB Statistics and Machine Learning Toolbox (for `normrnd`, `makedist`, etc.).

## Notes
- The function initializes a uniform planar array and simulates signal reception.
- Particle weights are updated based on likelihood estimates.
- Residual resampling ensures particles are effectively distributed.
- Time step and simulation length are predefined but can be adjusted.

## Author
- 鍾德融

## License
This code is provided for research and educational purposes only.

