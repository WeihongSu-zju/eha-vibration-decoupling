# eha-vibration-decoupling
MATLAB code for Decoupling of Multi-Source Excitation and Complex Vibration Analysis

This repository provides the MATLAB implementation used in the simulation study of the paper:
**"Decoupling of Multi-Source Excitation and Complex Vibration Analysis in Underwater Electro-Hydraulic Actuators"**

The code demonstrates the key signal-processing steps for analyzing **time-varying multi-source vibration signals**, including:
- instantaneous frequency extraction
- instantaneous amplitude estimation
The simulation validates the robustness of the proposed method under different noise conditions.

## Method Overview
The implemented framework corresponds to **Section 3 of the paper** and consists of two main stages.
### 1. Instantaneous Frequency Extraction
- Time-frequency analysis using **Short-Time Fourier Transform (STFT)**
- Frequency tracking using **dynamic programming**
- Robust estimation of instantaneous frequency from the spectrogram
### 2. Amplitude Estimation
- Harmonic component extraction using **Vold–Kalman filtering**
- Reconstruction of time-varying signal components
- Quantitative evaluation using correlation and MSE
  
Both **simulation signals** and **experimental vibration data** are provided.
## Experimental Dataset
The dataset contains vibration signals of the underwater EHA under different outlet pressures during a **0–2000 RPM run-up test**.
Pressure conditions included: 0 MPa, 2 MPa, 4 MPa, 6 MPa, 8 MPa, 10 MPa, 12 MPa, 14 MPa, 16 MPa, 18 MPa, 20 MPa.

## License
This code is provided for **academic research and educational purposes only**.
