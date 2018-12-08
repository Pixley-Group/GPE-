# GPE-
Core function library for numerically solving the Gross Pitaevskii Equation

The overall algorithm cycle is designed to evolve an intial QM wavefunction in imaginary time until it converges to the ground state solution of the Gross Pitaevskii Equation (GPE).

Basic Outline of the Algorithm:
(1) Input the "guess wavefunction" i.e. a constant or something like sin(x), and the potential energy of the system
(2) Input the step size in real space and imaginary time (determines how fine or course grained the approximation is)
(3) The initial wavefunction has is transformed to the momentum basis through the Fast Fourier Transform (FFT).
(4) The Kinetic Energy Operator is applied in real momentum space and imaginary time.
(5) The inverse FFT is applied to transform back to position space and the phase shift function is applied (fftshift).
(6) The Potential Energy Operator is applied to the test function in real space and imaginary time.
(7) The FFT transforms back to momentum space and the KE operator is applied in real momentum and imaginary time.
(8) The output of this cycle is then input to the next cycle (FFT, KE, iFFT, PE, FFT, KE...)
(9) This cycle is iterated over n times to evolve the initial wavefunction (eventually a convergence condition will be added).
(10) The converged solution is output
