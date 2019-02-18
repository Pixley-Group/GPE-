#Library of functions for solving the 2D GPE.

#input parameters--------------------------

#size of step in imaginary time
delt_t = -(10^-3)*im

#size of step in real space
delt_x = 10^-3
delt_y = 10^-3

#range over real space (-x_end:x_end)
x_end = 1
y_end = 1

#The number of steps in real space
N_x = 10^3
N_y = 10^3
x_array = collect(0:delt_x:x_end)
y_array = collect(0:delt_y:y_end)

#physical constants
k = 1
m = 1
#strength of particle-particle interactions (coupling parameter)
g = 5*10^2
#number of particles
N = 10^3
#---------------------------------------------------

# Write the guess wavefunction here
function psi_guess(x,y)
    return pi
end

#write the potential energy here
function pot(x,y)
    pot = ((k/2)*((x^2)+(y^2)))
end

#core Functions for the algortihm___________________________________

using LinearAlgebra, Statistics, Compat
#generates the array that is the discretization of the guess function in real space
function psi_guess_array()
    psi_guess_array = zeros(ComplexF64, N_x, N_y)

    for x in 1:N_x
        for y in 1:N_y
            psi_guess_array[x,y] = psi_guess(x,y)
        end
    end
    return psi_guess_array
end

#normalizes the wavefunction array
function normalizer(array)
    a = array.* conj(array)
    s = sqrt(sum(a))
    return (array./s).*sqrt(N)
end

#Generates the potential energy array
function pot_array(psi)
    pot_array = zeros(ComplexF64, N_x, N_y)

    for n in 1:N_x
        for m in 1:N_y
            pot_array[n,m] = pot(x_array[n], y_array[m])
        end
    end

    pot_array = pot_array.+g*(conj(psi).*psi)
    return pot_array

end

#calculates the x momentum for the fft minted momentum eignestate
function p_x(n)
    p_x = ((2*pi)*n)/N_x
end

#calculates the y momentum for the fft minted momentum eignestate
function p_y(m)
    p_y = ((2*pi)*m)/N_y
end

#calculates the total kinetic energy for each fft minted momentum eigenstate
function kin(p_x, p_y)
    return ((p_x^2)+(p_y^2))/2m
end

#Generates the Kinetic Energy array
function kin_array()
    kin_array = zeros(ComplexF64, N_x, N_y)

    for n in 1:N_x
        for m in 1:N_y
            kin_array[n,m] = kin(p_x(n), p_y(m))
        end
    end
    return kin_array
end

#exponentiates the pot energy operator elementwise
function e_V(psi)
    return exp.(-pot_array(psi).*((delt_t)*im))
end

#exponentiates the kin energy operator elementwise
function e_T()
    return e_T = exp.(-kin_array().*(delt_t/2)*im)
end

#evolves psi delt_t in time with the KE operator
function time_step_T(array)
    array.*e_T()
end

#evolves psi delt_t in time with the PE operator
function time_step_V(array)
    array.*e_V(array)
end

#evolves the guess function array one step in imaginary time
using FFTW
function time_evolve_step_im(array)
    psi_k_1 = fftshift(fft(array))
    psi_k_T_1 = time_step_T(psi_k_1)
    psi_x_T_1 = ifft(psi_k_T_1)
    psi_x_V = time_step_V(psi_x_T_1)
    psi_k = fftshift(fft(psi_x_V))
    psi_k_T = time_step_T(psi_k)
    psi_x_T = ifft(psi_k_T)
    return normalizer(psi_x_T)
end

#evolves the guess function array t steps in imaginary time
function time_evolve_im(array, t)
    evolved_array = reduce((x, y) -> time_evolve_step_im(x), 1:t, init=array)
    return normalizer(evolved_array)
end

#Finds the energy of psi
function Energy(array)
    ffts = fftshift(fft(array))
    kin = ifft(ffts.*kin_array())
    pot = psi.*pot_array(psi)
    energy = (kin.+pot).*conj(psi)
    return sum(energy)
end





#Plotters____________________________________________________
# using Plots
# #plots guess psi
# psi = psi_guess_array()
# time_step = time_evolve_step(psi)
# t = time_evolve(psi, 100)
# t_2 = conj.(t) .* t
#
# #plots a static plot of the real modulus of the evolved psi
# plot(x_array, real(t_2), title = "Psi Evolved")

#uncomment to create animations of the evolution of the guess function psi
# anim = @animate for i=1:800
#     t = time_evolve(psi, i)
#     plot(x_array, real(conj.(t) .* t), title = "Psi Evolved")
# end
# gif(anim, "C:/Users/Alucard/Desktop/julia/gifs/2D_test_fps10.gif", fps = 10)
