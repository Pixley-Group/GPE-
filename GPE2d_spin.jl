#Library of functions for solving the 2D GPE with SOC Hamiltonian.

#input parameters--------------------------

#size of step in real time

fib(n) = n < 2 ? n : fib(n-1) + fib(n-2)

#range over real space (-x_end:x_end)
#V HO "spring constant"
# const m = 1
# const omega = 0
#strength of particle-particle interactions (V coupling constant)
# const g = 0
# Energy scale for Pot QP
# W = .2
#phase for Pot QP

#momentum kinetic energy strength
# const t = 0
#number of particles
# const N = 1
#---------------------------------------------------

#----Enter the starting parameters here
# Write the guess wavefunction here
function psi_guess(x,y)
    return exp(-((x-45)^2) - ((y-45)^2))
end

#write the potential energy here
function pot_H_O(x,y,m,omega)
    pot = m*(omega^2)*(x^2)+(y^2)
end

#Write the Quasi-Periodic Potential Energy here
function pot_QP(x, y, W, phi_x, phi_y, n, m)
    return W*(cos((((2*pi)*(m/n))*(x))+phi_x) + cos((((2*pi)*(m/n))*(y))+phi_y))
end


#core Functions for the algortihm___________________________________

using LinearAlgebra, Statistics, Compat
#generates the array that is the discretization of the guess function in real space
function psi_guess_array(psi_guess_array, n)

    for x in 1:n
        for y in 1:n
            psi_guess_array[x,y,1] = psi_guess(x,y)
            psi_guess_array[x,y,2] = 0
        end
    end
    return normalizer(psi_guess_array)
end

#normalizes the wavefunction array
#array[:,:,1].* conj(array[:,:,1])) + (array[:,:,2].*conj(array[:,:,2])
function normalizer(array)
    s = sqrt(dot(array, array))
    return (array/s)
end

#-----------------POTENTIAL ENERGY

#Generates the potential energy array
function pot_array_H_O(psi, pot_array)

    for i in 1:n
        for j in 1:n
            pot_array[i,j] = pot_H_O(x_array[i], y_array[j])
        end
    end

    psi = psi[:,:,1]+psi[:,:,2]
    psi = reshape(psi, n, n)

    pot_array = pot_array+g*(conj(psi).*psi)
    return pot_array

end



function pot_array_QP(pot_array, W, phi_x, phi_y, n, m)

    for i in 1:n
        for j in 1:n
            pot_array[i,j] = pot_QP(i, j, W, phi_x, phi_y, n, m)
        end
    end

    return pot_array

end

pot_matrix_QP = pot_array_QP(zeros(ComplexF64, 89, 89), 0.1, 0, 0, 89, 34)

#---------------KINETIC ENERGY

#calculates the x momentum for the fft minted momentum eignestate
function p_x(x,n)
    return ((2*pi)*(x-1))/n
end

#calculates the y momentum for the fft minted momentum eignestate
function p_y(y,n)
    return ((2*pi)*(y-1))/n
end

#calculates the total kinetic energy for each fft minted momentum eigenstate
function kin_mom(p_x, p_y)
    return ((p_x)^2+(p_y)^2)
end


#Generates the momentum Kinetic Energy array
function kin_mom_array(kin_array, n)

    for i in 1:n
        for j in 1:n
            kin_array[i,j] = kin_mom(p_x(i,n), p_y(j,n))
        end
    end
    return kin_array
end

#generates the spin orbital coupling matrix


function kin_spin_matrix(spin_couple_matrix, n, del_t)
    A(x,y,n) = sin(p_x(x,n)) - im*sin(p_y(y,n))
    for x in 1:n
        for y in 1:n
            spin_couple_matrix[x,y,1,1] = cos(abs(A(x,y,n)) * del_t/2)
            spin_couple_matrix[x,y,1,2] = -im*exp(im * angle(A(x,y,n))) * sin(abs(A(x,y,n)) * del_t/2)
            spin_couple_matrix[x,y,2,1] = -im*exp(-im * angle(A(x,y,n))) * sin(abs(A(x,y,n)) * del_t/2)
            spin_couple_matrix[x,y,2,2] = cos(abs(A(x,y,n)) * del_t/2)
        end
    end
    return spin_couple_matrix
end

spin_matrix = kin_spin_matrix(zeros(ComplexF64, 89, 89, 2, 2), 89, 40^-1)

#exponentiates the harmonic oscilator pot energy operator elementwise
function e_V(psi)
    return exp.(((-m*(omega)^2)/2)*pot_array_H_O(psi).*((-del_t)*im))
end

#exponentiates the momentum kin energy operator elementwise
function e_T()
    return e_T = exp.(-t*kin_mom_array()*(del_t/2)*im)
end

#evolves psi delt_t in time with the KE operator
#evaluate this one time and store in memory

function time_step_T(array,n, gen_array)
    for x in 1:n
        for y in 1:n
            gen_array[x,y,1] = spin_matrix[x,y,1,1]*array[x,y,1] + spin_matrix[x,y,1,2]*array[x,y,2]
            gen_array[x,y,2] = spin_matrix[x,y,2,1]*array[x,y,1] + spin_matrix[x,y,2,2]*array[x,y,2]
        end
    end
    return gen_array
end

#evolves psi delt_t in time with the PE operator
function time_step_V(array, del_t)
    return array.*(exp.(pot_matrix_QP*(-im*del_t)))
end


using FFTW
using AbstractFFTs

function init_FFT(A, region)
    return plan_fft(A,region; flags=FFTW.PATIENT, timelimit=Inf)
end

function init_IFFT(A, region)
    return plan_ifft(A,region; flags=FFTW.PATIENT, timelimit=Inf)
end



function time_evolve_step(array, del_t, F_T, I_F_T)
    return I_F_T*(time_step_T(F_T*time_step_V(I_F_T*time_step_T(F_T*array,89,zeros(ComplexF64, 89, 89, 2)), del_t),89,zeros(ComplexF64, 89, 89, 2)))
end

#evolves the guess function array t steps in imaginary time
function time_evolve(array, t, F_T, I_F_T, del_t)
    evolved_array = reduce((x, y) -> time_evolve_step(x, del_t, F_T, I_F_T), 1:t, init=array)
    return evolved_array
end

using ProgressMeter

function time_evolve_fast(array, t, del_t, F_T, I_F_T)

    n_array = [array]

    @progress for i in 2:t
        push!(n_array, time_evolve_step(n_array[i-1], del_t, F_T, I_F_T))
    end
    return n_array
end

function T(x,y,n)
    return [[0 sin(p_x(x,n)) - (im* sin(p_y(y,n)))]; [sin(p_x(x,n)) + im * sin(p_y(y,n)) 0]]
end

function T_step(array, gen_array,n)
    for x in 1:n
        for y in 1:n
            gen_array[x,y,:] = T(x,y,n)*array[x,y,:]
        end
    end
    return gen_array
end

#Finds the energy of psi
function Energy(array, F_T, I_F_T)
    ffts = F_T*array
    kin = dot(array, I_F_T*T_step(ffts, zeros(ComplexF64, 89, 89, 2), 89))
    pot = dot(array, (pot_matrix_QP.*array))
    energy = kin + pot
    return kin
end

#spread auxillary functions

function r_2_array(r_2_array,n)
    for x in 1:n
        for y in 1:n
            r_2_array[x,y] = (x-45)^2 + (y-45)^2
        end
    end
    return r_2_array
end

function x_array(x_array, n)
    for x in 1:n
        for y in 1:n
            x_array[x,y] = (x-45)
        end
    end
    return x_array
end

function y_array(y_array, n)
    for x in 1:n
        for y in 1:n
            y_array[x,y] = (y-45)
        end
    end
    return y_array
end


r_2_matrix = r_2_array(zeros(89, 89), 89)
x_matrix = x_array(zeros(89, 89), 89)
y_matrix = y_array(zeros(89, 89), 89)

function expec_value(array, thing)
    return real(sum(conj(array[:,:,1]).*(thing.*array[:,:,1]))) + real(sum(conj(array[:,:,2]).*(thing.*array[:,:,2])))
end

#The spread of the wavefunction
function spread(array)
    return expec_value(array, r_2_matrix) - (expec_value(array, x_matrix))^2 - (expec_value(array, y_matrix))^2
end

function functionize(psi, x, y, s)
    return real(dot(psi[x,y,s], psi[x,y,s]))
end

#Plotters____________________________________________________
using Plots; pyplot()
function Plotter(psi, del_t, t)
    x = (1:10000)/40
    #array = load("C:/Users/Alucard/Desktop/julia/data_sets/spread_L_89_10000_1-40.jld", "data")
    array = spread.(time_evolve_fast(psi,t,del_t,init_FFT(zeros(ComplexF64, 89, 89, 2), 1:2),init_IFFT(zeros(ComplexF64, 89, 89, 2), 1:2)))
    plot!(x, array, xaxis = :log, yaxis = :log)
end

using JLD

function data(t, psi, del_t)
    array = time_evolve_fast(psi,t,del_t,init_FFT(zeros(ComplexF64, 89, 89, 2), 1:2),init_IFFT(zeros(ComplexF64, 89, 89, 2), 1:2))
    save("C:/Users/Alucard/Desktop/julia/data_sets/density_L_89_10000_1-40.jld", "data", array)
end

# using ProgressMeter
#
# prog = Progress(10000,1)
#
# array = load("C:/Users/Alucard/Desktop/julia/data_sets/density_L_89_10000_1-40.jld", "data")
# anim = @animate for i=1:10000
#     x=1:89
#     y=1:89
#     z(x,y) = functionize(array[i], x, y,1)
#     plot(x,y,z,st=:surface,camera=(-30,30))
#     next!(prog)
# end
# gif(anim, "C:/Users/Alucard/Desktop/julia/density_anim_AD/anim_L_89_10000_1-40.gif", fps = 30)

#data(10000,psi_guess_array(Array{ComplexF64}(undef, 89,89,2), 89), 40^-1)

#Plotter(psi_guess_array(Array{ComplexF64}(undef, 89,89,2), 89), 40^-1, 10000)


#function time_evolve(array, t, F_T, I_F_T)
#psi = psi_guess_array(Array{ComplexF64}(undef, 89,89,2),89)
#x = time_evolve(psi_guess_array(Array{ComplexF64}(undef, 89,89,2),89), 100, init_FFT(zeros(ComplexF64, 89, 89, 2), 1:2),init_IFFT(zeros(ComplexF64, 89, 89, 2), 1:2), 40^-1)
#Energy(psi, init_FFT(zeros(ComplexF64, 89, 89, 2), 1:2), init_IFFT(zeros(ComplexF64, 89, 89, 2), 1:2))
#Energy(x, init_FFT(zeros(ComplexF64, 89, 89, 2), 1:2), init_IFFT(zeros(ComplexF64, 89, 89, 2), 1:2))
