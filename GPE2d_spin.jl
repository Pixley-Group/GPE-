#Library of functions for solving the 2D GPE.

#input parameters--------------------------

#size of step in real time
del_t = (40^-1)

fib(n) = n < 2 ? n : fib(n-1) + fib(n-2)

#range over real space (-x_end:x_end)
n = fib(12)
m = fib(10)
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
function psi_guess(x,y,z)
    return exp(-((x-73)^2) - ((y-73)^2))
end

#write the potential energy here
function pot_H_O(x,y)
    pot = m*(omega^2)*(x^2)+(y^2)
end

#Write the Quasi-Periodic Potential Energy here
function pot_QP(x,y)
    W = 1
    phi_x = 0
    phi_y = 0
    Q = (2/pi)*(m/n)
    return W*(cos((Q*x)+phi_x) + cos((Q*y)+phi_y))
end


#core Functions for the algortihm___________________________________

using LinearAlgebra, Statistics, Compat
#generates the array that is the discretization of the guess function in real space
function psi_guess_array()
    psi_guess_array = zeros(ComplexF64, n, n, 2)

    for x in 1:n
        for y in 1:n
            psi_guess_array[x,y,1] = psi_guess(x,y,-1)
            psi_guess_array[x,y,2] = psi_guess(x,y,1)
        end
    end
    return normalizer(psi_guess_array)
end

#normalizes the wavefunction array
#array[:,:,1].* conj(array[:,:,1])) + (array[:,:,2].*conj(array[:,:,2])
function normalizer(array)
    s = sqrt(sum(array.*conj(array)))
    return (array/s)
end

#-----------------POTENTIAL ENERGY

#Generates the potential energy array
function pot_array_H_O(psi)
    pot_array = zeros(ComplexF64, n, n)

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



function pot_array_QP()
    pot_array = zeros(ComplexF64, n, n)

    for i in 1:n
        for j in 1:n
            pot_array[i,j] = pot_QP(i,j)
        end
    end

    return pot_array

end

pot_matrix_QP = pot_array_QP()

#---------------KINETIC ENERGY

#calculates the x momentum for the fft minted momentum eignestate
function p_x(n)
    return ((2*pi)*(n-1))/n
end

#calculates the y momentum for the fft minted momentum eignestate
function p_y(m)
    return ((2*pi)*(m-1))/n
end

#calculates the total kinetic energy for each fft minted momentum eigenstate
function kin_mom(p_x, p_y)
    return ((p_x)^2+(p_y)^2)
end


#Generates the momentum Kinetic Energy array
function kin_mom_array()
    kin_array = zeros(ComplexF64, n, n)

    for i in 1:n
        for j in 1:n
            kin_array[n,m] = kin_mom(p_x(i), p_y(j))
        end
    end
    return kin_array
end

#generates the spin orbital coupling matrix


function kin_spin_matrix()
    spin_couple_matrix = zeros(ComplexF64, n, n, 2, 2)
    A(x,y) = sin(p_x(x)) - im*sin(p_y(y))
    for x in 1:n
        for y in 1:n
            spin_couple_matrix[x,y,1,1] = cos(abs(A(x,y)) * del_t/2)
            spin_couple_matrix[x,y,1,2] = -im*exp(im * angle(A(x,y))) * sin(abs(A(x,y)) * del_t/2)
            spin_couple_matrix[x,y,2,1] = -im*exp(-im * angle(A(x,y))) * sin(abs(A(x,y)) * del_t/2)
            spin_couple_matrix[x,y,2,2] = cos(abs(A(x,y)) * del_t/2)
        end
    end
    return spin_couple_matrix
end

spin_matrix = kin_spin_matrix()

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
gen_array = zeros(ComplexF64, n, n, 2)
function time_step_T(array)
    for x in 1:n
        for y in 1:n
            gen_array[x,y,1] = spin_matrix[x,y,1,1]*array[x,y,1] + spin_matrix[x,y,1,2]*array[x,y,2]
            gen_array[x,y,2] = spin_matrix[x,y,2,1]*array[x,y,1] + spin_matrix[x,y,2,2]*array[x,y,2]
        end
    end
    return gen_array
end

#evolves psi delt_t in time with the PE operator
function time_step_V(array)
    return array.*(exp.(pot_matrix_QP*(-im*del_t)))
end

#evolves the guess function array one step in imaginary time
using FFTW
using AbstractFFTs

function init_FFT()
    A = zeros(ComplexF64, n, n, 2)
    region = 1:2
    return plan_fft(A,region; flags=FFTW.PATIENT, timelimit=Inf)
end

function init_IFFT()
    A = zeros(ComplexF64, n, n, 2)
    region = 1:2
    return plan_ifft(A,region; flags=FFTW.PATIENT, timelimit=Inf)
end

F_T = init_FFT()
I_F_T = init_IFFT()

function time_evolve_step(array)
    x = fftshift(F_T*array)
    x = time_step_T(x)
    x = I_F_T*x
    x = time_step_V(x)
    x = fftshift(F_T*x)
    x = time_step_T(x)
    return I_F_T*x
end

#evolves the guess function array t steps in imaginary time
function time_evolve(array, t)
    evolved_array = reduce((x, y) -> time_evolve_step(x), 1:t, init=array)
    return evolved_array
end

function time_evolve_fast(array,t)

    n_array = [array]

    for i in 2:t
        push!(n_array, time_evolve_step(n_array[i-1]))
    end
    return n_array
end

function T(x,y)
    return [[0 sin(p_x(x)) - (im* sin(p_y(y)))]; [sin(p_x(x)) + im * sin(p_y(y)) 0]]
end

function T_step(array)
    gen_array = zeros(ComplexF64, n, n, 2)
    for x in 1:n
        for y in 1:n
            gen_array[x,y,:] = T(x,y)*array[x,y,:]
        end
    end
    return gen_array
end

#Finds the energy of psi
function Energy(array)
    ffts = F_T*array
    kin = dot(array, I_F_T*T_step(ffts))
    pot = dot(array, (pot_matrix_QP.*array))
    energy = kin + pot
    return kin
end

#spread auxillary functions
function r(x,y)
    return r = sqrt((x*conj(x)) + (y*conj(y)))
end

function r_2(x,y)
    return r = (x*conj(x)) + (y*conj(y))
end

function r_array()
    r_array = zeros(n, n)
    for x in 1:n
        for y in 1:n
            r_array[x,y] = r(x,y)
        end
    end

    return r_array
end

function r_2_array()
    r_2_array = zeros(n, n)
    for x in 1:n
        for y in 1:n
            r_2_array[x,y] = r_2(x,y)
        end
    end

    return r_2_array
end

r_matrix = r_array()
r_2_matrix = r_2_array()

#The spread of the wavefunction
function spread(psi)

    psi_r = conj(psi).*r_matrix.*psi
    psi_r_2 = conj(psi).*r_2_matrix.*psi

    return real(sum(psi_r_2) - (sum(psi_r))^2)

end


#Plotters____________________________________________________
using Plots
function Plotter()
    t = 1:100
    psi = psi_guess_array()
    spread_array = zeros(100)
    @progress for t in 1:100
        spread_array[t] = spread(time_evolve_fast(psi,t)[t])
    end
    plot!(t, spread_array, xaxis =:log, yaxis =:log, legend = false)
end
