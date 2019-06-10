#Library of functions for solving the 2D GPE.

#input parameters--------------------------

#size of step in real time
delt_t = (10^-3)

fib(n) = n < 2 ? n : fib(n-1) + fib(n-2)

#range over real space (-x_end:x_end)
n = 12
N_x = fib(n)
N_y = fib(n)

#x and y real space array generator
x_array = collect(1:fib(n))
y_array = collect(1:fib(n))

#V HO "spring constant"
m = 1
omega = 0
#strength of particle-particle interactions (V coupling constant)
g = 0
# Energy scale for Pot QP
W = 1
#phase for Pot QP
phi_x = 0
phi_y = 0
#spin-orbital coupling constant
s = 1
#momentum kinetic energy strength
t = 0
#number of particles
N = 1
#magnitude of test Gaussian
A = 1
#---------------------------------------------------

#----Enter the starting parameters here
# Write the guess wavefunction here
function psi_guess(x,y,z)
    return A*exp(-(x^2) - (y^2))
end

#write the potential energy here
function pot_H_O(x,y)
    pot = m*(omega^2)*(x^2)+(y^2)
end

#Write the Quasi-Periodic Potential Energy here
function pot_QP(x,y)
    Q = (2/pi)*(fib(n-2)/fib(n))
    return pot_QP = W*(cos((Q*x)+phi_x) + cos((Q*y)+phi_y))
end


#core Functions for the algortihm___________________________________

using LinearAlgebra, Statistics, Compat
#generates the array that is the discretization of the guess function in real space
function psi_guess_array()
    psi_guess_array = zeros(ComplexF64, N_x, N_y, 2)

    for x in 1:N_x
        for y in 1:N_y
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
    return (array/s).*sqrt(N)
end

#-----------------POTENTIAL ENERGY

#Generates the potential energy array
function pot_array_H_O(psi)
    pot_array = zeros(ComplexF64, N_x, N_y)

    for n in 1:N_x
        for m in 1:N_y
            pot_array[n,m] = pot_H_O(x_array[n], y_array[m])
        end
    end

    psi = psi[:,:,1]+psi[:,:,2]
    psi = reshape(psi, N_x, N_y)

    pot_array = pot_array+g*(conj(psi).*psi)
    return pot_array

end

function pot_array_QP()
    pot_array = zeros(ComplexF64, N_x, N_y)

    for n in 1:N_x
        for m in 1:N_y
            pot_array[n,m] = pot_QP(x_array[n], y_array[m])
        end
    end

    return pot_array

end

#---------------KINETIC ENERGY

#calculates the x momentum for the fft minted momentum eignestate
function p_x(n)
    return p_x = ((2*pi)*n-1)/N_x
end

#calculates the y momentum for the fft minted momentum eignestate
function p_y(m)
    return p_y = ((2*pi)*m-1)/N_y
end

#calculates the total kinetic energy for each fft minted momentum eigenstate
function kin_mom(p_x, p_y)
    return ((p_x)^2+(p_y)^2)
end


#Generates the momentum Kinetic Energy array
function kin_mom_array()
    kin_array = zeros(ComplexF64, N_x, N_y)

    for n in 1:N_x
        for m in 1:N_y
            kin_array[n,m] = kin_mom(p_x(n), p_y(m))
        end
    end
    return kin_array
end

#generates the spin orbital coupling matrix
function kin_spin_matrix(x, y)
    spin_couple_matrix = zeros(ComplexF64,2,2)
    A = sin(p_x(x)) - im*sin(p_y(y))
    B = sin(p_x(x)) + im*sin(p_y(y))
    spin_couple_matrix[1,1] = s*cos(sqrt(A*B)*delt_t/2)
    spin_couple_matrix[1,2] = -s*im*sqrt(A)/sqrt(B)*sin(sqrt(A*B)*delt_t/2)
    spin_couple_matrix[2,1] = -s*im*sqrt(B)/sqrt(A)*sin(sqrt(A*B)*delt_t/2)
    spin_couple_matrix[2,2] = s*cos(sqrt(A*B)*delt_t/2)
    return spin_couple_matrix
end

#exponentiates the harmonic oscilator pot energy operator elementwise
function e_V(psi)
    return exp.(((-m*(omega)^2)/2)*pot_array_H_O(psi).*((-delt_t)*im))
end

#exponentiates the momentum kin energy operator elementwise
function e_T()
    return e_T = exp.(-t*kin_mom_array()*(delt_t/2)*im)
end

#evolves psi delt_t in time with the KE operator
function time_step_T(array)
    array3 = zeros(ComplexF64, N_x, N_y, 2)
    for x in 1:N_x
        for y in 1:N_y
            array3[x,y,1] = kin_spin_matrix(x,y)[1,1]*array[x,y,1] + kin_spin_matrix(x,y)[1,2]*array[x,y,2]
            array3[x,y,2] = kin_spin_matrix(x,y)[2,1]*array[x,y,1] + kin_spin_matrix(x,y)[2,2]*array[x,y,2]
        end
    end
    return array3
end

#evolves psi delt_t in time with the PE operator
function time_step_V(array)
    array[:,:,1] = array[:,:,1].*(exp.(pot_array_QP()*(-im*delt_t)))
    array[:,:,2] = array[:,:,2].*(exp.(pot_array_QP()*(-im*delt_t)))
    return array
end

#evolves the guess function array one step in imaginary time
using FFTW
using AbstractFFTs

function init_FFT()
    A = zeros(ComplexF64, N_x,N_y,2)
    region = 1:2
    F_T = plan_fft(A,region; flags=FFTW.PATIENT, timelimit=Inf)
    return F_T
end

function init_IFFT()
    A = zeros(ComplexF64, N_x,N_y,2)
    region = 1:2
    I_F_T = plan_ifft(A,region; flags=FFTW.PATIENT, timelimit=Inf)
    return I_F_T
end

F_T = init_FFT()
I_F_T = init_IFFT()

function time_evolve_step(array)
    x = F_T*array
    x = time_step_T(x)
    x = I_F_T*x
    x = time_step_V(x)
    x = F_T*x
    x = time_step_T(x)
    x = I_F_T*x
    return x
end

#evolves the guess function array t steps in imaginary time
function time_evolve(array, t)
    evolved_array = reduce((x, y) -> time_evolve_step(x), 1:t, init=array)
    return evolved_array
end

#Finds the energy of psi
# function Energy(array)
#     ffts = fftshift(F_T*array)
#     kin = I_F_T*(ffts.*kin_array())
#     pot = psi.*pot_array_H_O(psi)
#     energy = (kin.+pot).*conj(psi)
#     return sum(energy)
# end

#spread auxillary functions
function r(x,y)
    return r = sqrt((x*conj(x)) + (y*conj(y)))
end

function r_2(x,y)
    return r = (x*conj(x)) + (y*conj(y))
end

function r_array()
    r_array = zeros(N_x, N_y)
    for x in 1:N_x
        for y in 1:N_y
            r_array[x,y] = r(x,y)
        end
    end

    return r_array
end

function r_2_array()
    r_2_array = zeros(N_x, N_y)
    for x in 1:N_x
        for y in 1:N_y
            r_2_array[x,y] = r_2(x,y)
        end
    end

    return r_2_array
end

#The spread of the wavefunction
function spread(psi)

    psi_r_up = (conj(psi[:,:,1]).*r_array()).*psi[:,:,1]
    psi_r_down = (conj(psi[:,:,2]).*r_array()).*psi[:,:,2]
    psi_r_2_up = (conj(psi[:,:,1]).*r_2_array()).*psi[:,:,1]
    psi_r_2_down = (conj(psi[:,:,2]).*r_2_array()).*psi[:,:,2]
    psi_r = psi_r_up + psi_r_down
    psi_r_2 = psi_r_2_up + psi_r_2_down

    return spread = real(sum(psi_r_2) - (sum(psi_r))^2)

end


#Plotters____________________________________________________
using Plots

t = 1:50
psi = psi_guess_array()
spread_array = zeros(50)
for t in 1:50
    spread_array[t] = spread(time_evolve(psi,t))
end
plot!( spread_array, xaxis=:log, yaxis=:log)
