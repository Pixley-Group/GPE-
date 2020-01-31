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
    #return (x-72)-((x-72)^2)+(x-72)^3
end

#write the potential energy here
# function pot_H_O(x,y,M,omega)
#     return (M/2)*(omega^2)*(((x-116)^2)+((y-116)^2))
# end
#
# #Write the Quasi-Periodic Potential Energy here
# function pot_QP(x, y, W, phi_x, phi_y, n, m)
#     return W*(cos((((2*pi)*(m/n))*(x))+phi_x) + cos((((2*pi)*(m/n))*(y))+phi_y))
# end


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

# function psi_guess_array_dir(psi_guess_array, n)
#     for x in 1:n
#         for y in 1:n
#             psi_guess_array[x,y] = psi_guess(x,y)
#         end
#     end
#     return normalizer(psi_guess_array)
# end

#normalizes the wavefunction array
#array[:,:,1].* conj(array[:,:,1])) + (array[:,:,2].*conj(array[:,:,2])
function normalizer(array)
    return array/sqrt(dot(array, array))
end

# function prob_norm(array)
#     return array/sqrt(dot(array, array))
# end
#-----------------POTENTIAL ENERGY

#Generates the potential energy array
function pot_array_H_O(n, pot_array, M, omega)
    for i in 1:n
        for j in 1:n
            pot_array[i,j] = pot_H_O(i, j, M, omega)
        end
    end
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
# function kin_mom(p_x, p_y)
#     return ((p_x)^2+(p_y)^2)
# end


#Generates the momentum Kinetic Energy array
# function kin_mom_array(kin_array, n)
#
#     for i in 1:n
#         for j in 1:n
#             kin_array[i,j] = kin_mom(p_x(i,n), p_y(j,n))
#         end
#     end
#     return kin_array
# end

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



#exponentiates the harmonic oscilator pot energy operator elementwise
# function e_V(psi)
#     return exp.(((-m*(omega)^2)/2)*pot_array_H_O(psi).*((-del_t)*im))
# end
#
# #exponentiates the momentum kin energy operator elementwise
# function e_T()
#     return e_T = exp.(-t*kin_mom_array()*(del_t/2)*im)
# end

#evolves psi delt_t in time with the KE operator
#evaluate this one time and store in memory

function time_step_T(array,n, gen_array, spin_matrix)
    for x in 1:n
        for y in 1:n
            gen_array[x,y,1] = spin_matrix[x,y,1,1]*array[x,y,1] + spin_matrix[x,y,1,2]*array[x,y,2]
            gen_array[x,y,2] = spin_matrix[x,y,2,1]*array[x,y,1] + spin_matrix[x,y,2,2]*array[x,y,2]
        end
    end
    return gen_array
end

function interactions(g, array, del_t, L)
    exp.(g*(L^2)*(-im*del_t)*(abs2(array[:,:1])+abs2(array[:,:,2])))
end

#evolves psi delt_t in time with the PE operator -TOOK OUT QUASIPERIODIC POTENTIAL W=0
function time_step_V(array, del_t, g, L)
    # array .= array.*interactions(g,array,del_t,L)
    array[:,:,1] = array[:,:,1].*interactions(g, array, del_t, L)
    array[:,:,2] = array[:,:,2].*interactions(g, array, del_t, L)
    return array
end

using FFTW
using AbstractFFTs

function init_FFT(n, region)
    return plan_fft(zeros(ComplexF64, n, n, 2),region; flags=FFTW.PATIENT, timelimit=Inf)
end

function init_FFT_dir(n)
    return plan_fft(zeros(ComplexF64, n, n); flags=FFTW.PATIENT, timelimit=Inf)
end

function init_IFFT(n, region)
    return plan_ifft(zeros(ComplexF64, n, n, 2),region; flags=FFTW.PATIENT, timelimit=Inf)
end



function time_evolve_step(array, del_t, spin_matrix, g, FFT, IFFT, L)
    return IFFT*(time_step_T(FFT*time_step_V(IFFT*time_step_T(FFT*array, L, x_dummy(L), spin_matrix), del_t, g, L),L,x_dummy(L), spin_matrix))
end

#evolves the guess function array t steps in imaginary time
# function time_evolve(array, t, del_t)
#     evolved_array = reduce((x, y) -> time_evolve_step(x, del_t), 1:t, init=array)
#     return evolved_array
# end


function time_evolve_fast(array, t, del_t, spin_matrix, g, FFT, IFFT, L)

    n_array = [array]

    for i in 2:t
        push!(n_array, time_evolve_step(n_array[i-1], del_t, spin_matrix, g, FFT, IFFT, L))
    end
    return n_array
end

function spread_fast(array, t, del_t, spin_matrix, g, FFT, IFFT, L)

    n_array = [array]
    spread_array = [spread(array, L)]

    for i in 2:t
        #       time_evolve_step(array, del_t, spin_matrix, g, FFT, IFFT)
        array = time_evolve_step(n_array[1], del_t, spin_matrix, g, FFT, IFFT, L)
        push!(n_array, array)
        push!(spread_array, spread(array, L))
        popfirst!(n_array)
    end
    return spread_array
end

# function T(x,y,n)
#     return [[0 sin(p_x(x,n)) - (im* sin(p_y(y,n)))]; [sin(p_x(x,n)) + im * sin(p_y(y,n)) 0]]
# end
#
# function T_step(array, gen_array,n)
#     for x in 1:n
#         for y in 1:n
#             gen_array[x,y,:] = T(x,y,n)*array[x,y,:]
#         end
#     end
#     return gen_array
# end

#Finds the energy of psi
# function Energy(array)
#     k_array = init_FFT(zeros(ComplexF64, 233, 233, 2), 1:2)*array
#     kin = expec_value(k_array, kin_spin_matrix)
#     pot = expec_value(array, pot_matrix_QP)
#     return pot + kin
# end

#spread auxillary functions

function r_2_array(r_2_array,n)
    for x in 1:n
        for y in 1:n
            r_2_array[x,y] = (x)^2 + (y)^2
        end
    end
    return r_2_array
end

function x_array(x_array, n)
    for x in 1:n
        for y in 1:n
            x_array[x,y] = (x)
        end
    end
    return x_array
end

function y_array(y_array, n)
    for x in 1:n
        for y in 1:n
            y_array[x,y] = (y)
        end
    end
    return y_array
end

function expec_value(array, thing)
    #array = prob_norm(array)
    return real(sum(conj(array[:,:,1]).*(thing.*array[:,:,1]))) + real(sum(conj(array[:,:,2]).*(thing.*array[:,:,2])))
end

#The spread of the wavefunction
function spread(array, L)
    r_2_matrix = r_2_array(zeros(L, L), L)
    x_matrix = x_array(zeros(L, L), L)
    y_matrix = y_array(zeros(L, L), L)
    return expec_value(array, r_2_matrix) - (expec_value(array, x_matrix))^2 - (expec_value(array, y_matrix))^2
end

function functionize(psi, x, y, s)
    return real(dot(psi[x,y,s], psi[x,y,s]))
end

# function pot_error(t)
#
#     n_array = [psi]
#
#     for i in 2:t
#         push!(n_array, time_step_V(n_array[i-1], 10^-3))
#     end
#     return n_array
# end

#DIRECT INTEGRATION FUNCTIONS__________________________________________________

#Buildss the Hamiltonian-------------------------------------------
# function Ham_up(k_x, k_y, t)
#     cos(sqrt(sin(k_x)^2 + sin(k_y)^2)*t)
# end
#
# function Ham_down(k_x, k_y, t)
#     if sin(k_x)^2 + sin(k_y)^2 == 0
#         return 0
#     else
#         return ((sin(k_y) - im*sin(k_x))*sin(sqrt(sin(k_x)^2 + sin(k_y)^2)*t))/(sqrt(sin(k_x)^2 + sin(k_y)^2))
#     end
# end
#
# #Builds psi--------------------------------------------------------------------
# function psi_guess_array_dir(psi_guess_array, n)
#     for x in 1:n
#         for y in 1:n
#             psi_guess_array[x,y,1] = psi_guess(x,y)
#             psi_guess_array[x,y,2] = 0
#         end
#     end
#     return normalizer(psi_guess_array)
# end

function x_dummy(L)
    return zeros(ComplexF64, L, L, 2)
end
#psi_k_ = FFT*psi_guess_array_dir(x_dummy, 233)


# function psi_k_t(array, n, t)
#     for x in 1:n
#         for y in 1:n
#             array[x,y,1] = psi_k_[x,y,1]*Ham_up(p_x(x,n),p_y(y,n),t)
#
#             array[x,y,2] = psi_k_[x,y,1]*Ham_down(p_x(x,n),p_y(y,n),t)
#         end
#     end
#     return array
# end
#
# function psi_x_t(n, t)
#     return IFFT*psi_k_t(x_dummy, n, t)
# end
#______________________________________________________________________________

using DataFrames
using JLD
import MPI

using MPI
function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    Nproc = MPI.Comm_size(comm)

    data_amarel(500000, 89, rank, Nproc, (10^-2))

    MPI.Finalize()
end



#function hw(rank, Nproc)
    #print("Hello, world from processor $(rank), of $(Nproc)! \n")
#end

#function avg(rank, Nproc, comm)
#    sum = MPI.Reduce(rank, +, 0, comm)
 #   if rank == 0
  #      avg = sum/Nproc
   #     print("The sum of our $(Nproc) ranks is $(sum) \n")
    #    print("The average of all of our ranks is $(avg)")
    #end
#end

function data_amarel(t, L, rank, Nproc, del_t)

    x = [0,.2,.4,.6,.8,1]
    #pot_matrix_QP = pot_array_QP(zeros(ComplexF64, L, L), 0, 0, 0, L, 233)
    #pot_matrix_HO = pot_array_H_O(L, zeros(ComplexF64, L, L), 0, 0)
    spin_matrix = kin_spin_matrix(zeros(ComplexF64, L, L, 2, 2), L, del_t)
    FFT = init_FFT(L, 1:2)
    IFFT = init_IFFT(L, 1:2)
    #time_evolve_fast(array, t, del_t, spin_matrix, g, FFT, IFFT)
    #        spread_fast(array, t, del_t, spin_matrix, g, FFT, IFFT, L)
    #        time_evolve_fast(array, t, del_t, spin_matrix, g, FFT, IFFT)
    array = (time_evolve_fast(psi_guess_array(x_dummy(L), L), t, del_t, spin_matrix, x[rank+1], FFT, IFFT))
    save("/scratch/bomeisl/Spread_L_$(L)_W=0_$(del_t)_$(t)/L=$(L)_W=0_g=$(x[rank+1]).jld", "data", array)
end

# function data(t)
#     x = [.1,.15,.2,.25,.3,.35,.4,.45,.5,.52,.54,.56,.6,.65,.7,.75,.8,.85,.9]
#     for i in 1:19
#         pot_matrix_QP = pot_array_QP(zeros(ComplexF64, 233, 233), x[i], 0, 0, 233, 89)
#         pot_matrix_HO = pot_array_H_O(233, zeros(ComplexF64, 233, 233), 0, 0)
#         array = spread.(time_evolve_fast(psi_guess_array(x_dummy, 233), t, 10^-2,  pot_matrix_QP, pot_matrix_HO, 0.2))
#         save("D:/GPE_data/L_233_spread_g_0.2/W_$(x[i])_g_0.2.jld", "data", array)
#     end
# end


main()
