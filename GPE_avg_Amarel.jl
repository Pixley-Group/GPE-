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
function psi_guess(x,y, L)
    return exp(-((x-floor(Int,L/2))^2)/1 - ((y-floor(Int, L/2))^2))/1
    #return (x-72)-((x-72)^2)+(x-72)^3
end

#write the potential energy here
# function pot_H_O(x,y,M,omega)
#     return (M/2)*(omega^2)*(((x-116)^2)+((y-116)^2))
# end
#
# #Write the Quasi-Periodic Potential Energy here
# function pot_QP(x, y, phi_x, phi_y, n, m)
#     return (cos((((2*pi)*(m/n))*(x))+phi_x) + cos((((2*pi)*(m/n))*(y))+phi_y))
# end


#core Functions for the algortihm___________________________________

using LinearAlgebra, Statistics, Compat
#generates the array that is the discretization of the guess function in real space
function psi_guess_array(psi_guess_array, L)
    for x in 1:L
        for y in 1:L
            psi_guess_array[x,y,1] = psi_guess(x,y,L)
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

function pot_array_QP(pot_array, phi_x, phi_y, n, m)
    for i in 1:n
        for j in 1:n
            pot_array[i,j] = cos((((2*pi)*(m/n))*(i))+phi_x) + cos((((2*pi)*(m/n))*(j))+phi_y)
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


#SOC KE MATRIX
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

#DECOUPLED KE MATRIX
#function kin_spin_matrix_decoup(spin_couple_matrix, n, del_t)
    #A(x,y,n) = -cos(p_x(x,n)) - cos(p_y(y,n))
    #for x in 1:n
        #for y in 1:n
            #spin_couple_matrix[x,y,1,1] = exp(-im*A(x,y,n)*(del_t/2))
            #spin_couple_matrix[x,y,1,2] = 0
            #spin_couple_matrix[x,y,2,1] = 0
            #spin_couple_matrix[x,y,2,2] = exp(-im*A(x,y,n)*(del_t/2))
        #end
    #end
    #return spin_couple_matrix
#end


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

#evolves psi delt_t in time with the PE operator
function time_step_V(array, del_t, g, L, pot_matrix_QP, W)
    return array .= array .*  (exp.(g*(1)*(-im*del_t)*(abs2.(array[:,:,1]) + abs2.(array[:,:,2])))) .* exp.((-im*del_t)*W*pot_matrix_QP)
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



function time_evolve_step(array, del_t, spin_matrix, g, FFT, IFFT, L, pot_matrix_QP, W)
    return IFFT*(time_step_T(FFT*time_step_V(IFFT*time_step_T(FFT*array, L, x_dummy(L), spin_matrix), del_t, g, L, pot_matrix_QP, W),L,x_dummy(L), spin_matrix))
end

#evolves the guess function array t steps in imaginary time
# function time_evolve(array, t, del_t)
#     evolved_array = reduce((x, y) -> time_evolve_step(x, del_t), 1:t, init=array)
#     return evolved_array
# end


function time_evolve_fast_thin(array, t, del_t, spin_matrix, g, FFT, IFFT, L, thin_fac, pot_matrix_QP, W)
    return_array = [array]
    n_array = [array]
    for j in 1:(Int(t/thin_fac))
        for i in 2:thin_fac
            push!(n_array, time_evolve_step(n_array[1], del_t, spin_matrix, g, FFT, IFFT, L, pot_matrix_QP, W))
            deleteat!(n_array, 1)
        end
        push!(return_array, n_array[1])
    end
    return return_array
end

function time_evolve_fast(array, t, del_t, spin_matrix, g, FFT, IFFT, L, pot_matrix_QP, W)
    return_array = [array]
    for j in 2:t
        push!(return_array, time_evolve_step(return_array[j-1], del_t, spin_matrix, g, FFT, IFFT, L, pot_matrix_QP, W))
    end
    return return_array
end


function spread_fast(array, t, del_t, spin_matrix, g, FFT, IFFT, L, pot_matrix_QP, W)

    n_array = [array]
    spread_array = [spread(array, L)]

    for i in 2:t
        #       time_evolve_step(array, del_t, spin_matrix, g, FFT, IFFT)
        array = time_evolve_step(n_array[1], del_t, spin_matrix, g, FFT, IFFT, L, pot_matrix_QP, W)
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

# function g(g, del_t, array)
#     output = []
#     for i in 1:length(array)
#         push!(output, g*del_t*(abs2.(array[i][:,:,1]) + abs2.(array[i][:,:,2])))
#     end
#     return output
# end

function T_array(T_array, L)
    A(x,y,L) = sin(p_x(x,L)) + im*sin(p_y(y,L))
    for x in 1:L
        for y in 1:L
            T_array[x,y,1] = round(A(x,y,L), digits = 7)
            T_array[x,y,2] = round(conj(A(x,y,L)), digits = 7)
        end
    end
    return T_array
end

function Energy(array, T_array, L, g, FFT)
    T = round(dot((FFT*array), T_array.*(FFT*array)), digits = 6)
    U = (g/2) * (dot(array,array))^2
    return(T+U)
end

function IPR(psi,L)
    sum = 0
    for x in 1:L
        for y in 1:L
            sum += ((conj(psi[x,y,1])*psi[x,y,1]) + (conj(psi[x,y,2])*psi[x,y,2]))^2
        end
    end
    return sum
end

function PR(psi,L)
    sum = 0
    for x in 1:L
        for y in 1:L
            sum += ((conj(psi[x,y,1])*psi[x,y,1]) + (conj(psi[x,y,2])*psi[x,y,2]))^2
        end
    end
    return 1/sum
end


function PR_Spread_fast(array, t, del_t, spin_matrix, g, FFT, IFFT, L, pot_matrix_QP, W)

    n_array = [array]
    PR_array = [PR(array,L)]
    Spread_array = [spread(array, L)]

    for i in 2:t
        #       time_evolve_step(array, del_t, spin_matrix, g, FFT, IFFT)
        array = time_evolve_step(n_array[1], del_t, spin_matrix, g, FFT, IFFT, L, pot_matrix_QP, W)
        push!(n_array, array)
        push!(PR_array, PR(array,L))
        push!(Spread_array, spread(array, L))
        popfirst!(n_array)
    end
    return PR_array, Spread_array
end

function Energy_fast(array, T_array, t, del_t, spin_matrix, g, FFT, IFFT, L, pot_matrix_QP, W)

    n_array = [array]
    Energy_array = [Energy(array, T_array, L, g, FFT)]

    for i in 2:t
        #       time_evolve_step(array, del_t, spin_matrix, g, FFT, IFFT)
       	array = time_evolve_step(n_array[1], del_t, spin_matrix, g, FFT, IFFT, L, pot_matrix_QP, W)
        push!(n_array, array)
        push!(Energy_array, Energy(array, T_array, L, g, FFT))
        popfirst!(n_array)
    end
    return Energy_array
end

function Energy_load(array, T_array, L, g, FFT)
    Enerray = []

    for i in 1:101
        append!(Enerray, Energy(array[i], T_array, L, g, FFT))
    end

    return Enerray
end

function IPR_load(array, L)
    Enerray = []

    for i in 1:101
        append!(Enerray, IPR(array[i], L))
    end

    return Enerray
end

function PR_load(array, L)
    Enerray = []

    for i in 1:101
        append!(Enerray, (PR(array[i], L)))
    end

    return Enerray
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

#Builds the Hamiltonian-------------------------------------------
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
using Distributions
using Random
using MPI

function main()
    fib(n) = n < 2 ? n : fib(n-1) + fib(n-2)
    L = fib(11)
    del_t = 0.0005
    t = 400000
    spin_matrix = kin_spin_matrix(zeros(ComplexF64, L, L, 2, 2), L, del_t)
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    Nproc = MPI.Comm_size(comm)
    #Psi_L=144_W=0_0.01_100000_g*1
    PR, Spread = data_amarel(L, fib(9), t, rank, Nproc, 0.0005, 1, parse(Float64, ARGS[1]), parse(Float64, ARGS[2]), init_FFT(L, 1:2), init_IFFT(L, 1:2), spin_matrix)
    sum_PR = MPI.Reduce(PR, +, 0, comm)
    sum_Spread = MPI.Reduce(Spread, +, 0, comm)

    if rank == 0
        average_PR = sum_PR/Nproc
        save("/scratch/bomeisl/9-3-20/PR_L=$(L)_$(del_t)_$(t)_gauss-1_avg_SOC/PR_L=$(L)_g=$(parse(Float64, ARGS[1]))_W=$(parse(Float64, ARGS[2]))_gauss-1_avg_SOC.jld", "data", average_PR)
        average_Spread = sum_Spread/Nproc
        save("/scratch/bomeisl/9-3-20/Spread_L=$(L)_$(del_t)_$(t)_gauss-1_avg_SOC/Spread_L=$(L)_g=$(parse(Float64, ARGS[1]))_W=$(parse(Float64, ARGS[2]))_gauss-1_avg_SOC.jld", "data", average_Spread)
    end

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

function data_amarel(L, m, t, rank, Nproc, del_t, thin_fac, g, W, FFT, IFFT, spin_matrix)
    fib(n) = n < 2 ? n : fib(n-1) + fib(n-2)
    Random.seed!(rank)
    seed = rand(Uniform(0,2*pi),2)
    #g = [10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15]
    #x = [0.01,0.001, 10, 20]
    #pot_matrix_QP = pot_array_QP(zeros(ComplexF64, L, L), 0, 0, 0, L, L)
    #pot_matrix_HO = pot_array_H_O(L, zeros(ComplexF64, L, L), 0, 0)
    #               pot_array_QP(pot_array, W, phi_x, phi_y, n, m)
    pot_matrix_QP = pot_array_QP(zeros(ComplexF64, L, L), seed[1], seed[2], L, m)

    return PR_Spread_fast(psi_guess_array(x_dummy(L), L), t, del_t, spin_matrix, parse(Float64, ARGS[1]), FFT, IFFT, L, pot_matrix_QP, parse(Float64, ARGS[2]))
    #        time_evolve_fast(array, t, del_t, spin_matrix, g, FFT, IFFT, L, thin_fac)
    #return out =  PR_array = PR_fast(psi_guess_array(x_dummy(L), L), t, del_t, spin_matrix, g, FFT, IFFT, L, pot_matrix_QP, W)
    #E_array = Energy_fast(psi_guess_array(x_dummy(L), L), T_array(x_dummy(L), L), t, del_t, spin_matrix, g[rank+1], FFT, IFFT, L, pot_matrix_QP, W)
    #Psi_array = time_evolve_fast_thin(psi_guess_array(x_dummy(L), L), t, del_t, spin_matrix, g[rank+1], FFT, IFFT, L, 6000, pot_matrix_QP, W)

                                                       #Psi_L=244_W=0_0.0005_600000_gauss-1 #Psi_L=244_W=0_g=1000_gauss-1_data.jld                                                  #Energy(array, T_array, L, g, FFT)
    #E_array = Energy_load(load("/scratch/bomeisl/4-2-20/Psi_L=$(L)_W=$(W)_$(del_t)_$(t)_gauss-1/Psi_L=$(L)_W=$(W)_g=$(x[rank+1])_gauss-1_data.jld", "data"), T_array(zeros(ComplexF64, L, L, 2), L), L, x[rank+1], FFT)
    #PR_array = PR_load(load("/scratch/bomeisl/4-2-20/Psi_L=$(L)_W=$(W)_$(del_t)_$(t)_gauss-1/Psi_L=$(L)_W=$(W)_g=$(x[rank+1])_gauss-1_data.jld", "data"), L)

    #save("/scratch/bomeisl/5-24-20/W=$(W)/PR_L=$(L)_g=$(g)_W=$(W)_gauss-1.jld", "data", PR_array)
    #save("/scratch/bomeisl/5-18-20/Energy_L=$(L)_W=0_$(del_t)_$(t)_gauss-10/Energy_L=$(L)_W=0_g=$(x[rank+1])_gauss-10.jld", "data", E_array)
    #save("/scratch/bomeisl/5-18-20/Psi_L=$(L)_W=0_$(del_t)_$(t)_gauss-1_smallg/Psi_L=$(L)_W=0_g=$(x[rank+1])_gauss-1_smallg.jld", "data", Psi_array)
   #Psi_L=89_W=0_0.001_300000_DECOUP -- Psi_L=89_W=0_g=-1.0_DECOUP.jld
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
