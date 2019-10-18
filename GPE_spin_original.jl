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
    return exp(-((x-72)^2) - ((y-72)^2))
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

function psi_guess_array_dir(psi_guess_array, n)
    for x in 1:n
        for y in 1:n
            psi_guess_array[x,y] = psi_guess(x,y)
        end
    end
    return normalizer(psi_guess_array)
end

#normalizes the wavefunction array
#array[:,:,1].* conj(array[:,:,1])) + (array[:,:,2].*conj(array[:,:,2])
function normalizer(array)
    s = sqrt(dot(array, array))
    return (1000*array/s)
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

spin_matrix = kin_spin_matrix(zeros(ComplexF64, 144, 144, 2, 2), 144, -im*40^-1)

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
function time_step_V(array, del_t, pot_matrix_QP, g)
    array[:,:,1] = array[:,:,1].*(exp.((-im*del_t)*(pot_matrix_QP + (g*(conj(array[:,:,1]).*array[:,:,1])))))
    array[:,:,2] = array[:,:,2].*(exp.((-im*del_t)*(pot_matrix_QP + (g*(conj(array[:,:,2]).*array[:,:,2])))))
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

FFT = init_FFT(144, 1:2)

IFFT = init_IFFT(144, 1:2)

function time_evolve_step(array, del_t, pot_matrix_QP, g)
    return IFFT*(time_step_T(FFT*time_step_V(IFFT*time_step_T(FFT*array,144,x_dummy), del_t, pot_matrix_QP, g),144,x_dummy))
end

#evolves the guess function array t steps in imaginary time
function time_evolve(array, t, del_t)
    evolved_array = reduce((x, y) -> time_evolve_step(x, del_t), 1:t, init=array)
    return evolved_array
end

using ProgressMeter

function time_evolve_fast(array, t, del_t, pot_matrix_QP, g)

    n_array = [array]

    @progress for i in 2:t
        push!(n_array, time_evolve_step(n_array[i-1], del_t, pot_matrix_QP, g))
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
function Energy(array)
    k_array = init_FFT(zeros(ComplexF64, 144, 144, 2), 1:2)*array
    kin = expec_value(k_array, kin_spin_matrix)
    pot = expec_value(array, pot_matrix_QP)
    return pot + kin
end

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


r_2_matrix = r_2_array(zeros(144, 144), 144)
x_matrix = x_array(zeros(144, 144), 144)
y_matrix = y_array(zeros(144, 144), 144)

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

function pot_error(t)

    n_array = [psi]

    @progress for i in 2:t
        push!(n_array, time_step_V(n_array[i-1], 10^-2))
    end
    return n_array
end

#DIRECT INTEGRATION FUNCTIONS__________________________________________________

#Buildss the Hamiltonian-------------------------------------------
function Ham_up(k_x, k_y, t)
    cos(sqrt(sin(k_x)^2 + sin(k_y)^2)*t)
end

function Ham_down(k_x, k_y, t)
    if sin(k_x)^2 + sin(k_y)^2 == 0
        return 0
    else
        return ((sin(k_y) - im*sin(k_x))*sin(sqrt(sin(k_x)^2 + sin(k_y)^2)*t))/(sqrt(sin(k_x)^2 + sin(k_y)^2))
    end
end

#Builds psi--------------------------------------------------------------------
function psi_guess_array_dir(psi_guess_array, n)
    for x in 1:n
        for y in 1:n
            psi_guess_array[x,y,1] = psi_guess(x,y)
            psi_guess_array[x,y,2] = 0
        end
    end
    return normalizer(psi_guess_array)
end

x_dummy = zeros(ComplexF64, 144, 144, 2)
psi_k_ = FFT*psi_guess_array_dir(x_dummy, 144)


function psi_k_t(array, n, t)
    for x in 1:n
        for y in 1:n
            array[x,y,1] = psi_k_[x,y,1]*Ham_up(p_x(x,n),p_y(y,n),t)

            array[x,y,2] = psi_k_[x,y,1]*Ham_down(p_x(x,n),p_y(y,n),t)
        end
    end
    return array
end

function psi_x_t(n, t)
    return IFFT*psi_k_t(x_dummy, n, t)
end
#______________________________________________________________________________

#Plotters____________________________________________________
using Plots
function Plotter_dir(t)
    x = (1:t)/40
    array = Float64[]
    @progress for i in 1:t
        push!(array, spread(psi_x_t(144, i/40)))
    end
    #array = load("C:/Users/Alucard/Desktop/julia/data_sets/spread_L_144_10000_1-40.jld", "data")
    plot(x, array, xaxis = :log, yaxis = :log)
end

using DataFrames
using GLM
function Plotter(t)
    x = (1:t)/40
                                #(pot_array, W, phi_x, phi_y, n, m)
    pot_matrix_QP = pot_array_QP(zeros(ComplexF64, 144, 144), .1, 0, 0, 144, 55)
    #(array, t, del_t, pot_matrix_QP, g)
    array = spread.(time_evolve_fast(psi_guess_array(x_dummy, 144), t, pot_matrix_QP, -im*40^-1, 0))
    plot(x, array)
    #data = DataFrame(A = x[400:t], B = array[400:t])
    #linear = glm(@formula(log(B) ~ log(A)), data, Normal(), IdentityLink())
    # t = (1:10000)/40
    # x = [.1,.15,.2,.25,.3,.35,.4,.45,.5,.52,.54,.56,.6,.65,.7,.75,.8,.85,.9]
    # for i in 1:19
    #     array = load("D:/GPE_data/L_144_spread/spread_L_144_10000_1-40_W_$(x[i]).jld", "data")
    #     plot!(t, array, xaxis = :log, yaxis = :log)
    # end
end

using JLD

function data(t)
    x = [.1,.15,.2,.25,.3,.35,.4,.45,.5,.52,.54,.56,.6,.65,.7,.75,.8,.85,.9]
    for i in 1:19
        pot_matrix_QP = pot_array_QP(zeros(ComplexF64, 144, 144), x[i], 0, 0, 144, 55)
        array = spread.(time_evolve_fast(psi_guess_array(x_dummy, 144), t, -im*40^-1, pot_matrix_QP))
        save("D:/GPE_data/spread_L_144_10000_1-40_W_$(x[i]).jld", "data", array)
    end
end

function lin_data()
    x = [.1,.15,.2,.25,.3,.35,.4,.45,.5,.52,.54,.56,.6,.65,.7,.75,.8,.85,.9]
    t = (1:10000)/40
    m_array = Float64[]
    for i in 1:19
        array = load("D:/GPE_data/L_144_spread/spread_L_144_10000_1-40_W_$(x[i]).jld", "data")
        data = DataFrame(A = t[400:4000], B = array[400:4000])
        linear = glm(@formula(log(B) ~ log(A)), data, Normal(), IdentityLink())
        m = coef(linear)[2]
        push!(m_array,m)
        #save("D:/GPE_data/L_144_spread_lin_fit/spread_L_144_1-40_W_$(x[i])_lin_fit.jld", "data", array)
    end
    y = [1.82,1.75,1.66, 1.56, 1.47, 1.35, 1.25, 1.15, 1.12, 1.12, 1.12, 1.13, 1.17, 1.23, 1.3, 1.37, 1.41, 1.44, 1.46]
    scatter(x, m_array, title ="2/z vs. W Plot", xlabel = "W", ylabel = "2/z", markerstrokecolor = :black)
    plot!(x, m_array)
    plot!(x, y)
end

function Norm(array)
    return dot(array,array)
end

function E_Plot(t)
    x = 1:144
    y = 1:144
    #(pot_array, W, phi_x, phi_y, n, m)
    pot_matrix_QP = pot_array_QP(zeros(ComplexF64, 144, 144), 0, 0, 0, 144, 55)
    #(array, t, del_t, pot_matrix_QP, g)
    array = time_evolve_fast(psi_guess_array(x_dummy, 144), t, -im*40^-1,  pot_matrix_QP, 1000)[t]
    arra = array[:,:,1]
    z(x,y) = conj(arra[x,y]).*arra[x,y]

    surface(x,y,z)
end

# using ProgressMeter
#
# prog = Progress(10000,1)
#
# array = load("C:/Users/Alucard/Desktop/julia/data_sets/density_L_144_10000_1-40.jld", "data")
# anim = @animate for i=1:10000
#     x=1:144
#     y=1:144
#     z(x,y) = functionize(array[i], x, y,1)
#     plot(x,y,z,st=:surface,camera=(-30,30))
#     next!(prog)
# end
# gif(anim, "C:/Users/Alucard/Desktop/julia/density_anim_AD/anim_L_144_10000_1-40.gif", fps = 30)



#data(10000,psi_guess_array(Array{ComplexF64}(undef, 144,144,2), 144), -im*40^-1)



# x = [.1,.15,.2,.25,.3,.35,.4,.45,.5,.52,.54,.56,.6,.65,.7,.75,.8,.85,.9]
# y = [1.82,1.75,1.66, 1.56, 1.47, 1.35, 1.25, 1.15, 1.12, 1.12, 1.12, 1.13, 1.17, 1.23, 1.3, 1.37, 1.41, 1.44, 1.46]
# plot(x,y, title ="2/z vs. W Plot", xlabel = "W", ylabel = "2/z", markerstrokecolor = :black)
# scatter!(x,y)


# array = load("D:/GPE_data/L_144_spread/spread_L_144_10000_1-40_W_0.9.jld", "data")
# plot(ts, array, xaxis = :log, yaxis = :log, legend=false, title ="Spread Plot with Linear Fit for W = 0.9", xlabel = "t", ylabel = "spread")
# datas = DataFrame(A = ts[1000:4000], B = array[1000:4000])
# linear = glm(@formula(log(B) ~ log(A)), datas, Normal(), IdentityLink())
# ts = (1:10000)/40
# f(x) = ((x^m))
# plot!(f,1,250)
# m = coef(linear)[2]
# b = coef(linear)[1]
