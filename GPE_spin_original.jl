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
    #return (x-72)-((x-72)^2)+(x-72)^3
end

#write the potential energy here
function pot_H_O(x,y,M,omega)
    return (M/2)*(omega^2)*(((x-72)^2)+((y-72)^2))
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
    return array/sqrt(dot(array, array))
end

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

spin_matrix = kin_spin_matrix(zeros(ComplexF64, 233, 233, 2, 2), 233, (40^-1))

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

function interactions(g, array, del_t)
    exp.(g*(L^2)*(-im*del_t)*(abs2(array[:,:1])+abs2(array[:,:,2])))
end

#evolves psi delt_t in time with the PE operator
function time_step_V(array, del_t, pot_matrix_QP, pot_matrix_HO, g)
    array[:,:,1] = exp.(((-im*del_t)*(pot_matrix_QP)).*interactions(g, array, del_t).*array[:,:,1])
    array[:,:,2] = exp.(((-im*del_t)*(pot_matrix_QP)).*interactions(g, array, del_t).*array[:,:,2])
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

FFT = init_FFT(233, 1:2)

IFFT = init_IFFT(233, 1:2)

function time_evolve_step(array, del_t, pot_matrix_QP, pot_matrix_HO, g)
    return IFFT*(time_step_T(FFT*time_step_V(IFFT*time_step_T(FFT*array,233,x_dummy), del_t, pot_matrix_QP, pot_matrix_HO, g),233,x_dummy))
end

#evolves the guess function array t steps in imaginary time
function time_evolve(array, t, del_t)
    evolved_array = reduce((x, y) -> time_evolve_step(x, del_t), 1:t, init=array)
    return evolved_array
end

using ProgressMeter

function time_evolve_fast(array, t, del_t, pot_matrix_QP, pot_matrix_HO, g)

    n_array = [array]

    @progress for i in 2:t
        push!(n_array, time_evolve_step(n_array[i-1], del_t, pot_matrix_QP, pot_matrix_HO, g))
    end
    return n_array
end

function spread_fast(array, t, del_t, pot_matrix_QP, pot_matrix_HO, g)

    n_array = [array]
    spread_array = [spread(array)]

    @progress for i in 2:t
        array = time_evolve_step(n_array[1], del_t, pot_matrix_QP, pot_matrix_HO, g)
        push!(n_array, array)
        push!(spread_array, spread(array))
        popfirst!(n_array)
    end
    return spread_array
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
    k_array = init_FFT(zeros(ComplexF64, 233, 233, 2), 1:2)*array
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


r_2_matrix = r_2_array(zeros(233, 233), 233)
x_matrix = x_array(zeros(233, 233), 233)
y_matrix = y_array(zeros(233, 233), 233)

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
        push!(n_array, time_step_V(n_array[i-1], 10^-3))
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

x_dummy = zeros(ComplexF64, 233, 233, 2)
psi_k_ = FFT*psi_guess_array_dir(x_dummy, 233)


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
function Plotter_dir(t)
    x = (1:t)/40
    array = Float64[]
    @progress for i in 1:t
        push!(array, spread(psi_x_t(233, i/40)))
    end
    #array = load("C:/Users/Alucard/Desktop/julia/data_sets/spread_L_233_10000_1-40.jld", "data")
    plot(x, array, xaxis = :log, yaxis = :log)
end

using DataFrames
using GLM
function Plotter(t)
    x = (1:t)/40
    pot_matrix_QP = pot_array_QP(zeros(ComplexF64, 233, 233), 0.1, 0, 0, 233, 89)
    pot_matrix_HO = pot_array_H_O(233, zeros(ComplexF64, 233, 233), 0, 0)
    array = spread.(time_evolve_fast(psi_guess_array(x_dummy, 233), t, 40^-1,  pot_matrix_QP, pot_matrix_HO, 0))
    plot(x, array, xaxis = :log, yaxis = :log)
    #data = DataFrame(A = x[400:t], B = array[400:t])
    #linear = glm(@formula(log(B) ~ log(A)), data, Normal(), IdentityLink())
    # t = (1:10000)/40
    # x = [.1,.15,.2,.25,.3,.35,.4,.45,.5,.52,.54,.56,.6,.65,.7,.75,.8,.85,.9]
    # for i in 1:19
    #     array = load("D:/GPE_data/L_233_spread/spread_L_233_10000_1-40_W_$(x[i]).jld", "data")
    #     plot!(t, array, xaxis = :log, yaxis = :log)
    # end
end

using JLD

function data_amarel(t)
    x = [.1,.15,.2,.25,.3,.35,.4,.45,.5,.52,.54,.56,.6,.65,.7,.75,.8,.85,.9]
    for i in 1:19
        pot_matrix_QP = pot_array_QP(zeros(ComplexF64, 233, 233), 0, 0, 0, 233, 89)
        pot_matrix_HO = pot_array_H_O(233, zeros(ComplexF64, 233, 233), 0, 0)
        array = spread.(time_evolve_fast(psi_guess_array(x_dummy, 233), t, 40^-1,  pot_matrix_QP, pot_matrix_HO, x[i]))
        save("D:/GPE_data/L_233_spread_W_0_N_1/N-1_W-0_g_$(x[i]).jld", "data", array)
    end
end

function data(t)
    x = [.1,.15,.2,.25,.3,.35,.4,.45,.5,.52,.54,.56,.6,.65,.7,.75,.8,.85,.9]
    for i in 1:19
        pot_matrix_QP = pot_array_QP(zeros(ComplexF64, 233, 233), x[i], 0, 0, 233, 89)
        pot_matrix_HO = pot_array_H_O(233, zeros(ComplexF64, 233, 233), 0, 0)
        array = spread.(time_evolve_fast(psi_guess_array(x_dummy, 233), t, 40^-1,  pot_matrix_QP, pot_matrix_HO, 0.2))
        save("D:/GPE_data/L_233_spread_g_0.2/W_$(x[i])_g_0.2.jld", "data", array)
    end
end

function lin_data()
    x = [.1,.15,.2,.25,.3,.35,.4,.45,.5,.52,.54,.56,.6,.65,.7,.75,.8,.85,.9]
    t = (1:10000)/40
    m_array = Float64[]
    m_array2 = Float64[]
    for i in 1:19
        array = load("C:/Users/Kyle/Desktop/julia/data/L=233_N=233_W=0/N=234_W=0_g=$(x[i]).jld", "data")
        #array2 = load("D:/GPE_data/L_233_spread_g_0.1/W_$(x[i])_g_0.2.jld", "data")
        data = DataFrame(A = t[40:1000], B = array[40:4000])
        #data2 = DataFrame(A2 = t[40:1000], B2 = array2[40:1000])
        linear = glm(@formula(log(B) ~ log(A)), data, Normal(), IdentityLink())
        #linear2 = glm(@formula(log(B2) ~ log(A2)), data2, Normal(), IdentityLink())
        m = coef(linear)[2]
        m2 = coef(linear2)[2]
        push!(m_array,m)
        push!(m_array2,m2)
        #save("D:/GPE_data/L_233_spread_lin_fit/spread_L_233_1-40_W_$(x[i])_lin_fit.jld", "data", array)
    end
    y = [1.82,1.75,1.66, 1.56, 1.47, 1.35, 1.25, 1.15, 1.12, 1.12, 1.12, 1.13, 1.17, 1.23, 1.3, 1.37, 1.41, 1.44, 1.46]
    scatter(x, m_array, title ="2/z vs. W Plot", xlabel = "W", ylabel = "2/z", markerstrokecolor = :black)
    plot!(x, m_array)
    plot!(x, y)
    plot!(x, m_array2)
end

function spread(i)
    x = [.1,.15,.2,.25,.3,.35,.4,.45,.5,.52,.54,.56,.6,.65,.7,.75,.8,.85,.9]
    t = (1:10000)/40
    m_array = Float64[]
    m_array2 = Float64[]

    # array = load("D:/GPE_data/L_233_spread_g_0.2/W_$(x[i])_g_0.2.jld", "data")
    # array2 = load("D:/GPE_data/L_233_spread_g_0.1/W_$(x[i])_g_0.2.jld", "data")

    plot(t, array, xaxis = :log, yaxis = :log)
    plot!(t, array2, xaxis = :log, yaxis = :log)
end

function Norm(array)
    return dot(array,array)
end




using Plots
import PyPlot
using Distributions
using LinearAlgebra
#pygui(true)
#.54,.56,.6,.65,.7,.75,.8,.85,.9
function spread_data()
    y = ["20^-1", "40^-1", "10^-2", "20^-2"]
    x = ['1', '2', '3', '4']
    t = [(1:20000)/20,(1:40000)/40, (1:100000)/100, (1:200000)/200]
    for i in 1:4
        array = load("C:/Users/Kyle/Desktop/julia/data/Spread_L=233_N=1_W=0_$(y[i])/$(y[i])_N=1_W=0_g=0.2.jld", "data")
        t = t[i]
        plot!(array, fmt = :png, xlabel = "iterations", ylabel = "<del r^2>", label = "$(y[i])", title = "Spread Plot (log-log) del_t=20^-1", xaxis = :log, yaxis = :log)
    end
end

function E_Plot(t, mu)
    pyplot()

    x = (1:233)
    y = (1:233)

    #(pot_array, W, phi_x, phi_y, n, m)
    pot_matrix_QP = pot_array_QP(zeros(ComplexF64, 233, 233), 0.1, 0, 0, 233, 89)
    pot_matrix_HO = pot_array_H_O(233, zeros(ComplexF64, 233, 233), 0, 0)
    #(array, t, del_t, pot_matrix_QP, g)
    array = spread_fast(psi_guess_array(x_dummy, 233), t, 40^-1,  pot_matrix_QP, pot_matrix_HO, .1)[t]
    arra = normalizer(array[:,:,1])
    arra2 = normalizer(array[:,:,2])
    psi = psi_guess_array(x_dummy, 233)
    z = (conj.(arra[x,y]).*arra[x,y]) + (conj.(arra2[x,y]).*arra2[x,y])

    PyPlot.surf(z, rstride=2,edgecolors="k", cstride=2, alpha=1, linewidth=0.25)
    show()
    #z_2(x,y)= mu*ones(Float64, 233).- ((x).^2).-(y.^2)
    #plot!(x,y,z_2, st=:surface,colors = "green", camera=(90,30))
    # plot_surface(x, y, z_2, rstride=2,edgecolors="k", cstride=2,
    # cmap=ColorMap("gray"), alpha=0.8, linewidth=0.25)


    #surface!(x,y,z_2, alpha = 0.2)
end

function TF_Plot(mu)
    x = 1:233
    y = 1:233
    pot_matrix_HO = pot_array_H_O(233, zeros(Float64, 233, 233), 1, 1)

    z_2= (mu*ones(Float64, 233, 233) - pot_matrix_HO)/(300)
    analytic = PyPlot.surf(x,y,z_2, rstride=2,edgecolors="k", cstride=2, alpha=0.4, linewidth=0.25)
    show(analytic)
end



# pyplot()
# ioff()
# x = range(0; stop=2*pi, length=1000); y = sin.(3 * x + 4 * cos.(2 * x))
# thing = PyPlot.plot(x, y, color="blue", linewidth=2.0, linestyle="--")
# title("A new sinusoidally modulated sinusoid")
#
# show()

# using ProgressMeter
#
# prog = Progress(10000,1)
#
# array = load("C:/Users/Alucard/Desktop/julia/data_sets/density_L_233_10000_1-40.jld", "data")
# anim = @animate for i=1:10000
#     x=1:233
#     y=1:233
#     z(x,y) = functionize(array[i], x, y,1)
#     plot(x,y,z,st=:surface,camera=(-30,30))
#     next!(prog)
# end
# gif(anim, "C:/Users/Alucard/Desktop/julia/density_anim_AD/anim_L_233_10000_1-40.gif", fps = 30)



#data(10000,psi_guess_array(Array{ComplexF64}(undef, 233,233,2), 233), -im*40^-1)



# x = [.1,.15,.2,.25,.3,.35,.4,.45,.5,.52,.54,.56,.6,.65,.7,.75,.8,.85,.9]
# y = [1.82,1.75,1.66, 1.56, 1.47, 1.35, 1.25, 1.15, 1.12, 1.12, 1.12, 1.13, 1.17, 1.23, 1.3, 1.37, 1.41, 1.44, 1.46]
# plot(x,y, title ="2/z vs. W Plot", xlabel = "W", ylabel = "2/z", markerstrokecolor = :black)
# scatter!(x,y)

#ts = (1:100)/40
#array = load("D:/GPE_data/L_233_spread_g_0.1/W_0.54_g_0.1.jld", "data")
#plot(ts, array, legend=false, title ="Spread Plot with Linear Fit for W = 0.9", xlabel = "t", ylabel = "spread")
# datas = DataFrame(A = ts[1000:4000], B = array[1000:4000])
# linear = glm(@formula(log(B) ~ log(A)), datas, Normal(), IdentityLink())

# f(x) = ((x^m))
# plot!(f,1,250)
# m = coef(linear)[2]
# b = coef(linear)[1]

#E_Plot(2000, 10)
#TF_Plot(389)

#---
using Plots
using JLD
#y = ["20^-1", "40^-1", "10^-2", "20^-2", "40^-2"]
t = (1:499999)/100
f = [20, 40, 100, 200, 400]

g = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
j = 3

L = [89,144,233]

plot()
for i in 1:3
    array = load("C:/Users/Kyle/Desktop/julia/data/Spread_L=$(L[i])_W=0_0.01_500000/L=$(L[i])_W=0_g=$(g[j]).jld", "data")
    #array2 = load("C:/Users/Kyle/Desktop/julia/data/L=144_N=144_W=0_20^-1/20^-1_N=144^2_W=0_g=0.0.jld", "data")
    println(array[1])
    plot!(t, array[2:end].-array[1], fmt = :png, legend = :topleft, xlabel = "time", ylabel = "<del r^2>", label = "L = $(L[i])", title = "Spread Plot (log-log) W= 0, del_t=10^-2, g = $(g[j])", xaxis = :log, yaxis = :log)
end
#array2 = load("C:/Users/Kyle/Desktop/julia/data/L_377_spread_W=0_0.05_400000_0-0.2/0.05_N=377^2_W=0_g=0.0.jld", "data")
#plot!(t[j], array2, fmt = :png, legend = :bottomleft, xlabel = "time", ylabel = "<del r^2>", label = "g = 0", title = "Spread Plot (log-log) W= 0, del_t=20^-1, N = L^2= 377^2", xaxis = :log, yaxis = :log)
plot!()
#data_frame = DataFrame(A = t[j][400*f[j]:2500*f[j]], B = array[400*f[j]:2500*f[j]])
#linear = glm(@formula(B ~ A), data_frame, Normal(), IdentityLink())
#m = coef(linear)[2]
#b = coef(linear)[1]
# plot!(t[j][1000:10000*f[j]], ((t[j][1000:10000*f[j]]*exp(b/m)).^m) , label = "2/z = $(m)", legend = :topleft)
#plot!(t[j][3000:10000*f[j]], (t[j][3000:10000*f[j]].*m.+b), label = "2/z = $(m)", legend = :topleft)
#---
