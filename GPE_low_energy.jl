using LinearAlgebra, Statistics, Compat

#normalizes the function
function normalizer(array)
    s = sqrt(dot(array, array))
    return (array/s)
end

#Fourier Transform setup
using FFTW
using AbstractFFTs

function init_FFT(n, region)
    return plan_fft(zeros(ComplexF64, n, n),region; flags=FFTW.PATIENT, timelimit=Inf)
end

function init_IFFT(n, region)
    return plan_ifft(zeros(ComplexF64, n, n, 2),region; flags=FFTW.PATIENT, timelimit=Inf)
end

FFT = init_FFT(89, 1:2)
IFFT = init_IFFT(89, 1:2)

#Generate psi_0--------------------------------------------------
function psi_guess(x,y)
    return exp(-((x-45)^2) - ((y-45)^2))
end

function psi_guess_array_dir(psi_guess_array, n)
    for x in 1:n
        for y in 1:n
            psi_guess_array[x,y] = psi_guess(x,y)
        end
    end
    return normalizer(psi_guess_array)
end

x_dummy = zeros(ComplexF64, 89, 89, 2)
x_dummy_dir = zeros(ComplexF64, 89, 89)
psi_k_ = FFT*psi_guess_array_dir(x_dummy_dir, 89)

#------------------------------------------------------------------

#calculates the x momentum for the fft minted momentum eignestate
function p_x(x,n)
    return ((2*pi)*(x-1))/n
end

#calculates the y momentum for the fft minted momentum eignestate
function p_y(y,n)
    return ((2*pi)*(y-1))/n
end


#Buildss the Hamiltonian-------------------------------------------


function Ham_up(x, y, t, n)
    cos(sqrt(p_x(x,n)^2 + p_y(y,n)^2)*t)
end

function Ham_down(x, y, n, t)
    if p_x(x,n)^2 + p_y(y,n)^2 == 0
        return 0
    else
        return ((p_y(y,n) - im*p_x(x,n))*sin(sqrt(p_x(x,n)^2 + p_y(y,n)^2)*t))/(sqrt(p_x(x,n)^2 + p_y(y,n)^2))
    end
end

#--------------------------------------------------------------------
function psi_k_t(array, n, t)
    for x in 1:n
        for y in 1:n
            array[x,y,1] = psi_k_[x,y]*Ham_up(x,y,t,n)

            array[x,y,2] = psi_k_[x,y]*Ham_down(x,y,t,n)
        end
    end
    return array
end

function psi_x_t(n, t)
    return IFFT*psi_k_t(x_dummy, n, t)
end

#expec_value-------------------------------------------------------------
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

function r_2_array(r_2_array,n)
    for x in 1:n
        for y in 1:n
            r_2_array[x,y] = (x)^2 + (y)^2
        end
    end
    return r_2_array
end

r_2_matrix = r_2_array(x_dummy_dir, 89)
x_matrix = x_array(x_dummy_dir, 89)
y_matrix = y_array(x_dummy_dir, 89)

function expec_value(array, thing)
    return real(sum(conj(array[:,:,1]).*(thing.*array[:,:,1]))) + real(sum(conj(array[:,:,2]).*(thing.*array[:,:,2])))
end

function spread(array)
    return expec_value(array, r_2_matrix) - (expec_value(array, x_matrix))^2 - (expec_value(array, y_matrix))^2
end

using ProgressMeter
using Plots
function Plotter_low(t)
    x = (1:t)/40
    y = 1:89
    array = Float64[]
    @progress for i in 1:t
        push!(array, spread(psi_x_t(89, i/40)))
    end
    #array = load("C:/Users/Alucard/Desktop/julia/data_sets/spread_L_89_10000_1-40.jld", "data")
    plot(x, array)
    plot(p_x(y,89))
    plot(p_y(y,89))
end

#Plotter(1000)

function Ham_low(x, y, n, t)
    return exp(-im*t*[0 (im*p_x(x,n)-p_y(y,n)); (p_x(x,n)+im*p_y(y,n)) 0 ])
end

Plotter_low(1000)
