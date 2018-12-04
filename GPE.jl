#input parameters--------------------------

#size of step in imaginary time
delt_t = -(10^-3)*im

#size of step in real space
delt_x = 10^-3

#range over real space (-x_end:x_end)
x_end = 1

#The number of steps in real space
#N_x = Int(floor(2*x_end/delt_x)) # calculated value (default)
N_x = (10^3) #manual value
x_array = LinRange(-x_end,x_end,N_x)

#physical constants
k = 1
m = 1
#strength of particle-particle interactions (coupling parameter)
g = 5*10^5
#number of particles
N = 10^3
#---------------------------------------------------


#write the guess wavefunction here
function psi_guess(x)
    return pi
end


#-----------------------------------------------------------

#core Functions for the algortihm___________________________________

#generates the array that is the discretization of the guess function in real space
function psi_guess_array()
    psi_guess_array = Float64[]

    for x in x_array
        push!(psi_guess_array,psi_guess(x))
    end
    return normalize(psi_guess_array)
end

#The potential energy function
function pot(x, psi)
    pot = ((k/2)*(x^2)+g*(conj(psi)*psi))
    return pot
end

#The potential energy operator
function e_V(x, psi)
    return exp.(-pot(x, psi)*(delt_t)*im)
end

#The kinetic energy operator with standard dispersion in momentum space
function e_T(n)
    p = ((2*pi)*n)/N_x
    return exp(-((p^2)/2m)*(delt_t/2)*im)
end

#applies the Kinetic Energy operator once
function time_step_T(array)
    psi_k_T = ComplexF64[]

    for n = 1:N_x
        k_T = array[n]*e_T(n)
        push!(psi_k_T, k_T)
    end

    return psi_k_T
end

#applies the Potential Energy operator once
function time_step_V(array)
    return_array = ComplexF64[]
    for x = 1:N_x
        k_V = array[x]*e_V(x_array[x], array[x])
        push!(return_array, k_V)
    end
    return return_array
end


#evolves the guess function array one step in imaginary time
using FFTW
function time_evolve_step(array)
    psi_k_1 = fftshift(fft(array))
    psi_k_T_1 = time_step_T(psi_k_1)
    psi_x_T_1 = ifft(psi_k_T_1)
    psi_x_V = time_step_V(psi_x_T_1)
    psi_k = fftshift(fft(psi_x_V))
    psi_k_T = time_step_T(psi_k)
    psi_x_T = ifft(psi_k_T)
    return normalize(psi_x_T)
end

#evolves the guess function in imaginary time t times (over t iteration)
function time_evolve(array, t)
    evolved_array = reduce((x, y) -> time_evolve_step(x), 1:t, init=array)
    return normalize(evolved_array)
end

#normalizes an array
function normalize(array)
    a = array .* conj.(array)
    s = sqrt(sum(a))
    return (array/s)*sqrt(N)
end
##################################

#Plotters____________________________________________________
using Plots
#plots guess psi
psi = psi_guess_array()
time_step = time_evolve_step(psi)
t = time_evolve(psi, 10000)
t_2 = conj.(t) .* t

#plots a static plot of the real modulus of the evolved psi
plot(x_array, real(t_2), title = "Psi Evolved")

#uncomment to create animations of the evolution of the guess function psi
# anim = @animate for i=1:800
#     t = time_evolve(psi, i)
#     plot(x_array, real(conj.(t) .* t), title = "Psi Evolved")
# end
# gif(anim, "C:/Users/Alucard/Desktop/julia/gifs/SE_e^-x_fps10.gif", fps = 10)
