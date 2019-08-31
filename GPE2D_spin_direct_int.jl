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
#------------------------------------------------------------------

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

#--------------------------------------------------------------------

function psi_k(n)
    return init_FFT(n, 1:2)*psi_guess_array_dir(zeros(ComplexF64, n, n), n)
end

function psi_k_t(array, n, t)
    for x in 1:n
        for y in 1:n
            array[x,y,1] = psi_k(n)[x,y]*Ham_up(x,y,t)

            array[x,y,2] = psi_k(n)[x,y]*Ham_down(x,y,t)
        end
    end
    return array
end

function psi_x_t(n, t)
    return init_IFFT(n, 1:2)*psi_k_t(zeros(ComplexF64, 89, 89, 2), n, t)
end

#expec_value-------------------------------------------------------------
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

function r_2_array(r_2_array,n)
    for x in 1:n
        for y in 1:n
            r_2_array[x,y] = (x-45)^2 + (y-45)^2
        end
    end
    return r_2_array
end

function expec_value(array, thing)
    return real(sum(conj(array[:,:,1]).*(thing.*array[:,:,1]))) + real(sum(conj(array[:,:,2]).*(thing.*array[:,:,2])))
end

function spread(array)
    return expec_value(array, r_2_matrix) - (expec_value(array, x_matrix))^2 - (expec_value(array, y_matrix))^2
end
