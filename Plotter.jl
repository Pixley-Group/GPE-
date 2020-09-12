using Plots

using JLD

# plot()
#
# G = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
# w = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.0]
# for i = 1:9
#
#     for j = 1:11
#         plot!(
#             (1:600000) * 0.0005,
#             real(load(
#             "/home/kyle/GPE-Data/6-13-20/PR_L=233_0.0005_600000_gauss-1_avg_SOC/PR_L=233_g=$(G[j])_W=$(w[i])_gauss-1_avg.jld", "data"
#             )),
#             legend = :topleft,
#             title = "PR vs. Time SOC, L=144, Variance = 1, W=$(w[i]), del_t=0.005",
#             xlabel = "Time",
#             ylabel = "PR",
#             label = "g = $(G[j])",
#             xaxis = :log,
#             yaxis = :log,
#             )
#         end
#      savefig("PR-SOC_W=$(w[i])_100")
#      plot()
# end

# plot()
# j=11
#
# G = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
# w = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#
# for i = 1:10
#     plot!(
#          (1:600000) * 0.0005,
#          real(load(
#          "/home/kyle/GPE-Data/6-13-20/Spread_L=233_0.0005_600000_gauss-1_avg_SOC/Spread_L=233_g=$(G[j])_W=$(w[i])_gauss-1_avg.jld", "data"
#          )),
#          legend = :topleft,
#          title = "Spread vs. Time SOC, L=233, Variance = 1, g=$(G[j]), del_t=0.0005",
#          xlabel = "Time",
#          ylabel = "PR",
#          label = "W = $(w[i])",
#          xaxis = :log,
#          yaxis = :log,
#          )
# end
# savefig("Spread-SOC_g=$(G[j])_100")

using DataFrames
using GLM

# Plots Linear Fit and Curves

# Spread

G = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
w = [0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70]
#w = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2]
leng = length(w)

plot(size = (1500, 1000))

for i = 1:9

    for j = 1:7

        ts = (2:599999) * 0.0005
        array = real(load("/home/kyle/GPE-Data/9-9-20/Spread_L=144_g=$(G[i])_W=$(w[j])_gauss-1_avg_SOC.jld", "data"))
        array_ = array .- array[1]
        plot!(size = (1500, 1000), ts, array_[2:599999], title ="<\\delta r^{2} > Plot with Linear Fit for g = $(G[i])", xlabel = "time", ylabel = "<\\delta r^{2}>", legend = :topleft, label = "W = $(w[j])", xaxis= :log, yaxis= :log, titlefontsize= 20 , tickfontsize = 17, legendfontsize = 17, guidefontsize = 17)
        datas = DataFrame(A = ts[115000:499998], B = array_[115000:499998])
        linear = glm(@formula(log(B) ~ log(A)), datas, Normal(), IdentityLink())
        m = coef(linear)[2]
        b = coef(linear)[1]
        plot!(ts[20000:599998], (exp(b)*(ts[20000:599998]).^m), label = "slope = $(m)", linestyle = :dash, titlefontsize= 20 , tickfontsize = 17, legendfontsize = 17, guidefontsize = 17)
    end

    

    savefig("/home/kyle/GPE-Plots/9-9-20/<\\delta r^2>_SOC_subtr_Lin_Fit_g=$(G[i])_I_8-27-20_Fine")

    plot(size = (1500, 1000))
    
    for j = 7:15
    
        ts = (2:599999) * 0.0005
        array = real(load("/home/kyle/GPE-Data/9-9-20/Spread_L=144_g=$(G[i])_W=$(w[j])_gauss-1_avg_SOC.jld", "data"))
        array_ = array .- array[1]
        plot!(size = (1500, 1000), ts, array_[2:599999], title ="<\\delta r^{2} > Plot with Linear Fit for g = $(G[i])", xlabel = "time", ylabel = "<\\delta r^{2}>", legend = :topleft, label = "W = $(w[j])", xaxis= :log, yaxis= :log, titlefontsize= 20 , tickfontsize = 17, legendfontsize = 17, guidefontsize = 17)
        datas = DataFrame(A = ts[115000:499998], B = array_[115000:499998])
        linear = glm(@formula(log(B) ~ log(A)), datas, Normal(), IdentityLink())
    
        m = coef(linear)[2]
        b = coef(linear)[1]
        plot!(ts[20000:599998], (exp(b)*(ts[20000:599998]).^m), label = "slope = $(m)", linestyle = :dash,size = (1200, 1000), titlefontsize= 20 , tickfontsize = 17, legendfontsize = 17, guidefontsize = 17)
    end
    
    savefig("/home/kyle/GPE-Plots/9-9-20/<\\delta r^2>_SOC_subtr_Lin_Fit_g=$(G[i])_II_8-27-20_Fine")

    plot(size = (1500, 1000))
    
    for j = 15:leng
    
        ts = (2:599999) * 0.0005
        array = real(load("/home/kyle/GPE-Data/9-9-20/Spread_L=144_g=$(G[i])_W=$(w[j])_gauss-1_avg_SOC.jld", "data"))
        array_ = array .- array[1]
        plot!(size = (1500, 1000), ts, array_[2:599999], title ="<\\delta r^{2} > Plot with Linear Fit for g = $(G[i])", xlabel = "time", ylabel = "<\\delta r^{2}>", legend = :topleft, label = "W = $(w[j])", xaxis= :log, yaxis= :log, titlefontsize= 20 , tickfontsize = 17, legendfontsize = 17, guidefontsize = 17)
        datas = DataFrame(A = ts[115000:499998], B = array_[115000:499998])
        linear = glm(@formula(log(B) ~ log(A)), datas, Normal(), IdentityLink())
    
        m = coef(linear)[2]
        b = coef(linear)[1]
        plot!(ts[20000:599998], (exp(b)*(ts[20000:599998]).^m), label = "slope = $(m)", linestyle = :dash, size = (1200, 1000), titlefontsize= 20 , tickfontsize = 17, legendfontsize = 17, guidefontsize = 17)
    end
    
    savefig("/home/kyle/GPE-Plots/9-9-20/<\\delta r^2>_SOC_subtr_Lin_Fit_g=$(G[i])_III_8-27-20_Fine")
    
    plot(size = (1500, 1000))

end



## PR
# for j = 1:5
#
#     ts = (2:599999) * 0.0005
#     array = real(load("/home/kyle/GPE-Data/8-21-20/PR_L=144_g=$(G[j])_W=$(w)_gauss-1_avg_SOC.jld", "data"))
#     array_ = array .- array[1]
#     plot!(ts, array_[2:599999], title ="PR Plot with Linear Fit for W = $(w[i])", xlabel = "time", ylabel = "PR", legend = :topleft, label = "g = $(G[j])", xaxis= :log, yaxis= :log)
#     datas = DataFrame(A = ts[20000:200000], B = array_[20000:200000])
#     linear = glm(@formula(log(B) ~ log(A)), datas, Normal(), IdentityLink())
#
#     m = coef(linear)[2]
#     b = coef(linear)[1]
#     plot!(ts[2000:599998], (exp(b)*(ts[2000:599998]).^m), label = "slope = $(m)", linestyle = :dash)
# end
#
# savefig("/home/kyle/GPE-Plots/8-21-20/PR_SOC_subtr_Lin_Fit_W=$(w)_I_8-21-20")
#
# plot()
#
# for j = 5:len
#
#     ts = (2:599999) * 0.0005
#     array = real(load("/home/kyle/GPE-Data/8-21-20/PR_L=144_g=$(G[j])_W=$(w)_gauss-1_avg_SOC.jld", "data"))
#     array_ = array .- array[1]
#     plot!(ts, array_[2:599999], title ="PR Plot with Linear Fit for W = $(w[i])", xlabel = "time", ylabel = "PR", legend = :topleft, label = "g = $(G[j])", xaxis= :log, yaxis= :log)
#     datas = DataFrame(A = ts[20000:200000], B = array_[20000:200000])
#     linear = glm(@formula(log(B) ~ log(A)), datas, Normal(), IdentityLink())
#
#     m = coef(linear)[2]
#     b = coef(linear)[1]
#     plot!(ts[2000:599998], (exp(b)*(ts[2000:599998]).^m), label = "slope = $(m)", linestyle = :dash)
# end
#
# savefig("/home/kyle/GPE-Plots/8-21-20/PR_SOC_subtr_Lin_Fit_W=$(w)_II_8-21-20")




##Reading out the exponents
#
# G = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
# w = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#
# j = 11
# for i = 1:10
#
#     ts = (1:600000) * 0.0005
#     array = real(load("/home/kyle/GPE-Data/6-13-20/Spread_L=233_0.0005_600000_gauss-1_avg_SOC/Spread_L=233_g=$(G[j])_W=$(w[i])_gauss-1_avg.jld", "data"))
#     array_ = array = array[1]
#
#     datas = DataFrame(A = ts[20000:200000], B = array_[20000:200000])
#     linear = glm(@formula(log(B) ~ log(A)), datas, Normal(), IdentityLink())
#
#     m = coef(linear)[2]
#     b = coef(linear)[1]
#     plot!(ts[2000:600000], (exp(b)*(ts[2000:600000]).^m), label = "slope = $(m)")
# end
#

#

## g vs. W

# G = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
# w = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
# slope_del = []
# slope_PR = []
# plot()
#
# for i = 21
#     for j = 1:10
#         ts = (2:299999) * 0.0005
#         array = real(load("/home/kyle/GPE-Data/8-13-20/Spread_L=89_g=$(G[i])_W=$(w[j])_gauss-1_avg_SOC.jld", "data"))
#         array_ = array .- array[1]
#         datas = DataFrame(A = ts[20000:200000], B = array_[20000:200000])
#         linear = glm(@formula(log(B) ~ log(A)), datas, Normal(), IdentityLink())
#
#         m = coef(linear)[2]
#
#         append!(slope_del, m)
#     end
#
#     plot(w, slope_del, title = "Long Time Power Law Fits of <\\delta r^{2} > for g = $(G[i])", xlabel = "W", ylabel = "Power Law Fit")
#
#     savefig("/home/kyle/GPE-Plots/8-13-20/<\\delta r^2>_SOC_W_vs_Lin_Fit_g=$(G[i])_8-13-20")
#
#     for j = 1:10
#         ts = (2:299999) * 0.0005
#         array = real(load("/home/kyle/GPE-Data/8-13-20/PR_L=89_g=$(G[i])_W=$(w[j])_gauss-1_avg_SOC.jld", "data"))
#         array_ = array .- array[1]
#         datas = DataFrame(A = ts[20000:200000], B = array_[20000:200000])
#         linear = glm(@formula(log(B) ~ log(A)), datas, Normal(), IdentityLink())
#
#         m = coef(linear)[2]
#         append!(slope_PR, m)
#     end
#
#     plot(w, slope_PR, title = "Long Time Power Law Fits of PR for g = $(G[i])", xlabel = "W", ylabel = "Power Law Fit")
#
#     savefig("/home/kyle/GPE-Plots/8-13-20/PR_SOC_W_vs_Lin_Fit_g=$(G[i])_8-13-20")
# end

### Self Trapping Phase Diagram
# plot()
#
# G = [3.9, 4.75,5.45,6.15,6.9, 7.55,8.15,8.7,9.25, 9.75, 10.25]
# w = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#
# scatter(w, G, title = "Phase Diagram of PR Self Trapping Boundary", xlabel = "W", ylabel = "g\\ _{st} (W)", label = "g\\ _{st} (W)", legend = :bottomright)
# datas = DataFrame(B = G, A = w)
# linear = glm(@formula(B ~ A), datas, Normal(), IdentityLink())
# m = coef(linear)[2]
# b = coef(linear)[1]
# plot!(w, m*w .+ b, label = "g\\ _{st} (W)=g\\ _{st} (0)+$(m)*W")


#
# ### Diffusive Phase Diagram
# plot()
#
# G = []
# w = [0.0, 0.1, 0.2, 0.2, 0.3, 0.4, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9]
#
# scatter(w, G, title = "Phase Diagram of PR Self Trapping Boundary", xlabel = "W", ylabel = "g\\ _{st} (W)", label = "g\\ _{st} (W)", legend = :bottomright)
# datas = DataFrame(B = G, A = w .+ 3.9)
# linear = glm(@formula(B ~ A), datas, Normal(), IdentityLink())
# m = coef(linear)[2]
# plot!(w, m*w .+ 3.9, label = "g\\ _{st} (W)=g\\ _{st} (0)+$(m)*W")

# w = [ 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#
# i=1
#
#
# for j = 1:7
#
#         ts = (2:299999) * 0.0005
#         array = real(load("/home/kyle/GPE-Data/8-13-20/Spread_L=89_g=$(G[i])_W=$(w[j])_gauss-1_avg_SOC.jld", "data"))
#         array_ = array .- array[1]
#         plot!(size = (1500, 1000), ts, array_[2:299999], title ="<\\delta r^{2} > Plot with Linear Fit for g = $(G[i])", xlabel = "time", ylabel = "<\\delta r^{2}>", legend = :topleft, label = "W = $(w[j])", xaxis= :log, yaxis= :log, titlefontsize= 20 , tickfontsize = 17, legendfontsize = 17, guidefontsize = 17)
#         #datas = DataFrame(A = ts[60000:200000], B = array_[60000:200000])
#         #linear = glm(@formula(log(B) ~ log(A)), datas, Normal(), IdentityLink())
#         #m = coef(linear)[2]
#         #b = coef(linear)[1]
#         #plot!(ts[2000:299996], (exp(b)*(ts[2000:299996]).^m), label = "slope = $(m)", linestyle = :dash, titlefontsize= 20 , tickfontsize = 17, legendfontsize = 17, guidefontsize = 17)
# end
# plot!()
