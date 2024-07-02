using TypedTables, Dates
using Unitful, Formatting, LaTeXStrings, Measures
using Plots, StatsBase, PropDicts
using LegendHDF5IO, LegendDSP, LegendSpecFits
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using ArraysOfArrays
using Distributed
using RadiationDetectorDSP
ENV["JULIA_DEBUG"] = Main # enable debug

plotlyjs(size=(800, 600))

dsp_folder = joinpath(datafolder, "run-202311-011", "dsp2")
figures_folder = joinpath("figures", "run-202311-011")
mkpath(figures_folder)


data = fast_flatten([
            LHDataStore(
                ds -> begin
                    # @debug "Reading from \"$(ds.data_store.filename)\""
                    ds["dsp"][:]
                end,
                joinpath(dsp_folder, fk)
            ) for fk in readdir(dsp_folder) if occursin(".lh5", fk)
        ])

energy_config = readprops("config/run-202311-011/energy.json").default
th228_lines = Vector{Float64}(energy_config.th228_lines)
th228_names = Symbol.(energy_config.th228_names)
th228_names_dict  = Dict{Symbol, Float64}(Symbol.(energy_config.th228_names) .=> th228_lines)
window_sizes = Vector{Tuple{Float64, Float64}}([(l,r) for (l,r) in zip(Vector{Float64}(energy_config.left_window_sizes), Vector{Float64}(energy_config.right_window_sizes))])
n_bins = energy_config.n_bins
quantile_perc = energy_config.quantile_perc

# qc_config = readprops("config/r004/qc_config.json").default

# qc = qc_cal_energy(data, qc_config)

# @info format("Surrival Fraction: {:.2f}%", count(qc)/length(qc)*100)

data_afterQC = data[qc]
data_afterQC = data
stephist(data_afterQC.e_trap, yscale=:log10)
result_simple, report_simple = simple_calibration(data_afterQC.e_trap[data_afterQC.e_trap .< 5100], th228_lines, window_sizes,; n_bins=n_bins, quantile_perc=quantile_perc)

plot(report_simple, margin=5mm, yformatter=:plain, thickness_scaling=1.5, cal=true, title=format("ASIC BEGe Simple Calibration"))
savefig(joinpath(figures_folder, "simple_calibration_FC.png"))

m_cal_simple = result_simple.c

result_fit, report_fit = fit_peaks(result_simple.peakhists, result_simple.peakstats, th228_names)
gr()
peak_fit_plot = plot.(values(report_fit), titleloc=:center, titlefont=font(family="monospace",halign=:center, pointsize=20), ticks=:native, right_margin=10mm, top_margin=5mm, legend=false; show_label=true)
for (peak_name, p) in zip(keys(report_fit), peak_fit_plot)
    xticks!(p, convert(Int, round(xlims(p)[1], digits=0)):10:convert(Int, round(xlims(p)[2], digits=0)))
    title!(p, string(round(th228_names_dict[peak_name], digits=2)) * " keV")
end
plot(
    peak_fit_plot...,
    framestyle=:box,
    legend=:outerright,
    layout=(5,3),
    thickness_scaling=2,
    grid=true, gridalpha=0.2, gridcolor=:black, gridlinewidth=0.5,
    xguidefont=font(family="monospace",halign=:center, pointsize=18),
    yguidefont=font(family="monospace",halign=:center, pointsize=18),
    xtickfontsize=10,
    ytickfontsize=10,
    size=(4000, 4000),
    margins=1mm
)
savefig(joinpath(figures_folder, "peak_fit.png"))

μ = [result_fit[p].μ for p in th228_names] ./ m_cal_simple
μ_err = [result_fit[p].err.μ for p in th228_names] ./ m_cal_simple

m_calib, n_calib = fit_calibration(μ, th228_lines)

gr(size=(1200, 1000))
scatter(μ, th228_lines, yerror=μ_err, ms=5, color=:black, framestyle=:box, markershape= :x, layout = @layout[grid(2, 1, heights=[0.8, 0.2])], link=:x, label="Peak Positions", xlabel="Energy (ADC)", xlabelfontsize=10, ylabel="Energy (keV)", ylabelfontsize=10, legend=:topleft, legendfontsize=8, legendfont=font(8), legendtitlefontsize=8, legendtitlefont=font(8), xlims = (0, 5000), xticks = (0:500:12000), margin=2mm, thickness_scaling=1.5, xformatter=:plain)
plot!(ylims = (0, 3000), yticks = (200:250:3000), subplot=1, xlabel="", xticks = :none, bottom_margin=-4mm)
plot!(0:1:20000, x -> m_calib* x + n_calib, label="Best Fit", line_width=2, color=:red, subplot=1, xformatter=_->"")
plot!(μ, ((m_calib .* μ .+ n_calib) .- th228_lines) ./ th228_lines .* 100 , label="", ylabel="Residuals (%)", line_width=2, color=:black, st=:scatter, ylims = (-0.21, 0.21), markershape=:o, subplot=2, legend=:topleft, top_margin=0mm, framestyle=:box)
plot!(legend = :topleft, title="TUM/Polimi ASIC ("* L"R_F"*" external) - FC Calibration Curve", subplot=1, titlefontvalign=:bottom, thickness_scaling=1.8)
savefig(joinpath(figures_folder, "calibration_curve_FC.png"))

th228_names_fit_fwhm = Symbol.(["Tl208a", "Bi212a", "Tl208b", "Bi212FEP", "Tl208SEP", "Tl208FEP"])
th228_names_fit_fwhm = Symbol.(["Tl208a", "Bi212a", "Tl208c", "Bi212d", "Tl208b", "Bi212c", "Bi212FEP", "Tl208FEP"])
fwhm     = ([result_fit[p].fwhm for p in th228_names] ./ m_cal_simple) .* m_calib
fwhm_err = ([result_fit[p].err.σ * 2.355 for p in th228_names] ./ m_cal_simple) .* m_calib

result_fwhm, report_fwhm = fit_fwhm(getindex.([th228_names_dict], th228_names_fit_fwhm), ([result_fit[p].fwhm for p in th228_names_fit_fwhm] ./ m_cal_simple) .* m_calib)
result_fwhm, report_fwhm = fit_fwhm(th228_lines, ([result_fit[p].fwhm for p in th228_names] ./ m_cal_simple) .* m_calib)

scatter(th228_lines, fwhm, yerror=fwhm_err, ms=5, color=:black, framestyle=:box, markershape=:hline, layout = @layout[grid(2, 1, heights=[0.8, 0.2])], link=:x, label="Peak FWHMs", xlabel="Energy (keV)", xlabelfontsize=10, ylabel="FWHM (keV)", ylabelfontsize=10, legend=:topleft, legendfontsize=8, legendfont=font(8), legendtitlefontsize=8, legendtitlefont=font(8), xlims = (0, 3000), xticks = (convert(Int, 0):300:convert(Int, round(3000, digits=0))), margin=5mm, thickness_scaling=1.5, ylims=(1.1, 4.0), yticks=(1:0.5:5.0))
# plot!(0:0.1:3000, x -> report_fwhm.f_fit(x), label="Best Fit: Sqrt($(round(report_fwhm.v[1], digits=2)) + x*$(round(report_fwhm.v[2]*100, digits=2))e-3)", line_width=2, color=:red, subplot=1, xlabel="", xticks=:none, bottom_margin=-4mm)
plot!(0:0.1:3000, x -> report_fwhm.f_fit(x), label="Best Fit", line_width=4, color=:red, subplot=1, xlabel="", xticks=:none, bottom_margin=-4mm)
hline!([result_fwhm.qbb], label="Qbb/keV: $(round(result_fwhm.qbb, digits=2))+-$(round(result_fwhm.err.qbb, digits=2))", color=:green)
hspan!([result_fwhm.qbb - result_fwhm.err.qbb, result_fwhm.qbb + result_fwhm.err.qbb], color=:green, alpha=0.2, label="")
plot!(th228_lines, ((report_fwhm.f_fit.(th228_lines) .- fwhm) ./ fwhm) .* 100 , label="", ylabel="Residuals (%)", line_width=2, color=:black, st=:scatter, ylims = (-5, 5), markershape=:x, legend=:topleft, subplot=2, framestyle=:box, top_margin=0mm, yticks=(-6:2:6))
plot!(legend = :topleft, title=ftitle="TUM/Polimi ASIC ("* L"R_F"*" external) - FC FWHM Curve", subplot=1, margin=5mm, thickness_scaling=2.0)
savefig(joinpath(figures_folder, "resolution_curve_MC.png"))

m_calib, n_calib = 0.5844066009240668, -0.8235788185329779
stephist(data.e_trap .* m_calib .+ n_calib, bins=100:1.0:3000, xlabel="Energy (keV)", ylabel="Counts / keV", yscale=:log10, size=(1500, 700), thickness_scaling=2, label="Energy")
plot!(title="TUM/Polimi ASIC ("* L"R_F"*" external) - Th-228 FC Energy Spectrum", titlefontvalign=:bottom, ylims=(1, 1.5e4), margin=2mm, framestyle=:box, legend=:topright, dpi=300, font=font(20,"Computer Modern"))
ylims!(1, 5e5)
xaxis!(minorticks=6, minorgrid=true, xticks=(0:300:3000), xlims=(0, 3100))
savefig(joinpath(figures_folder, "FC_e_spectrum.png"))