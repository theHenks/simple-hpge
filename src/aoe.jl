using TypedTables, Dates
using Unitful, Formatting, LaTeXStrings, Measures
using Plots, StatsBase, PropDicts, LsqFit
using LegendHDF5IO, LegendDSP, LegendSpecFits
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using ArraysOfArrays
using Distributed
using RadiationDetectorDSP
ENV["JULIA_DEBUG"] = Main # enable debug

plotlyjs(size=(700, 400))

datafolder = "/mnt/artemis02/projects/LEGEND/CUBE3/ASIC-HPGe/"
dsp_folder = joinpath(datafolder, "run-202311-011", "dsp_shortSG")

data = fast_flatten([
            LHDataStore(
                ds -> begin
                    # @debug "Reading from \"$(ds.data_store.filename)\""
                    ds["dsp"][:]
                end,
                joinpath(dsp_folder, fk)
            ) for fk in readdir(dsp_folder) if occursin(".lh5", fk)
        ])

m_calib, n_calib = 0.5844066009240668, -0.8235788185329779

psd_config = readprops("config/run-202311-011/psd_config.json").default
compton_bands  = Vector{Float64}(psd_config.compton_bands)
compton_window = psd_config.compton_window
p_value        = psd_config.p_value

a = data.a
e = data.e_trap .* m_calib .+ n_calib

aoe = a ./ e
qc_cut = isfinite.(aoe) .&& isfinite.(e) .&& isfinite.(a)

aoe, e, a = aoe[qc_cut], e[qc_cut], a[qc_cut]

gr(size=(800, 500))
histogram2d(e, aoe, nbins=(0:0.5:3000, 0.005:1e-4:0.1), xlims=(0, 3000), ylims=(0.04, 0.1), size=(1000, 600), color=cgrad(:magma), colorbar_scale=:log10, legend=:topleft, xlabel="Energy (keV)", ylabel="A/E (a.u.)", margin=5mm)
xticks!(0:250:3000)
plot!(title="TUM/Polimi ASIC ("* L"R_F"*" external) - Th-228 A/E Spectrum")
plot!(margin=1mm, thickness_scaling=1.6, dpi=600, size=(1200, 800))

compton_band_peakhists = generate_aoe_compton_bands(aoe, e, compton_bands, compton_window)

result_fit, report_fit = fit_aoe_compton(compton_band_peakhists.peakhists, compton_band_peakhists.peakstats, compton_bands,; uncertainty=true)

plot(report_fit[1290.0], legend=:topleft)

p_values = [result_fit[band].p_value for band in compton_bands]
μ = [result_fit[band].μ for band in compton_bands]
μ_err = [result_fit[band].err.μ for band in compton_bands]
σ = [result_fit[band].σ for band in compton_bands]
σ_err = [result_fit[band].err.σ for band in compton_bands]

aoe_corrections = fit_aoe_corrections(compton_bands, μ, σ)


scatter(aoe_corrections.e, aoe_corrections.μ, yerr=μ_err, ms=5, color=:black, layout = @layout[grid(2, 1, heights=[0.8, 0.2])], xlabel=L"Energy\ (keV)", ylabel=L"\mu_{A/E} (a.u.)", label=L"\mu_{SCS}", xticks = (minimum(compton_bands):200:maximum(compton_bands)), xlims=(minimum(compton_bands)-50, maximum(compton_bands)+50))
plot!(ylims=(0.9*median(aoe_corrections.μ), 1.05*median(aoe_corrections.μ)), subplot=1, xlabel="", xticks = :none)
plot!(0.0:1500:3000, x -> aoe_corrections.f_μ_scs(x), label="Best Fit", line_width=3.5, color=:red, subplot=1, xformatter=_->"")
plot!(aoe_corrections.e, (aoe_corrections.f_μ_scs.(aoe_corrections.e) .- aoe_corrections.μ) ./ aoe_corrections.μ .* 100 , label="", ylabel=L"Residuals (\%)", line_width=2, color=:black, st=:scatter, ylims = (-5.0, 5.4), markershape=:x, subplot=2)
plot!(legend = :topright, title="TUM/Polimi ASIC ("* L"R_F"*" external) - Th-228 A/E "*L"\mu", subplot=1)
plot!(margin=-1mm, thickness_scaling=1.6, dpi=600, size=(1200, 700))


f_aoe_sigma(x, p) = p[1] .+ p[2]*exp.(-p[3]./x)
σ_scs = curve_fit(f_aoe_sigma, compton_bands, σ, [median(σ)^2, 1, 1])

f_σ_scs = x -> Base.Fix2(f_aoe_sigma, σ_scs.param)(x)

scatter(aoe_corrections.e, aoe_corrections.σ, yerr=σ_err, ms=5, color=:black, layout = @layout[grid(2, 1, heights=[0.8, 0.2])], xlabel=L"Energy\ (keV)", ylabel=L"\sigma_{A/E} (a.u.)", label=L"\sigma_{SCS}", xticks = (minimum(compton_bands):200:maximum(compton_bands)), xlims = (minimum(compton_bands)-50, maximum(compton_bands)+50))
plot!(ylims=(0.1*f_σ_scs(maximum(compton_bands)), 2*f_σ_scs(minimum(compton_bands))), subplot=1, xlabel="", xticks = :none)
plot!(minimum(compton_bands)-50:0.1:maximum(compton_bands)+50, x -> f_σ_scs(x), label=format("Best Fit: sqrt({:.2E}+({:.2E}/x^2)", aoe_corrections.σ_scs...), line_width=3.5, color=:red, subplot=1, xformatter=_->"")
plot!(aoe_corrections.e, (f_σ_scs.(aoe_corrections.e) .- aoe_corrections.σ) ./ aoe_corrections.σ .* 100 , label="Residuals", ylabel="Residuals (%)", line_width=2, color=:black, st=:scatter, ylims = (-50.0, 50.0), markershape=:x, subplot=2)
plot!(legend = :topright, title="TUM/Polimi ASIC ("* L"R_F"*" external) - Th-228 A/E "*L"\sigma", subplot=1)

aoe_norm = aoe .- Base.Fix2(LegendSpecFits.f_aoe_mu, aoe_corrections.μ_scs).(e)
aoe_norm = aoe_norm ./ f_σ_scs.(e) 


gr()
histogram2d(e, aoe_norm, nbins=(0:0.5:3000, -20:0.05:20), xlims=(0, 3000), ylims=(-20, 20), size=(1000, 600), color=cgrad(:magma), colorbar_scale=:log10, legend=:topleft, xlabel=L"Energy\ (keV)", ylabel=L"A/E\ (\sigma_{A/E})")
plot!(margin=1mm, thickness_scaling=1.6, dpi=600, size=(1500, 800))
plot!(title="TUM/Polimi ASIC ("* L"R_F"*" external) - Th-228 Normalized A/E Spectrum")
savefig("AoE_classifier.png")


psd_peaks = Float64.(psd_config.psd_peaks)
psd_peaks_window_sizes = Vector{Tuple{Float64, Float64}}([(l,r) for (l,r) in zip(Vector{Float64}(psd_config.psd_peaks_windows_left), Vector{Float64}(psd_config.psd_peaks_windows_right))])
psd_peak_names = Symbol.(psd_config.psd_peaks_names)
psd_peak_dict = Dict(zip(psd_peak_names, psd_peaks))
qbb =  psd_config.qbb
qbb_window = psd_config.qbb_window
sigma_high_sided = 5.0
e_type = Symbol(psd_config.energy_type)

# plotlyjs()
# dep=1592.53
# window=[12.0, 10.0]
# stephist(aoe[dep-first(window).< e .< dep+last(window)], bins=0.06:1e-4:0.07)
# xlims!(0.06, 0.07)
# band = compton_bands[1]
# stephist!(aoe[band.< e .< band+compton_window], bins=0.06:1e-4:0.07)

# stephist(aoe_norm[dep-first(window).< e .< dep+last(window)], bins=-250:0.5:50)
# xlims!(0.06, 0.07)
# band = compton_bands[1]
# stephist(aoe_norm[band.< e .< band+compton_window],bins=-100.0:5e-1:50.0)
# plot(report_fit[band], legend=:topleft)

result_cut = get_psd_cut(aoe_norm, e,; cut_search_interval=(-10.0, 0.0), window=[20.0, 20.0], rtol=1e-5, bin_width_window=5.0, fixed_position=false)            

result_peaks, report_peaks = get_peaks_surrival_fractions(aoe_norm, e, psd_peaks, psd_peak_names, psd_peaks_window_sizes, result_cut.cut,; bin_width_window=3.0, low_e_tail=false, sigma_high_sided=3.0)

peak_sf_plot = plot.([rep.after for rep in values(report_peaks)], titleloc=:left, titlefont=font(8), ticks=:native, legend=:bottomright; show_label=true, show_fit=true)
for (p, rep_before) in zip(peak_sf_plot, [rep.before for rep in values(report_peaks)])
    plot!(p, rep_before,; show_label=true, show_fit=false)
    p.series_list[1][:label] = "After"
    p.series_list[2][:label] = "Before"
end
for (p, peak_name, res) in zip(peak_sf_plot, keys(result_peaks), values(result_peaks))
    xticks!(p, convert(Int, round(xlims(p)[1], digits=0)):10:convert(Int, round(xlims(p)[2], digits=0)))
    title!(p, format("{} ({} keV) - SF: {:.2f} ± {:.2f}%", string(peak_name), psd_peak_dict[peak_name], res.sf*100, res.err.sf*100))
end
plot(
    peak_sf_plot...,
    layout = @layout[grid(2, 2)], 
    size=(2000, 1200), legend=:bottomright,
    framestyle=:box,
    grid=true, minor=true, gridalpha=0.2, gridcolor=:black, gridlinewidth=0.5,
    xlabel="Energy (keV)", ylabel="Counts",
    dpi = 300, thickness_scaling = 2,
    yformatter=:plain, titlefont=12,
    fontfamily=:sansserif
)
plot!(margin=1mm, thickness_scaling=1.2, dpi=600, size=(1200, 900))
plot!(plot_title="TUM/Polimi ASIC - Th-228 Peak SF")

qbb_result = get_continuum_surrival_fraction(aoe_norm, e, qbb, qbb_window, result_cut.cut,; sigma_high_sided=sigma_high_sided)



data_peaks = LHDataStore(peaks_file, "r")

data_sep = data_peaks["Tl208SEP"][:]

plot(u"µs", NoUnits)
plot!(data_sep.waveform[1:10])
xlims!(34, 36)