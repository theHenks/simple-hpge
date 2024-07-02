using TypedTables, Dates
using Unitful, Formatting, LaTeXStrings, Measures, Measurements
using Measurements: value as mvalue, uncertainty as muncertainty
using Plots, StatsBase, PropDicts
using LegendHDF5IO, LegendDSP, LegendSpecFits
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using ArraysOfArrays
using Distributed
using RadiationDetectorDSP
using LegendDataManagement: readlprops, writelprops


# enable debug
ENV["JULIA_DEBUG"] = "Main,LegendSpecFits,LegendDSP" 


# You can decide here on the plotting backend
# for interactive plotting use plotljyjs, for static plots use gr
plotlyjs(size=(700, 500))
# gr(size=(700, 500))

# read datafolder and run name from global config
data_pd = readprops("data_config.json")

datafolder = data_pd.datafolder
raw_folder = joinpath(datafolder, data_pd.run, "raw")
config_folder = data_pd.configfolder
figures_folder = mkpath(joinpath(data_pd.figurefolder, data_pd.run))
# get pars
pars_folder = mkpath(joinpath(data_pd.parsfolder, data_pd.run))
pars_file = joinpath(pars_folder, "pars_optimization.json")
pars_db = 
    if isfile(pars_file)
        readlprops(pars_file)
    else
        PropDict()
    end

# get decay time
tau = mvalue(pars_db.pz.τ)

# get peak files
peaks_folder = joinpath(datafolder, data_pd.run, "peaks")

peaks_file = joinpath(peaks_folder, only(readdir(peaks_folder)))

# open peaks file to load waveforms from specfic peaks
wvfs_fep = lh5open(peaks_file, "r")["Tl208FEP"].waveform[:]
plot(wvfs_fep[1:20], label="Tl208FEP", legend=:outertopright, xunit=u"µs")

# get DSP configuration data --> Can be modified in .json file
dsp_meta = readlprops(joinpath(config_folder, "dsp_config.json"))
dsp_config = DSPConfig(dsp_meta.default)

# config for optimization procedure
# in principle, different energy filters can be used for this step. 
# --> At the moment, only the trap, cusp and zac filter are implemented
opti_config = dsp_meta.optimization.default.trap

# minima and maximal baseline ENC noise values for the histogram in ADC
min_enc, max_enc = opti_config.min_enc, opti_config.max_enc

# simple DSP to extract ENC baseline values for different rise time combinations
enc_trap_grid = dsp_trap_rt_optimization(wvfs_fep, dsp_config, tau,; ft=1.0u"µs")

i = 11
enc_i = flatview(enc_trap_grid)[i,:]
histogram(enc_i, bins=opti_config.nbins_enc_sigmas, label="ENC Distribution", xlabel="ENC [e-]", title="RT: $(dsp_config.e_grid_rt_trap[i])", legend=:outertopright, xlims=(min_enc, max_enc))

# fit the ENC distribution with a truncated gaussian to get optimal RT value
result_rt, report_rt = fit_enc_sigmas(enc_trap_grid, dsp_config.e_grid_rt_trap, opti_config.min_enc, opti_config.max_enc, opti_config.nbins_enc_sigmas, opti_config.rel_cut_fit_enc_sigmas)

# plot resulting ENC rt dependency
plot(report_rt)
savefig(joinpath(figures_folder, "rt_optimization.png"))

@info "Optimal RT value: $(round(u"µs", result_rt.rt, digits=2))"

# simple DSP to extract FWHM values for different rise time combinations
e_trap_grid = dsp_trap_ft_optimization(wvfs_fep, dsp_config, tau, mvalue(result_rt.rt))

# fit the FWHM values with a gamma peakshape to get optimal FT value
result_ft, report_ft = fit_fwhm_ft_fep(e_trap_grid, dsp_config.e_grid_ft_trap, mvalue(result_rt.rt), opti_config.min_e_fep, opti_config.max_e_fep, opti_config.nbins_e_fep, opti_config.rel_cut_fit_e_fep)

# plot resulting FWHM ft dependency
plot(report_ft, ylims=(2, 4))
savefig(joinpath(figures_folder, "ft_fwhm_fit.png"))

@info "Optimal FT value: $(round(u"µs", result_ft.ft, digits=2))"

# add pars to pars database while overwriting existing pars
pars_db.trap = merge(result_rt, result_ft)

# write pars to file
writelprops(pars_file, pars_db)