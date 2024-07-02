using TypedTables, Dates
using Unitful, Formatting, LaTeXStrings, Measurements, Measures
using Plots, StatsBase, PropDicts
using LegendHDF5IO, LegendDSP, LegendSpecFits
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using LegendDataManagement: readlprops, writelprops


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

peaks_folder = joinpath(datafolder, data_pd.run, "peaks")

peaks_file = joinpath(peaks_folder, only(readdir(peaks_folder)))

# open peaks file to load waveforms from specfic peaks
wvfs_fep = lh5open(peaks_file, "r")["Tl208FEP"].waveform[:]
plot(wvfs_fep[1:20], label="Tl208FEP", legend=:outertopright, xunit=u"µs")

# get DSP configuration data --> Can be modified in .json file
dsp_meta = readlprops(joinpath(config_folder, "dsp_config.json"))
dsp_config = DSPConfig(dsp_meta.default)

decay_times = dsp_decay_times(wvfs_fep, dsp_config)


# get configuration for decay time extraction
pz_config = dsp_meta.pz.default
min_τ, max_τ = pz_config.min_tau, pz_config.max_tau
nbins = pz_config.nbins
rel_cut_fit = pz_config.rel_cut_fit

# plot decay time distribution
histogram(decay_times[min_τ .< decay_times .< max_τ], bins=:fd, label="Decay Time Distribution", xlabel="Decay Time")
savefig(joinpath(figures_folder, "decay_times.png"))

# at first define cut around peak top to get a better fit
cuts_τ = cut_single_peak(decay_times, min_τ, max_τ,; n_bins=nbins, relative_cut=rel_cut_fit)

# fit the decay time distribution at the peak top with a truncated gaussian
result, report = fit_single_trunc_gauss(decay_times, cuts_τ)

# plot resulting ditribution
plot(report, decay_times, cuts_τ, xlabel="Decay Time [µs]", title="Decay Time Distribution")
savefig(joinpath(figures_folder, "decay_time_fit.png"))

# add pars to pars database while overwriting existing pars
pars_db.pz.fit = result
pars_db.pz.τ = result.μ

# write pars to file
writelprops(pars_file, pars_db)

