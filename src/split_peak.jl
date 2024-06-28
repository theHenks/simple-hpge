using TypedTables, Dates
using Unitful, Formatting, LaTeXStrings, Measures, IntervalSets
using Plots, StatsBase, PropDicts
using LegendHDF5IO, LegendDSP, LegendSpecFits
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using LegendDataTypes, HDF5

# You can decide here on the plotting backend
# for interactive plotting use plotljyjs, for static plots use gr
plotlyjs(size=(700, 500))
# gr(size=(700, 500))

# read datafolder and run name from global config
data_pd = readprops("data_config.json")

datafolder = data_pd.datafolder
raw_folder = joinpath(datafolder, data_pd.run, "raw")

# define energy windows around peaks of interest
energy_windows = IdDict(
        :Tl208a => 558u"keV"..608u"keV",
        :Bi212a => 702u"keV"..752u"keV",
        :Tl208b => 836u"keV"..886u"keV",
        :Tl208DEP_Bi212FEP => 1568u"keV"..1646u"keV",
        :Tl208SEP => 2079u"keV"..2129u"keV",
        :Tl208FEP => 2590u"keV"..2640u"keV"
)

# load raw DAQ online energy from FADC into memore by flattening over all run files
E_raw = fast_flatten([
    lh5open(
        ds -> begin
            @info "Reading DAQ energy from \"$(ds.data_store.filename)\""
            ds["FCEvent/daqenergy"][:]
        end,
        joinpath(raw_folder, fk)
    ) for fk in readdir(raw_folder) if occursin(".lh5", fk)
])

# plot raw ADC energy histogram
stephist(E_raw, bins=10000, yscale=:log10)

# set max energy in ADC
# in case you expect events with higher energies than the Tl-208 FEP you should set a max ADC here
e_max = maximum(E_raw)
e_max = 7000
f_calib, diagnostics = autocal_energy(E_raw[E_raw .< e_max])
e_cal_simple = f_calib(E_raw)

### You should check here if the simple calibrated histogram looks good
stephist(e_cal_simple[e_cal_simple .< 4000u"keV"], bins=10000, yscale=:log10)

# split original .lh5 data into energy windows
slim_data = flatten_by_key([LHDataStore(joinpath(raw_folder,fk)) do ds
    @info "Filtering $(fk)"
    filter_raw_data_by_energy(ds["FCEvent"][:], f_calib, energy_windows)
end for fk in readdir(raw_folder) if occursin(".lh5", fk)])

# create folder for peak file
peaks_folder = joinpath(datafolder, data_pd.run, "peaks")
mkpath(peaks_folder)

# save peak files
output_filename = joinpath(peaks_folder, format("{}_{}_{}.lh5", split(readdir(raw_folder)[1], "_")[1], split(readdir(raw_folder)[1], "_")[2], "peaks"))
h5open(output_filename, "w") do output
    for label in sort(collect(keys(slim_data)))
        LegendDataTypes.writedata(output, "$label", slim_data[label])
    end
end

# Sanity check: Load FEP waveforms and plot them
data_fep = lh5open(output_filename, "r")["Tl208FEP"][:]

wvfs_fep = data_fep.waveform
ts = unix2datetime.(ustrip.(data_fep.timestamp))

# select index range to be plotted
idx = 1:10
plot(wvfs_fep[idx], label=permutedims(ts[idx]), legend=:bottomright, xlabel="Time", ylabel="Singal (ADC)", title="FEP waveforms", xunit=u"µs")