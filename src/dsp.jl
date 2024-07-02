using TypedTables, Dates
using Unitful, Formatting, LaTeXStrings, Measures, Measurements
using Measurements: value as mvalue, uncertainty as muncertainty
using Plots, StatsBase, PropDicts
using LegendHDF5IO, LegendDSP, LegendSpecFits
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using LegendDSP: get_fltpars
using ArraysOfArrays
using Distributed
using RadiationDetectorDSP
using LegendDataManagement
using LegendDataManagement: readlprops, writelprops
import Pkg

# enable debug
ENV["JULIA_DEBUG"] = "Main,LegendSpecFits,LegendDSP"
# set CPU target to generic
ENV["JULIA_CPU_TARGET"] = "generic"


# You can decide here on the plotting backend
# for interactive plotting use plotljyjs, for static plots use gr
plotlyjs(size=(700, 500))
# gr(size=(700, 500))

# read datafolder and run name from global config
data_pd = readprops("data_config.json")

datafolder = data_pd.datafolder
raw_folder = joinpath(datafolder, data_pd.run, "raw")
dsp_folder = mkpath(joinpath(datafolder, data_pd.run, "dsp_tst"))

config_folder = mkpath(data_pd.configfolder)
pars_folder   = mkpath(data_pd.parsfolder)
log_folder    = mkpath(data_pd.logfolder)

reprocess = true

filenames = readdir(raw_folder)
filenames = filenames[endswith.(filenames, ".lh5")]


dsp_meta = readlprops(joinpath(config_folder, "dsp_config.json"))
dsp_config = DSPConfig(dsp_meta.default)

# get pars
pars_folder = mkpath(joinpath(data_pd.parsfolder, data_pd.run))
pars_file = joinpath(pars_folder, "pars_optimization.json")
pars_db = 
    if isfile(pars_file)
        get_values(readlprops(pars_file))
    else
        PropDict()
    end

# kill all workers at exist
function kill_sessions()
    @info "Kill all sessions"
    # kill all sessions
    rmprocs(workers()...)
    run(`pkill -u $USER -f worker`)
end
atexit(kill_sessions)

# number of workers to run dsp with in paralell
n_workers = 5
# config parallel processing
addprocs(n_workers, exeflags=`--threads=4 --project=$(dirname(Pkg.project().path)) --heap-size-hint=10G`, topology=:master_worker)

# Sanity check to check for correct setup of workers
worker_ids = Distributed.remotecall_fetch.(Ref(Distributed.myid), Distributed.workers())
@assert length(worker_ids) == Distributed.nworkers()

@info "$(Distributed.nworkers()) Julia worker processes active."

# activate all packages on all possible workers
@everywhere begin
    using LegendHDF5IO, LegendDSP, LegendSpecFits, LegendDataTypes
    using IntervalSets, TypedTables, StatsBase, PropDicts
    using Unitful, Formatting, Printf, Measures, Dates
    using Distributed, ProgressMeter, TimerOutputs
    using RadiationDetectorDSP
    using ArraysOfArrays

    using HDF5
    using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
    using LegendDSP: get_fltpars
    using TypedTables, Dates
    using Measurements: value as mvalue, uncertainty as muncertainty
    using LegendDataManagement: readlprops, writelprops
    import Pkg
    # free memory
    GC.gc()
end

@everywhere begin
    raw_folder = $raw_folder
    dsp_folder = $dsp_folder
    dsp_config = $dsp_config
    pars_db = $pars_db
end


###################
# Main DSP function
###################
@everywhere function dsp(data::Q, config::LegendDSP.DSPConfig, τ::Quantity{T}, pars_filter::PropDict) where {Q <: Table, T<:Real}
    # get config parameters
    bl_window                = config.bl_window
    t0_threshold             = config.t0_threshold
    tail_window              = config.tail_window
    inTraceCut_std_threshold = config.inTraceCut_std_threshold
    sg_flt_degree            = config.sg_flt_degree
    current_window           = config.current_window
    qdrift_int_length        = config.qdrift_int_length
    lq_int_length            = config.lq_int_length
    
    # get optimal filter parameters
    trap_rt, trap_ft = get_fltpars(pars_filter, :trap, config)
    cusp_rt, cusp_ft = get_fltpars(pars_filter, :cusp, config)
    zac_rt, zac_ft = get_fltpars(pars_filter, :zac, config)
    sg_wl   = get_fltpars(pars_filter, :sg, config)

    # get waveform data 
    wvfs = data.waveform
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # get CUSP and ZAC filter length and flt scale
    flt_length_zac              = config.flt_length_zac
    zac_scale                   = ustrip(NoUnits, flt_length_zac/step(wvfs[1].time))
    flt_length_cusp             = config.flt_length_cusp
    cusp_scale                  = ustrip(NoUnits, flt_length_cusp/step(wvfs[1].time))

    # get number of samples the waveform is saturated at low and high of FADC range
    bit_depth = config.kwargs_pars.fc_bit_depth # of FlashCam FADC
    sat_low, sat_high = 0, 2^bit_depth - bit_depth
    sat_stats = saturation.(wvfs, sat_low, sat_high)

    # set tau for CUSP filter to very high number to switch of CR filter
    τ_cusp = 10000000.0u"µs"
    τ_zac = 10000000.0u"µs"

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # pretrace difference 
    pretrace_diff = flatview(wvfs.signal)[1, :] - bl_stats.mean

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)
    
    # get wvf maximum
    wvf_max = maximum.(wvfs.signal)
    wvf_min = minimum.(wvfs.signal)

    # extract decay times
    tail_stats = tailstats.(wvfs, leftendpoint(tail_window), rightendpoint(tail_window))

    # deconvolute waveform 
    # --> wvfs = wvfs_pz
    deconv_flt = InvCRFilter(τ)
    wvfs = deconv_flt.(wvfs)

    # get tail mean, std and slope
    pz_stats = signalstats.(wvfs, leftendpoint(tail_window), rightendpoint(tail_window))

    # t0 determination
    t0 = get_t0(wvfs, t0_threshold; flt_pars=config.kwargs_pars.t0_flt_pars, mintot=config.kwargs_pars.t0_mintot)

    # if all waveforms are saturated set threshold to 1.0 to avoid numerical problems
    # replace!(wvf_max, zero(wvf_max[1]) => one(wvf_max[1]))

    # get threshold points in rise
    t10 = get_threshold(wvfs, wvf_max .* 0.1; mintot=config.kwargs_pars.tx_mintot)
    t50 = get_threshold(wvfs, wvf_max .* 0.5; mintot=config.kwargs_pars.tx_mintot)
    t80 = get_threshold(wvfs, wvf_max .* 0.8; mintot=config.kwargs_pars.tx_mintot)
    t90 = get_threshold(wvfs, wvf_max .* 0.9; mintot=config.kwargs_pars.tx_mintot)
    t99 = get_threshold(wvfs, wvf_max .* 0.99; mintot=config.kwargs_pars.tx_mintot)
    
    drift_time = uconvert.(u"ns", t90 - t0)

    # get Q-drift parameter
    qdrift = get_qdrift(wvfs, t0, qdrift_int_length; pol_power=config.kwargs_pars.int_interpolation_order, sign_est_length=config.kwargs_pars.int_interpolation_length)

    # get LQ parameter
    lq  = get_qdrift(wvfs, t80, lq_int_length; pol_power=config.kwargs_pars.int_interpolation_order, sign_est_length=config.kwargs_pars.int_interpolation_length)

    # robust energy reconstruction with long, middle and short rise and flat-top times
    uflt_10410 = TrapezoidalChargeFilter(10u"µs", 4u"µs")
    e_10410  = maximum.((uflt_10410.(wvfs)).signal)

    uflt_535 = TrapezoidalChargeFilter(5u"µs", 3u"µs")
    e_535  = maximum.((uflt_535.(wvfs)).signal)

    uflt_313 = TrapezoidalChargeFilter(3u"µs", 1u"µs")
    e_313  = maximum.((uflt_313.(wvfs)).signal)

    # signal estimator for precise energy reconstruction
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))

    # get trap energy of optimized rise and flat-top time
    uflt_trap_rtft = TrapezoidalChargeFilter(trap_rt, trap_ft)

    e_trap = signal_estimator.(uflt_trap_rtft.(wvfs), t50 .+ (trap_rt + trap_ft/2))

    # get cusp energy of optimized rise and flat-top time
    uflt_cusp_rtft = CUSPChargeFilter(cusp_rt, cusp_ft, τ_cusp, flt_length_cusp, cusp_scale)

    e_cusp = signal_estimator.(uflt_cusp_rtft.(wvfs), t50 .+ (flt_length_cusp /2))

    # get zac energy of optimized rise and flat-top time
    uflt_zac_rtft = ZACChargeFilter(zac_rt, zac_ft, τ_zac, flt_length_zac, zac_scale)

    e_zac = signal_estimator.(uflt_zac_rtft.(wvfs), t50 .+ (flt_length_zac /2))

    # extract current with optimal SG filter length with second order polynominal and first derivative
    wvfs_sgflt_deriv = SavitzkyGolayFilter(sg_wl, sg_flt_degree, 1).(wvfs)
    a_sg = get_wvf_maximum.(wvfs_sgflt_deriv, leftendpoint(current_window), rightendpoint(current_window))

    a_60 = get_wvf_maximum.(SavitzkyGolayFilter(60u"ns", sg_flt_degree, 1).(wvfs), leftendpoint(current_window), rightendpoint(current_window))
    a_100 = get_wvf_maximum.(SavitzkyGolayFilter(100u"ns", sg_flt_degree, 1).(wvfs), leftendpoint(current_window), rightendpoint(current_window))
    a_raw = get_wvf_maximum.(DifferentiatorFilter(1).(wvfs), leftendpoint(current_window), rightendpoint(current_window))

    # get in-trace pile-up
    inTrace_pileUp = get_intracePileUp(wvfs_sgflt_deriv, inTraceCut_std_threshold, bl_window; mintot=config.kwargs_pars.intrace_mintot)
    
    # get position of current rise
    thres = maximum.(wvfs_sgflt_deriv.signal) .* 0.5
    # replace!(thres, zero(thres[1]) => one(thres[1]))

    t50_current = get_threshold(wvfs_sgflt_deriv, thres; mintot=config.kwargs_pars.tx_mintot)

    # invert waveform for DC tagging
    # wvfs --> wvfs_pz_inv
    wvfs = multiply_waveform.(wvfs, -1.0)

    # get inverted waveform maximum for long and short filter
    e_10410_max_inv  = maximum.(uflt_10410.(wvfs).signal)

    e_313_max_inv  = maximum.(uflt_313.(wvfs).signal)

    # t0 determination
    t0_inv = get_t0(wvfs, t0_threshold; mintot=config.kwargs_pars.t0_mintot)

    # output Table 
    return TypedTables.Table(blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
        tailmean = pz_stats.mean, tailsigma = pz_stats.sigma, tailslope = pz_stats.slope, tailoffset = pz_stats.offset,
        t0 = t0, t10 = t10, t50 = t50, t80 = t80, t90 = t90, t99 = t99,
        t50_current = t50_current, 
        drift_time = drift_time,
        tail_τ = tail_stats.τ, tail_mean = tail_stats.mean, tail_sigma = tail_stats.sigma,
        e_max = wvf_max, e_min = wvf_min,
        e_10410 = e_10410, e_535 = e_535, e_313 = e_313,
        e_10410_inv = e_10410_max_inv, e_313_inv = e_313_max_inv,
        t0_inv = t0_inv,
        e_trap = e_trap, e_cusp = e_cusp, e_zac = e_zac, 
        qdrift = qdrift, lq = lq,
        a_sg = a_sg, a_60 = a_60, a_100 = a_100, a_raw = a_raw,
        blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
        pretrace_diff = pretrace_diff, 
        inTrace_intersect = inTrace_pileUp.intersect, inTrace_n = inTrace_pileUp.n,
        n_sat_low = sat_stats.low, n_sat_high = sat_stats.high, n_sat_low_cons = sat_stats.max_cons_low, n_sat_high_cons = sat_stats.max_cons_high
    )
end

@everywhere function single_file_dsp(idx::Int64)
    dsp_timer = TimerOutput()
    filename    = joinpath(raw_folder, filenames[idx])
    outfilename = joinpath(dsp_folder, replace(filenames[idx], ".lh5" => "_dsp.lh5"))
    @timeit dsp_timer "DSP $(basename(filename))" begin

        @info "Processing file: $(basename(filename))"
        data    = lh5open(filename, "r")
        @info "Using output file: $(basename(outfilename))"
        # if reprocess && isfile(outfilename)
        #     @info "Reprocess $(basename(outfilename)), remove old DSP."
        #     rm(outfilename)
        # else
        #     try 
        #         lh5open(outfilename, "cw")
        #     catch e
        #         @warn "LoadError: $e"
        #         @warn "Filename $(basename(outfilename)) seems broken, remove it."
        #         rm(outfilename)
        #     end
        # end
        outdata = lh5open(outfilename, "cw")
    end

    @info "Start DSP"
    @timeit dsp_timer "DSP $(basename(filename))" begin

        # load data from HDF5
        data_ch = data["FCEvent"][:]
        outdata["dsp"]  = dsp(data_ch, dsp_config, pars_db.pz.tau, pars_db)
        # free memory
        GC.gc()
    end
    @info "Finished processing file: $(basename(filename))"

    close(data)
    close(outdata)

    total_time      = canonicalize(Dates.Nanosecond(TimerOutputs.tottime(dsp_timer)))
    total_allocated = Base.format_bytes(TimerOutputs.totallocated(dsp_timer))
    @info "Total time: $total_time, Total allocated: $total_allocated"
    return (timer = dsp_timer, )
end

result_dsp = @showprogress pmap(eachindex(filenames), batch_size = 1, on_error=identity) do idx
    try
        idx => single_file_dsp(idx)
    catch e
        @error "Error in file: $(filenames[idx]): $e"
        idx => (timer = TimerOutput(), )
    end
end

rmprocs(workers()...)

total_dsp_timer = TimerOutput()
for (idx, res) in result_dsp
    # merge timer into total timer
    merge!(total_dsp_timer, res.timer)
end
@info "Finished processing all files."

main_log = """
# Total Timing
```
$total_dsp_timer
```
"""
@info main_log
log_filename = joinpath(log_folder, "main_log.log")
open(log_filename, "w+") do file
    write(file, main_log)
end

@info "Finished DSP processing."