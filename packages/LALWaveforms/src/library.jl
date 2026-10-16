# Locating, loading and caching symbols from liblal and liblalsimulation.
#
# Loading is lazy (first use) because `__init__` must never throw: the package has
# to load, and pass Aqua, on machines without LAL.

const LIBLAL_SYMBOLS = (
    :XLALCreateDict,
    :XLALDestroyDict,
    :XLALDestroyValue,
    :XLALDestroyCOMPLEX16FrequencySeries,
    :XLALGetErrnoPtr,
    :XLALGetBaseErrno,
    :XLALClearErrno,
    :XLALErrorString,
    :XLALSetSilentErrorHandler,
    :XLALClobberDebugLevel,
)

const LIBSIM_SYMBOLS = (
    :XLALSimInspiralChooseFDWaveform,
    :XLALSimInspiralFD,
    :XLALSimInspiralImplementedFDApproximants,
    :XLALSimInspiralGetApproximantFromString,
    :XLALSimInspiralTransformPrecessingNewInitialConditions,
    :XLALSimInspiralfLow2fStart,
    :XLALSimInspiralWaveformParamsInsertPNSpinOrder,
    :XLALSimInspiralWaveformParamsInsertPNTidalOrder,
    :XLALSimInspiralWaveformParamsInsertPNPhaseOrder,
    :XLALSimInspiralWaveformParamsInsertPNAmplitudeOrder,
    :XLALSimInspiralWaveformParamsInsertTidalLambda1,
    :XLALSimInspiralWaveformParamsInsertTidalLambda2,
    :XLALSimInspiralWaveformParamsInsertModeArray,
    :XLALSimInspiralCreateModeArray,
    :XLALSimInspiralModeArrayActivateMode,
)

# One `Ptr{Cvoid}` field per symbol, so `lal().XLALFoo` is a type-stable field read.
@eval struct LALFunctions
    $((:($s::Ptr{Cvoid}) for s in (LIBLAL_SYMBOLS..., LIBSIM_SYMBOLS...))...)
end

const STATE = Ref{Union{Nothing, LALFunctions}}(nothing)
const LOADED_FROM = Ref{String}("")
const LOAD_LOCK = ReentrantLock()

const LIBDIR_ENV = "LALWAVEFORMS_LIBDIR"
const LIBDIR_PREF = "libdir"

"""
    LALError(func, errno, msg)

Error raised for a failed LAL call (`func` is the XLAL function, `errno` its base
errno) or for a problem locating the LAL libraries (`func == :libdir`).
"""
struct LALError <: Exception
    func::Symbol
    errno::Int
    msg::String
end

function Base.showerror(io::IO, e::LALError)
    print(io, "LALError: ", e.func)
    e.errno != 0 && print(io, " (errno ", e.errno, ")")
    return print(io, ": ", e.msg)
end

const SETUP_HELP = """
LALWaveforms needs the LAL shared libraries (liblal, liblalsimulation) but no
library directory is configured. Either
  * set the environment variable $LIBDIR_ENV to the directory containing them, or
  * call `LALWaveforms.set_libdir!(path)` (saved as a preference), or
    `LALWaveforms.set_libdir_from_python!()` to take them from a `lalsuite` wheel
    installed in the current Python environment."""

# `liblal.20.dylib`, `liblal-9eb973fe.so.20.5.1`; never `liblalsimulation*` or `liblalburst*`.
const LIBLAL_REGEX = r"^liblal([.-].*)?\.(dylib|so)(\.\d+)*$"
const LIBSIM_REGEX = r"^liblalsimulation([.-].*)?\.(dylib|so)(\.\d+)*$"

function find_library(dir::AbstractString, rx::Regex)
    matches = filter(name -> occursin(rx, name), readdir(dir))
    isempty(matches) && return nothing
    return joinpath(dir, last(sort!(matches)))
end

"""
    find_libraries(dir) -> (liblal, liblalsimulation)

Paths of the two libraries inside `dir`, or `nothing` if either is missing. Handles
the macOS wheel layout (`lal/.dylibs/`) and the Linux auditwheel layout
(`lalsuite.libs/`, hash-suffixed names).
"""
function find_libraries(dir::AbstractString)
    isdir(dir) || return nothing
    lal = find_library(dir, LIBLAL_REGEX)
    sim = find_library(dir, LIBSIM_REGEX)
    return (lal === nothing || sim === nothing) ? nothing : (lal, sim)
end

function resolve_libdir()
    dir = get(ENV, LIBDIR_ENV, "")
    isempty(dir) || return dir
    pref = @load_preference(LIBDIR_PREF)
    pref === nothing || return String(pref)
    throw(LALError(:libdir, 0, SETUP_HELP))
end

function load_library!(dir::AbstractString)
    paths = find_libraries(dir)
    paths === nothing && throw(
        LALError(:libdir, 0, "no liblal / liblalsimulation found in \"$dir\"; $SETUP_HELP")
    )
    # liblalsimulation resolves its liblal dependency through the global namespace.
    liblal = dlopen(paths[1], RTLD_GLOBAL)
    libsim = dlopen(paths[2], RTLD_GLOBAL)
    ptrs = (
        (dlsym(liblal, s) for s in LIBLAL_SYMBOLS)...,
        (dlsym(libsim, s) for s in LIBSIM_SYMBOLS)...,
    )
    fns = LALFunctions(ptrs...)
    STATE[] = fns
    LOADED_FROM[] = String(dir)
    silent_errors!(fns)
    return fns
end

"""
    lal() -> LALFunctions

Cached LAL symbol table, loading the libraries on first use.
"""
function lal()::LALFunctions
    s = STATE[]
    s === nothing || return s
    return lock(LOAD_LOCK) do
        s = STATE[]
        s === nothing ? load_library!(resolve_libdir()) : s
    end
end

# XLAL error handlers are per-thread, so this is re-applied at each entry point. The
# silent handler only stops LAL from aborting; the "XLAL Error - ..." line printed at
# the failure site is controlled by the debug level, so that is zeroed as well. Both
# are process-wide LAL settings that other LAL users in the same process will see.
function silent_errors!(fns::LALFunctions = lal())
    ccall(fns.XLALSetSilentErrorHandler, Ptr{Cvoid}, ())
    ccall(fns.XLALClobberDebugLevel, Cvoid, (Cint,), 0)
    return nothing
end

"""
    set_libdir!(path)

Save `path` (a directory holding `liblal` and `liblalsimulation`) as the `libdir`
preference and load the libraries from it. The `LALWAVEFORMS_LIBDIR` environment
variable, when set, still takes precedence in later sessions. Libraries cannot be
unloaded, so changing the directory after a successful load only takes effect in
the next Julia session.
"""
function set_libdir!(path::AbstractString)
    dir = abspath(path)
    find_libraries(dir) === nothing && throw(
        ArgumentError("no liblal / liblalsimulation found in \"$dir\"")
    )
    @set_preferences!(LIBDIR_PREF => dir)
    lock(LOAD_LOCK) do
        if STATE[] === nothing
            load_library!(dir)
        elseif LOADED_FROM[] != dir
            @warn "LAL already loaded from $(LOADED_FROM[]); $dir will be used after restarting Julia"
        end
    end
    return dir
end

# Keep in sync with test/python/libdir.py.
const PYTHON_LIBDIR_SNIPPET = """
import glob, os, sys
import lal, lalsimulation
site = os.path.dirname(os.path.dirname(lal.__file__))
for d in (os.path.join(site, "lal", ".dylibs"), os.path.join(site, "lalsuite.libs"),
          os.path.join(sys.prefix, "lib")):
    if glob.glob(os.path.join(d, "liblalsimulation*")):
        print(d)
        break
else:
    sys.exit("could not find liblalsimulation next to the lal Python package")
"""

"""
    set_libdir_from_python!(python = "python3")

Find the LAL libraries bundled with (or installed next to) the `lalsuite` Python
package visible to `python`, then call [`set_libdir!`](@ref) on that directory.
"""
function set_libdir_from_python!(python::AbstractString = "python3")
    dir = strip(read(`$python -c $PYTHON_LIBDIR_SNIPPET`, String))
    return set_libdir!(dir)
end
