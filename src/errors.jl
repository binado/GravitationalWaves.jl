# XLAL functions report failure by a nonzero return value (or NULL) plus a
# thread-local errno. With the silent handler installed nothing is printed, so we
# translate the errno into a `LALError` ourselves.

clear_errno!(fns::LALFunctions = lal()) = ccall(fns.XLALClearErrno, Cint, ())

"""
    check(rc, func)

Throw a [`LALError`](@ref) if LAL's errno is set or `rc != 0`, clearing the errno.
Must run on the same thread, with no yield in between, as the failing call.
"""
function check(rc::Integer, func::Symbol, fns::LALFunctions = lal())
    errno = unsafe_load(ccall(fns.XLALGetErrnoPtr, Ptr{Cint}, ()))
    (errno == 0 && rc == 0) && return nothing
    base = ccall(fns.XLALGetBaseErrno, Cint, ())
    msg = unsafe_string(ccall(fns.XLALErrorString, Cstring, (Cint,), base))
    clear_errno!(fns)
    throw(LALError(func, Int(base), msg))
end
