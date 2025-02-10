function follow_branch end
function plot_1D_solutions_branch end
function steady_state_sweep end

function plot_spaghetti end
function plot_eigenvalues end
function plot_rotframe_jacobian_response end
function plot_phase_diagram end
function plot_linear_response end

# ## Method error handling
# We also inject a method error handler, which
# prints a suggestion if the Proj extension is not loaded.
function _error_hinter(package, extension, func)
    return function _dump_error_hinter(io, exc, argtypes, kwargs)
        if isnothing(Base.get_extension(HarmonicBalance, extension)) && exc.f == func
            print(
                io,
                "\n\nThe `$(string(func))` method requires the $package.jl package to be explicitly loaded.\n",
            )
            print(io, "You can do this by simply typing ")
            printstyled(io, "using $package"; color=:cyan, bold=true)
            println(
                io, " in your REPL, \nor otherwise loading $package.jl via using or import."
            )
        else # this is a more general error
            nothing
        end
    end
end
