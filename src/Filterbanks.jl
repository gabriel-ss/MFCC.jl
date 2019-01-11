module Filterbanks

"""
    Filterbank(length::Int, windows::Vector{Vector{Number}}, windowstarts::Vector{Int})

A representation of a filterbank to be applied to a signal with `length`
samples. The `windows` is a vector of vectors, each representing the nonzero
part of a filter and the `windowstarts` is a vector containing the start of
the nonzero portion of the filter.

The `*` operator is overloaded to apply the filterbank to a signal as if the
filterbank were a matrix and an multiplication was performed, but in a more
optimized way.


See also: [`getfilters`](@ref), [`applyfilters`](@ref)
"""
struct Filterbank
	length::Int
	windows::Vector{Vector{Number}}
	windowstarts::Vector{Int}
end


Base.:*(filterbank::Filterbank, signal::AbstractArray) = begin

	return map((window, winstart)->
	window'*signal[winstart:winstart + length(window) - 1]
	, filterbank.windows, filterbank.windowstarts)

end


"""
    getfilters(filterbank::Filterbank)

Return a vector containing vectors filled with the values of each filter of the
given filterbank, including the zero portion.

See also: [`applyfilters`](@ref)
"""
function getfilters(filterbank::Filterbank)

	return hcat(map((window, winstart)->
		[
			zeros(winstart)
			window
			zeros(filterbank.length - (winstart + length(window)))
		], filterbank.windows, filterbank.windowstarts)...)

end


function filtermesh(filterbank::Filterbank)

	filtermesh = zeros(filterbank.length)

	foreach((window, winstart)->
		filtermesh[winstart:winstart + length(window) - 1] += window
	, filterbank.windows, filterbank.windowstarts)

	return filtermesh

end


"""
    applyfilters(signal, filterbank::Filterbank)

Apply a filterbank to a signal and return the result.

See also: [`getfilters`](@ref)
"""
function applyfilters(signal, filterbank::Filterbank)

	return hcat(map((window, winstart)->
		[
			zeros(winstart)
			window.*signal[winstart:winstart + length(window) - 1]
			zeros(filterbank.length - (winstart + length(window)))
		], filterbank.windows, filterbank.windowstarts)...)

end


export Filterbank, getfilters, filtermesh, applyfilters
end  # module Filterbanks
