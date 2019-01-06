module Filterbanks


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
