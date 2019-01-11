module Mel


"""
    mel(frequency::Number)

Transform the input to mel scale
"""
mel(frequency::Number)::Number = 1127*log(1 + frequency/700)


"""
    imel(mel::Number)

Transform the input back from mel scale
"""
imel(mel::Number)::Number = 700*(exp(mel/1127) - 1)


"""
    melspace(start::Number, length::Int, stop::Number)::Array

Return an array of mel spaced points
"""
melspace(start::Number, length::Int, stop::Number)::Array =
	imel.(range(mel(start), length=length, mel(stop)))


"""
    mffilters(noffilters::Int, filterlength::Int, sr::Real) -> (windows, starts)

Project a filterbank composed by triangular and even-spaced filters in the mel
scale.


The first element of the tuple is an array containing arrays, each one filled
with the non-zero part of a filter from a filterbank. The second element is an
array containing the start of the non-zero part of each filter. The return can
be used to create a [`Filterbank`](@ref) defined in the Filterbank module.


**Parameters:**

`noffilters`: The number of filters in the filterbank.

`filterlenght`: The length of the filters, should match the lenght of the lenght
of the signal to which the filterbank will be applied.

`sr`: The sampling rate to be used, should match the original sampling rate of
the sequence.


See also: [`Filterbank`](@ref)
"""
function mffilters(noffilters::Int, filterlength::Int, sr::Real)

	melpoints = ceil.(Int, melspace(1, noffilters + 2, sr)*filterlength/sr)

	return (map(m->	begin

		start, apex, stop = melpoints[m:m + 2]
		w = zeros(stop - start)

		for k in 1:(apex - start)
			w[k] = k/(apex - start)
		end

		for k in apex:(stop - 1)
			w[k - start + 1] = (stop - k)/(stop - apex)
		end

		return 2w/(stop - start)	#Normalized so that the triangle has area 1

	end, 1:noffilters),
	melpoints[1:end - 2])

end

export mel, imel, melspace, mffilters
end  # module Mel
