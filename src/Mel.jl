module Mel


mel(frequency::Number) = 1127*log(1 + frequency/700)

imel(mel) = 700*(exp(mel/1127) - 1)

melspace(start, length, stop) = imel.(range(mel(start), length=length, mel(stop)))

function mffilters(noffilters, filterlength, sr)

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
