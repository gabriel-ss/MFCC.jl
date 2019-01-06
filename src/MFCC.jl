module MFCC


using DSP
include("Mel.jl")
include("Filterbanks.jl")
using MFCC.Mel
using MFCC.Filterbanks


function mfcc(signal, window = ones, windowsize = 400, overlap = 160, sr = 44100, numcep = 26)

	#Each column is a frame of stft
	sigstft = stft(signal, windowsize, overlap, window=window)

	#If the signal is real then the absolute value of it's spectrum is symmetric
	eltype(signal)<:Real && (windowsize = windowsize ÷ 2 + 1; sr ÷ 2)

	#Evaluate the periodogram of the signal
	periodgrm = map(frame->(abs2.(frame)/windowsize), sigstft)

	#Create the filterbank
	windows, windowstarts = mffilters(numcep, windowsize, sr)
	melfilterbank = Filterbank(windowsize, windows[1:end], windowstarts[1:end])

	#Apply the filterbank and sum the result of each filter
	filteredsignal = mapslices(frame->melfilterbank*frame, periodgrm, dims = 1)

	#Take the Type-III DCT of the log for each coeficient
	return mapslices(frame->idct(log.(frame)), filteredsignal, dims=1)

end

export mfcc

end # module
