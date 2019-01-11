module MFCC


using DSP, FFTW
include("Mel.jl")
include("Filterbanks.jl")
using MFCC.Mel
using MFCC.Filterbanks


"""
    mfcc(signal::AbstractVector, window::Function = ones, windowsize = 400, overlap = 160, sr = 44100, numcep = 26)

Compute the Mel-Frequency Cepstral Coefficients of the given signal and return
an Array where each column contains the coefficients of a given frame obtained
from the stft.

**Parameters:**

`signal`: The signal on which the coefficients will be computed.

`window`: The window function to be applied on the signal by the STFT.

`windowsize`: The size of the STFT window.

`overlap`: The overlap of two adjacent windows of the STFT in samples.

`sr`: Sample rate of the signal.

`numcep`: The amount of cepstral coefficients to be extracted from the signal.


See also: [`Filterbank`](@ref)
"""
function mfcc(signal::AbstractVector,
	window::Function = ones,
	windowsize::Int = 400,
	overlap::Int = 160,
	sr::Int = 44100,
	numcep::Int = 26
)

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
