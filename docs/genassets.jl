using MFCC, MFCC.Mel, MFCC.Filterbanks, Plots



windows, windowstarts = mffilters(10, 22050, 22050)
filterbank = Filterbank(22050, windows, windowstarts)
melscale = mel.(range(1, length=22050, 22050))

# Plot the curve of the hertz scale by the mel scale
melhz = plot(melscale, xlabel="Frequency(Hz)", ylabel="Frequency(Mel)")
savefig(melhz, "docs/assets/mel_hz_mapping.svg")

# Plot the filterbank in the mel scale
melfb = plot(melscale, getfilters(filterbank), xlabel="Frequency(Mel)")
savefig(melfb, "docs/assets/mel_scale_filterbank.svg")

# Plot the filterbank in the hertz scale
hzfb = plot(getfilters(filterbank), xlabel="Frequency(Hz)")
savefig(hzfb, "docs/assets/hz_scale_filterbank.svg")
