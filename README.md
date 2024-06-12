AVFFT v1.11 - Optimized SIMD vector capable real/complex fft transform (forward/backward)

This is a well tuned templated real/complex implementation of FFT as C++ class based on AV MPEG C version,
it works well and I am using it actively in my current multimedia project. I think it is even a bit 
faster than PFFFTD when used in simd vector mode, it sounds more precise and sounds better than PFFFT
to me at least. I have a very refined hearing, so I used this in a IR convolver, and I could tell significant
difference from PFFFT.
