AVFFT v1.21 - Optimized templatable SIMD vector capable forward and reverse real/complex fft transform.

This is a well tuned templated C++ real/complex implementation of FFT based on AV MPEG C version,
it works well and I am using it actively in my current multimedia project. I think it is even a bit 
faster than PFFFT when used in simd vector mode, it sounds more precise and sounds better than PFFFT
to me at least. I have a very refined hearing, so I used this in a IR convolver, and I could tell significant
difference from PFFFT.

I tested it on iOS and OS X.
