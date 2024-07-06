
AVFFT v1.6 - Optimized and very precise C++ templated SIMD vector forward and inverse real/complex FFT transform for iOS, OS X, and other platforms
----------------------------------------------------------------------------------------------------------------

NEW:  v1.6 now includes my useful utility class - const1.h, which is required for compiling SIMD AVFFT 1.6.
I fixed major bugs again in real_fft(): it is now computing correctly and with the most efficiency. I think this 
will be the last update for a while since I feel that I accomplished what I wanted to do with vectorization of FFT.
It is fast and the sound quality is outstanding prodced by this AV FFT when used in my audio project.

This is a fast and well tuned templated C++ real/complex implementation of FFT transform based on 
AV MPEG C version. It works well and compiles in the most recent version of XCode 15.4 and iOS 17.5. 

I think it is even a bit faster than PFFFT when used in simd vector mode, it sounds more precise 
to me than the PFFFT, but I suppose that could be subjective, though I have a very refined 
sound engineer ear.

The input/output is 100% compatible with PFFFT so you can just drop replace PFFFT if you like. 
I made Zita Convolver work with AVFFT in simd_float8 mode, which was a great accomplishment too,
and I will post it as soon as I double check everything. 

To do: Implement assembly ARM optimizations for double precision, and vectors.

enjoy.

Dimitry (not from Paris ;)
