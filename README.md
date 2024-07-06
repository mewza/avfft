AVFFT v1.6 - Optimized C++ templated SIMD-vector capable forward and reverse real/complex fft transform

NEW:  v1.6 now includes my useful const1.h which is required for compiling SIMD AVFFT 1.6. I fixed major bugs
again in real_fft() it is now computing it correctly and with the most efficiency. I think this 
will be the last update for a while since I feel that I accomplished what I wanted to do with vectorization of FFT.
It is working really good finally! it is fast and the sound quality is outstanding prodced by this AV FFT!

This is a fast and well tuned templated C++ real/complex implementation of FFT transform based on 
AV MPEG C version. It works well and I am using it actively in my current multimedia project. 

I think it is even a bit faster than PFFFT when used in simd vector mode, it sounds more precise 
to me than the PFFFT, but I suppose that could be subjective, though I have a very refined 
sound engineer ear.

The input/output is 100% compatible with PFFFT so you can just drop replace PFFFT if you like. 
I made Zita Convolver work with AVFFT in simd_float8 mode, which was a great accomplishment too,
and I will post it as soon as I double check everything. 

enjoy.

Dimitry (not from Paris ;)
