AVFFT v1.5 - Optimized C++ templated SIMD-vector capable forward and reverse real/complex fft transform

This is a fast and well tuned templated C++ real/complex implementation of FFT transform based on 
AV MPEG C version. It works well and I am using it actively in my current multimedia project. 

I think it is even a bit faster than PFFFT when used in simd vector mode, it sounds more precise 
to me at least, but I suppose that's subjective, though I do have a very refined sound engineer hearing.

The input/output is 100% compatible with PFFFT so you can just drop replace PFFFT if you like. 
I made Zita Convolver work in simd_float8 mode, which was another great accomplishment. 

So what's new in 1.5? and why a sudden version jump? I revamped the whole thing, cleaned it up, 
added my nifty cmplxT<T> class that I made, and added ARM64 asm optimizations but they will only 
be applied if you defined zfloat as a float and not double. I will adapt assembler code to work 
w/ doubles later... enjoy.

Dimitry (not from Paris ;)
