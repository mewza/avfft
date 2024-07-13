
AVFFT v1.7 - Optimized and very precise C++ templated SIMD vector forward and inverse real/complex FFT transform for iOS, OS X, and other platforms
----------------------------------------------------------------------------------------------------------------

NEW:  v1.7 This is definitely IT. It sounds amaaazing! Yes, that good. Clean, clear, bassy, and detailed.
      Numerous improvements throughout, added sine/cosine tables and redone the real_fft(). Besides optimizations,
      this is a perfection so far, and yes, it is still compatible w/ PFFFT, in fact, I made sure that the two sounded
      as close as they can, in terms of data levels.

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

D
