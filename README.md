
AVFFT v1.75 - Optimized and very precise C++ templated SIMD vector forward and inverse real/complex FFT transform for iOS, OS X, and other platforms
----------------------------------------------------------------------------------------------------------------

NEW 9/28/24: v1.75 I've moved N into template parameter, easier that way, but arrays used in the inner core
             are still allocated, that is due to assembler compatibility that indexes into class structure by offset.

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
DMT
