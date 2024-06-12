AVFFT v1.21 - Optimized C++ templated SIMD-vector capable forward and reverse real/complex fft transform.

This is a fast and well tuned templated C++ real/complex implementation of FFT transform based on 
AV MPEG C version. It works well and I am using it actively in my current multimedia project. 
I think it is even a bit faster than PFFFT when used in simd vector mode, it sounds more precise 
and sounds better than PFFFT to me at least. I have a very refined hearing, so I used this in a 
IR convolver (Zita Convolver which I prefer over others), and I could tell significant improvement
from PFFFT drop-replace by AVFFT. I am trying to make a vectorized version of Zita Convolver using
this AVFFT, I think it should signnificantly improve the performance of Zita, from the serialized
version which I used with PFFFT.

I tested it on iOS and OS X.
