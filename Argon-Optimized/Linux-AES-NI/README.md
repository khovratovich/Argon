Argon Optimized for Linux with AES-NI support
=====

This implementation is optimized for x86/x64 architectures with AES-NI support. The AES instructions and 128-bit registers are used within intrinsics, some of which are Windows-specific.

Can be compiled with: (gcc-4.7.3)
	g++ -O3 -maes -mavx -std=c++11 -o Argon-Optimized ./Argon-Optimized.cpp -lpthread
