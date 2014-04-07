Argon Optimized for Windows with AES-NI support
=====

This implementation is optimized for x86/x64 architectures with AES-NI support. The AES instructions and 128-bit registers are used within intrinsics, some of which are Windows-specific.

The implementation also supports multiple threads with the Boost library.
