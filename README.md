Argon
=====

These are implementations of the password hashing scheme Argon.

Argon was created by Alex Biryukov and Dmtry Khovratovich from University of Luxembourg.

Argon's main feature is tough resistance to tradeoff attacks. An adversary, who wants to compute the hash using less memory (for instance, to run a password cracker on fast but memory-unfriendly architectures like GPUs), is penalized significantly in computational costs. For example, using 64 MB instead of 128 MB gives a penalty factor of 160 to the entire hash computation.

Argon has its own webpage at https://www.cryptolux.org/index.php/Argon
